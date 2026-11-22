/* Copyright (c) 1997- Miller Puckette and others.
* For information on usage and redistribution, and for a DISCLAIMER OF ALL
* WARRANTIES, see the file, "LICENSE.txt," in this distribution.  */

/* --------- Pd interface to FFTW library; imitate Mayer API ---------- */

/* changes and additions for FFTW3 by Thomas Grill                      */

#include "m_pd.h"
#include "m_imp.h"
#include <string.h>
#include <fftw3.h>
#if PDTHREADS && defined(PDINSTANCE)
#include <pthread.h>
#endif

int ilog2(int n);

#define MINFFT 0
#define MAXFFT 30

/* from the FFTW website:
 #include <fftw3.h>
     ...
     {
         fftw_complex *in, *out;
     fftw_plan p;
     ...
         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
         p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
     ...
         fftw_execute(p);
     ...
     fftw_destroy_plan(p);
         fftw_free(in); fftw_free(out);
     }

FFTW_FORWARD or FFTW_BACKWARD, and indicates the direction of the transform you
are interested in. Alternatively, you can use the sign of the exponent in the
transform, -1 or +1, which corresponds to FFTW_FORWARD or FFTW_BACKWARD
respectively. The flags argument is either FFTW_MEASURE

*/

/* each Pd instance has its own plans and buffers, so that instances in
different threads can compute FFTs at the same time.  They are stored as
per-instance data of a class that exists just for this. */

#define NFFT (MAXFFT+1 - MINFFT)

    /* plugdata: FFTW_ESTIMATE avoids audio stalls when planning on demand.
    Real FFT planning took <1 ms with ESTIMATE, versus ~70 ms at 4096
    points and ~270 ms at 65536 with MEASURE. */
#define FFTW_PLANNER_FLAGS FFTW_ESTIMATE

    /* plugdata: Use FFTW_MEASURE for small FFTs during setup. Measuring
    costs ~30 ms once per process and avoids up to 1.6x longer execution
    times with FFTW_ESTIMATE. */
#define CFFT_MEASURE_MAXN 32
#define RFFT_HALFCOMPLEX_MAXN 64

typedef struct {
    fftwf_plan plan;
    fftwf_complex *in,*out;
} cfftw_info;

    /* plugdata: Larger real FFTs use r2c/c2r for SIMD speedups
    of 2-3x over halfcomplex plans. */
typedef struct {
    fftwf_plan plan;
    float *real;            /* n samples */
    float *spec;            /* spectrum: n halfcomplex values (r2r), or
                               n/2+1 interleaved complex bins (r2c/c2r) */
    int halfcomplex;
} rfftw_info;

typedef struct {
    cfftw_info c_fwd[NFFT], c_bwd[NFFT];
    rfftw_info r_fwd[NFFT], r_bwd[NFFT];
} t_fftwplans;

static t_class *fftwplans_class;

    /* of FFTW's routines only fftwf_execute() is thread-safe */
#if PDTHREADS && defined(PDINSTANCE)
static pthread_mutex_t planner_mutex = PTHREAD_MUTEX_INITIALIZER;
#define planner_lock() pthread_mutex_lock(&planner_mutex)
#define planner_unlock() pthread_mutex_unlock(&planner_mutex)
#else
#define planner_lock()
#define planner_unlock()
#endif

static void fftwplans_free(void *z);

static t_fftwplans *fftwplans_get(void)
{
    t_fftwplans *x = (t_fftwplans *)class_getinstancedata(fftwplans_class);
    if (!x)
    {
        x = (t_fftwplans *)getbytes(sizeof(*x));
        class_setinstancedata(fftwplans_class, x, fftwplans_free);
    }
    return (x);
}

/* complex stuff */

    /* premeasure: called from fftw_instance_setup(), not the audio thread */
static cfftw_info *cfftw_doplan(int n, int fwd, int premeasure)
{
    cfftw_info *info;
    t_fftwplans *x;
    int logn = ilog2(n);
    if (logn < MINFFT || logn > MAXFFT)
        return (0);
    x = fftwplans_get();
    info = (fwd?x->c_fwd:x->c_bwd)+(logn-MINFFT);
    if (!info->plan)
    {
        planner_lock();
        info->in =
            (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n);
        info->out =
            (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n);
        info->plan = fftwf_plan_dft_1d(n, info->in, info->out,
            fwd?FFTW_FORWARD:FFTW_BACKWARD,
            (premeasure && n <= CFFT_MEASURE_MAXN) ?
                FFTW_MEASURE : FFTW_PLANNER_FLAGS);
        planner_unlock();
    }
    return info;
}

static cfftw_info *cfftw_getplan(int n, int fwd)
{
    return cfftw_doplan(n, fwd, 0);
}

static void cfftw_term(t_fftwplans *x)
{
    int i, j;
    cfftw_info *cinfo[2];

    for (i = 0; i < NFFT; i++)
    {
      cinfo[0] = &x->c_fwd[i];
      cinfo[1] = &x->c_bwd[i];

      for (j = 0; j < 2; j++)
      {
        if (cinfo[j]->plan)
        {
          fftwf_destroy_plan(cinfo[j]->plan);
          fftwf_free(cinfo[j]->in);
          fftwf_free(cinfo[j]->out);
          cinfo[j]->plan = 0;
          cinfo[j]->in = 0;
          cinfo[j]->out = 0;
        }
      }
    }
}


/* real stuff */

static rfftw_info *rfftw_doplan(int n, int fwd, int premeasure)
{
    rfftw_info *info;
    t_fftwplans *x;
    unsigned hcflags = premeasure ? FFTW_MEASURE : FFTW_PLANNER_FLAGS;
    int logn = ilog2(n);
    if (logn < MINFFT || logn > MAXFFT)
        return (0);
    x = fftwplans_get();
    info = (fwd?x->r_fwd:x->r_bwd)+(logn-MINFFT);
    if (!info->plan)
    {
        planner_lock();
        info->real = (float*) fftwf_malloc(sizeof(float) * n);
            /* n+2 floats hold either layout of the spectrum */
        info->spec = (float*) fftwf_malloc(sizeof(float) * (n + 2));
        info->halfcomplex = (n <= RFFT_HALFCOMPLEX_MAXN);
        if (info->halfcomplex)
            info->plan = fwd ?
                fftwf_plan_r2r_1d(n, info->real, info->spec, FFTW_R2HC,
                    hcflags) :
                fftwf_plan_r2r_1d(n, info->spec, info->real, FFTW_HC2R,
                    hcflags);
        else
            info->plan = fwd ?
                fftwf_plan_dft_r2c_1d(n, info->real,
                    (fftwf_complex *)info->spec, FFTW_PLANNER_FLAGS) :
                fftwf_plan_dft_c2r_1d(n, (fftwf_complex *)info->spec,
                    info->real, FFTW_PLANNER_FLAGS);
        planner_unlock();
    }
    return info;
}

static rfftw_info *rfftw_getplan(int n, int fwd)
{
    return rfftw_doplan(n, fwd, 0);
}

static void rfftw_term(t_fftwplans *x)
{
    int i, j;
    rfftw_info *rinfo[2];

    for (i = 0; i < NFFT; i++)
    {
      rinfo[0] = &x->r_fwd[i];
      rinfo[1] = &x->r_bwd[i];

      for (j = 0; j < 2; j++)
      {
        if (rinfo[j]->plan)
        {
          fftwf_destroy_plan(rinfo[j]->plan);
          fftwf_free(rinfo[j]->real);
          fftwf_free(rinfo[j]->spec);
          rinfo[j]->plan = 0;
          rinfo[j]->real = 0;
          rinfo[j]->spec = 0;
        }
      }
    }
}

    /* called when a Pd instance is freed, or by class_free() below */
static void fftwplans_free(void *z)
{
    t_fftwplans *x = (t_fftwplans *)z;
    planner_lock();
    cfftw_term(x);
    rfftw_term(x);
    planner_unlock();
    freebytes(x, sizeof(*x));
}

    /* plugdata: Preplan common sizes at Pd instance creation to avoid
    planning them on the audio thread. */
void fftw_instance_setup(void)
{
    int i;
    sys_lock();
    for (i = 0; i < 12; i++)
    {
        cfftw_doplan(1 << i, 1, 1);
        cfftw_doplan(1 << i, 0, 1);
        rfftw_doplan(1 << i, 1, 1);
        rfftw_doplan(1 << i, 0, 1);
    }
    sys_unlock();
}

static int mayer_refcount = 0;

void mayer_init(void)
{
    if (mayer_refcount++ == 0)
        fftwplans_class = class_new(gensym("fftwplans"), 0, 0,
            sizeof(t_pd), CLASS_PD, 0);
}

void mayer_term(void)
{
    if (--mayer_refcount == 0)
    {
        class_free(fftwplans_class);
        fftwplans_class = 0;
    }
}


void mayer_fht(t_sample *fz, int n)
{
    post("FHT: not yet implemented");
}

static void mayer_do_cfft(int n, t_sample *fz1, t_sample *fz2, int fwd)
{
    int i;
    float *fz;
    cfftw_info *p = cfftw_getplan(n, fwd);
    if (!p)
        return;

    for (i = 0, fz = (float *)p->in; i < n; i++)
        fz[i*2] = fz1[i], fz[i*2+1] = fz2[i];

    fftwf_execute(p->plan);

    for (i = 0, fz = (float *)p->out; i < n; i++)
        fz1[i] = fz[i*2], fz2[i] = fz[i*2+1];
}

void mayer_fft(int n, t_sample *fz1, t_sample *fz2)
{
    mayer_do_cfft(n, fz1, fz2, 1);
}

void mayer_ifft(int n, t_sample *fz1, t_sample *fz2)
{
    mayer_do_cfft(n, fz1, fz2, 0);
}

/*
    in the following the sign flips are done to
    be compatible with the mayer_fft implementation,
    but it's probably the mayer_fft that should be corrected...
*/

    /* copy between Pd's sample vectors and FFTW's single-precision buffers.
    memcpy is noticeably faster than an element loop for small FFTs. */
static void samples_to_float(float *dst, const t_sample *src, int n)
{
#if PD_FLOATSIZE == 32
    memcpy(dst, src, sizeof(float) * n);
#else
    int i;
    for (i = 0; i < n; i++)
        dst[i] = src[i];
#endif
}

static void float_to_samples(t_sample *dst, const float *src, int n)
{
#if PD_FLOATSIZE == 32
    memcpy(dst, src, sizeof(float) * n);
#else
    int i;
    for (i = 0; i < n; i++)
        dst[i] = src[i];
#endif
}

    /* Pd's layout: the real parts of bins 0..n/2 in fz[0..n/2], then the
    negated imaginary part of bin n-i in fz[i] for i > n/2. This is FFTW's
    R2HC halfcomplex order with the second half negated. */
void mayer_realfft(int n, t_sample *fz)
{
    int i, n2 = n/2;
    rfftw_info *p = rfftw_getplan(n, 1);
    if (!p)
        return;

    samples_to_float(p->real, fz, n);
    fftwf_execute(p->plan);
    if (p->halfcomplex)
    {
        const float *s = p->spec;
        for (i = 0; i <= n2; i++)
            fz[i] = s[i];
        for (; i < n; i++)
            fz[i] = -s[i];
    }
    else
    {
        const fftwf_complex *c = (const fftwf_complex *)p->spec;
        for (i = 0; i <= n2; i++)
            fz[i] = c[i][0];
        for (i = n2 + 1; i < n; i++)
            fz[i] = -c[n - i][1];
    }
}

void mayer_realifft(int n, t_sample *fz)
{
    int i, n2 = n/2;
    rfftw_info *p = rfftw_getplan(n, 0);
    if (!p)
        return;

    if (p->halfcomplex)
    {
        float *s = p->spec;
        for (i = 0; i <= n2; i++)
            s[i] = fz[i];
        for (; i < n; i++)
            s[i] = -fz[i];
    }
    else
    {
        fftwf_complex *c = (fftwf_complex *)p->spec;
            /* the DC bin, and for even n the Nyquist bin, have no imaginary
            part */
        c[0][0] = fz[0];
        c[0][1] = 0;
        for (i = 1; i <= n2; i++)
        {
            c[i][0] = fz[i];
            c[i][1] = (n - i > i) ? -fz[n - i] : 0;
        }
    }
    fftwf_execute(p->plan);
    float_to_samples(fz, p->real, n);
}

    /* ancient ISPW-like version, used in fiddle~ and perhaps other externs
    here and there. */
void pd_fft(t_float *buf, int npoints, int inverse)
{
    cfftw_info *p = cfftw_getplan(npoints, !inverse);
    int i;
    float *fz;
    for (i = 0, fz = (float *)(p->in); i < 2 * npoints; i++)
        *fz++ = buf[i];
    fftwf_execute(p->plan);
    for (i = 0, fz = (float *)(p->out); i < 2 * npoints; i++)
        buf[i] = *fz++;
}

