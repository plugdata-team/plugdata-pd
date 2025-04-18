/* Copyright (c) 1997-1999 Miller Puckette.
* For information on usage and redistribution, and for a DISCLAIMER OF ALL
* WARRANTIES, see the file, "LICENSE.txt," in this distribution.  */

/* sinusoidal oscillator and table lookup; see also tabosc4~ in d_array.c.
*/

#include "m_pd.h"
#include "math.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define BIGFLOAT 1.0e+19
#define UNITBIT32 1572864.  /* 3*2^19; bit 32 has place value 1 */

#include "m_private_utils.h"

#if BYTE_ORDER == LITTLE_ENDIAN
# define HIOFFSET 1
# define LOWOFFSET 0
#else
# define HIOFFSET 0    /* word offset to find MSB */
# define LOWOFFSET 1    /* word offset to find LSB */
#endif


union tabfudge
{
    double tf_d;
    int32_t tf_i[2];
};

/* the table size before 0.55 was hard-coded to 512 points.  Now we let it
default to 2048, but we can compile with COSTABSIZE defined to choose another
table size.  For the Pd application we also include a 512-point version for
back compatibility.  Since the code (on Intel arch at least) runs 10% faster
with COSTABSIE hard-coded, we have to make the choices at compile time.  If
the compatiblity version is needed, we end up compiling two version of the
same code.  This is done by including d_osc.h, twice if two table sizes must
be defined, to define the performance routines for cos~, osc~, and vcf~.  The
inclusions are done at the end of the file (after the instance structs are
defined) but are declared here. */

struct _cos;
struct _osc;
struct _vcfctl;

#ifndef COSTABLESIZE
#define COSTABLESIZE 2048
#define OLDTABSIZE 512  /* size of pre-0.55 compatibility table */
static t_int *cos_perform_old(t_int *w);
static t_int *osc_perform_old(t_int *w);
static t_int *sigvcf_perform_old(t_int *w);
#endif
static t_int *cos_perform(t_int *w);
static t_int *osc_perform(t_int *w);
static t_int *sigvcf_perform(t_int *w);
#define COSTABLENAME cos_newtable   /* keep cos_table back-compatibile */


/* -------------------------- phasor~ ------------------------------ */
static t_class *phasor_class;

#if 1   /* in the style of R. Hoeldrich (ICMC 1995 Banff) */

typedef struct _phasor
{
    t_object x_obj;
    double x_phase;
    t_float x_conv;
    t_float x_f;						// scalar frequency
} t_phasor;

static void *phasor_new(t_floatarg f)
{
    t_phasor *x = (t_phasor *)pd_new(phasor_class);
    x->x_f = f;
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_float, gensym("ft1"));
    x->x_phase = 0;
    x->x_conv = 0;
    outlet_new(&x->x_obj, gensym("signal"));
    return (x);
}

static t_int *phasor_perform(t_int *w)
{
    t_phasor *x = (t_phasor *)(w[1]);
    t_sample *in = (t_float *)(w[2]);
    t_sample *out = (t_float *)(w[3]);
    int n = (int)(w[4]);
    double dphase = x->x_phase + (double)UNITBIT32;
    union tabfudge tf;
    int normhipart;
    t_float conv = x->x_conv;

    tf.tf_d = UNITBIT32;
    normhipart = tf.tf_i[HIOFFSET];
    tf.tf_d = dphase;

    while (n--)
    {
        tf.tf_i[HIOFFSET] = normhipart;
        dphase += *in++ * conv;
        *out++ = tf.tf_d - UNITBIT32;
        tf.tf_d = dphase;
    }
    tf.tf_i[HIOFFSET] = normhipart;
    x->x_phase = tf.tf_d - UNITBIT32;
    return (w+5);
}

static void phasor_dsp(t_phasor *x, t_signal **sp)
{
    x->x_conv = 1./sp[0]->s_sr;
    dsp_add(phasor_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, (t_int)sp[0]->s_n);
}

static void phasor_ft1(t_phasor *x, t_float f)
{
    x->x_phase = (double)f;
}

static void phasor_setup(void)
{
    phasor_class = class_new(gensym("phasor~"), (t_newmethod)phasor_new, 0,
        sizeof(t_phasor), 0, A_DEFFLOAT, 0);
    CLASS_MAINSIGNALIN(phasor_class, t_phasor, x_f);
    class_addmethod(phasor_class, (t_method)phasor_dsp,
        gensym("dsp"), A_CANT, 0);
    class_addmethod(phasor_class, (t_method)phasor_ft1,
        gensym("ft1"), A_FLOAT, 0);
}

#endif  /* Hoeldrich version */

/* ------------------------ cos~ ----------------------------- */

float *cos_table;
static float *cos_newtable;

static t_class *cos_class;

typedef struct _cos
{
    t_object x_obj;
    t_float x_f;			// scalar frequency
} t_cos;

static void *cos_new(t_floatarg f)
{
    t_cos *x = (t_cos *)pd_new(cos_class);
    outlet_new(&x->x_obj, gensym("signal"));
    x->x_f = f;
    return (x);
}

static void cos_dsp(t_cos *x, t_signal **sp)
{
    signal_setmultiout(&sp[1], sp[0]->s_nchans);
#ifdef OLDTABSIZE
    if (pd_compatibilitylevel < 55)
        dsp_add(cos_perform_old, 3, sp[0]->s_vec, sp[1]->s_vec,
            (t_int)(sp[0]->s_length * sp[0]->s_nchans));
    else
#endif
    dsp_add(cos_perform, 3, sp[0]->s_vec, sp[1]->s_vec,
        (t_int)(sp[0]->s_length * sp[0]->s_nchans));
}

static void cos_maketable(void)
{
    union tabfudge tf;
    int i;
       /* here we check at startup whether the byte alignment
            is as we declared it.  If not, the code has to be
            recompiled the other way. */
    tf.tf_d = UNITBIT32 + 0.5;
    if ((unsigned)tf.tf_i[LOWOFFSET] != 0x80000000)
        bug("cos~: unexpected machine alignment");
    if (!cos_newtable)
    {
        cos_newtable = (float *)getbytes(sizeof(float) * (COSTABLESIZE+1));
        for (i = 0; i < COSTABLESIZE + 1; i++)
            cos_newtable[i] = cos(2*M_PI*i / (double)COSTABLESIZE);
            /* fill in true 1s and 0s for 1/4-cycle points.  These should
            be possible to address exactly in cos~ object.  The rest are
            irrational anyway and hence will never be exact so we leave them
            as computed by the lobrary cosine function. */
        cos_newtable[0] = cos_newtable[COSTABLESIZE] = 1;
        cos_newtable[COSTABLESIZE/4] = cos_newtable[3*COSTABLESIZE/4] = 0;
        cos_newtable[COSTABLESIZE/2] = -1;
    }
#ifdef OLDTABSIZE
    if (!cos_table)
    {
        float *fp, phase, phsinc = (2. * 3.14159) / OLDTABSIZE;
        cos_table = (float *)getbytes(sizeof(float) * (OLDTABSIZE+1));
        for (i = OLDTABSIZE + 1, fp = cos_table, phase = 0; i--;
            fp++, phase += phsinc)
                *fp = cos(phase);
    }
#endif
}

static void cos_cleanup(t_class *c)
{
#ifdef OLDTABSIZE
    if (cos_table)
        freebytes(cos_table, sizeof(float) * (OLDTABSIZE+1));
    cos_table = 0;
#endif
    if (cos_newtable)
        freebytes(cos_newtable, sizeof(float) * (COSTABLESIZE+1));
    cos_newtable = 0;
}

static void cos_setup(void)
{
    cos_class = class_new(gensym("cos~"), (t_newmethod)cos_new, 0,
        sizeof(t_cos), CLASS_MULTICHANNEL, A_DEFFLOAT, 0);
    class_setfreefn(cos_class, cos_cleanup);
    CLASS_MAINSIGNALIN(cos_class, t_cos, x_f);
    class_addmethod(cos_class, (t_method)cos_dsp, gensym("dsp"), A_CANT, 0);
    cos_maketable();
}

/* ------------------------ osc~ ----------------------------- */

static t_class *osc_class, *scalarosc_class;

typedef struct _osc
{
    t_object x_obj;
    double x_phase;
    t_float x_conv;
    t_float x_f;						// scalar frequency
} t_osc;

static void *osc_new(t_floatarg f)
{
    t_osc *x = (t_osc *)pd_new(osc_class);
    x->x_f = f;
    outlet_new(&x->x_obj, gensym("signal"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_float, gensym("ft1"));
    x->x_phase = 0;
    x->x_conv = 0;
    return (x);
}

static void osc_dsp(t_osc *x, t_signal **sp)
{
    x->x_conv = COSTABLESIZE/sp[0]->s_sr;
#ifdef OLDTABSIZE
    if (pd_compatibilitylevel < 55)
    {
        x->x_conv = OLDTABSIZE/sp[0]->s_sr;
        dsp_add(osc_perform_old, 4, x, sp[0]->s_vec, sp[1]->s_vec,
            (t_int)sp[0]->s_n);
    }
    else
#endif
        dsp_add(osc_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec,
            (t_int)sp[0]->s_n);
}

static void osc_ft1(t_osc *x, t_float f)
{
#ifdef OLDTABSIZE
    if (pd_compatibilitylevel < 55)
        x->x_phase = OLDTABSIZE * f;
    else
#endif
        x->x_phase = COSTABLESIZE * f;
}

static void osc_setup(void)
{
    osc_class = class_new(gensym("osc~"), (t_newmethod)osc_new, 0,
        sizeof(t_osc), 0, A_DEFFLOAT, 0);
    CLASS_MAINSIGNALIN(osc_class, t_osc, x_f);
    class_addmethod(osc_class, (t_method)osc_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(osc_class, (t_method)osc_ft1, gensym("ft1"), A_FLOAT, 0);

    cos_maketable();
}

/* ---- vcf~ - resonant filter with audio-rate center frequency input ----- */

typedef struct _vcfctl
{
    t_float c_re;
    t_float c_im;
    t_float c_q;
    t_float c_isr;
} t_vcfctl;

typedef struct sigvcf
{
    t_object x_obj;
    t_vcfctl x_cspace;
    t_float x_f;
} t_sigvcf;

t_class *sigvcf_class;

static void *sigvcf_new(t_floatarg q)
{
    t_sigvcf *x = (t_sigvcf *)pd_new(sigvcf_class);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("ft1"));
    outlet_new(&x->x_obj, gensym("signal"));
    outlet_new(&x->x_obj, gensym("signal"));
    x->x_cspace.c_re = 0;
    x->x_cspace.c_im = 0;
    x->x_cspace.c_q = q;
    x->x_cspace.c_isr = 0;
    x->x_f = 0;
    return (x);
}

static void sigvcf_ft1(t_sigvcf *x, t_float f)
{
    if(f < 0.) f = 0.;
    if(f > BIGFLOAT) f = BIGFLOAT;
    x->x_cspace.c_q = f;
}

static void sigvcf_dsp(t_sigvcf *x, t_signal **sp)
{
    x->x_cspace.c_isr = 6.28318f/sp[0]->s_sr;
#ifdef OLDTABSIZE
    if (pd_compatibilitylevel < 55)
        dsp_add(sigvcf_perform_old, 6,
            sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec,
                &x->x_cspace, (t_int)sp[0]->s_n);
    else
#endif
        dsp_add(sigvcf_perform, 6,
            sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec,
                &x->x_cspace, (t_int)sp[0]->s_n);
}

static void sigvcf_clear(t_sigvcf *x)
{
    x->x_cspace.c_re = 0;
    x->x_cspace.c_im = 0;
}

static
void sigvcf_setup(void)
{
    sigvcf_class = class_new(gensym("vcf~"), (t_newmethod)sigvcf_new, 0,
        sizeof(t_sigvcf), 0, A_DEFFLOAT, 0);
    CLASS_MAINSIGNALIN(sigvcf_class, t_sigvcf, x_f);
    class_addmethod(sigvcf_class, (t_method)sigvcf_clear,
        gensym("clear"), 0);
    class_addmethod(sigvcf_class, (t_method)sigvcf_dsp,
        gensym("dsp"), A_CANT, 0);
    class_addmethod(sigvcf_class, (t_method)sigvcf_ft1,
        gensym("ft1"), A_FLOAT, 0);
}

/* -------------------------- noise~ ------------------------------ */
static t_class *noise_class;

typedef struct _noise
{
    t_object x_obj;
    int x_val;
} t_noise;

static void *noise_new(void)
{
    t_noise *x = (t_noise *)pd_new(noise_class);
        /* seed each instance differently.  Once in a blue moon two threads
        could grab the same seed value.  We can live with that. */
    static int init = 307;
    x->x_val = (init *= 1319);
    outlet_new(&x->x_obj, gensym("signal"));
    return (x);
}

static t_int *noise_perform(t_int *w)
{
    t_sample *out = (t_sample *)(w[1]);
    int *vp = (int *)(w[2]);
    int n = (int)(w[3]);
    int val = *vp;
    while (n--)
    {
        *out++ = ((t_sample)((val & 0x7fffffff) - 0x40000000)) *
            (t_sample)(1.0 / 0x40000000);
        val = val * 435898247 + 382842987;
    }
    *vp = val;
    return (w+4);
}

static void noise_dsp(t_noise *x, t_signal **sp)
{
    dsp_add(noise_perform, 3, sp[0]->s_vec, &x->x_val, (t_int)sp[0]->s_n);
}

static void noise_float(t_noise *x, t_float f)
{
    /* set the seed */
    x->x_val = (int)f;
}

static void noise_setup(void)
{
    noise_class = class_new(gensym("noise~"), (t_newmethod)noise_new, 0,
        sizeof(t_noise), 0, 0);
    class_addmethod(noise_class, (t_method)noise_dsp,
        gensym("dsp"), A_CANT, 0);
    class_addmethod(noise_class, (t_method)noise_float,
        gensym("seed"), A_FLOAT, 0);
}


/******************** tabosc4~ ***********************/

static t_class *tabosc4_tilde_class;

typedef struct _tabosc4_tilde
{
    t_object x_obj;
    t_float x_fnpoints;
    t_float x_finvnpoints;
    t_word *x_vec;
    t_symbol *x_arrayname;
    t_float x_f;
    double x_phase;
    t_float x_conv;
} t_tabosc4_tilde;

static void *tabosc4_tilde_new(t_symbol *s)
{
    t_tabosc4_tilde *x = (t_tabosc4_tilde *)pd_new(tabosc4_tilde_class);
    x->x_arrayname = s;
    x->x_vec = 0;
    x->x_fnpoints = 512.;
    x->x_finvnpoints = (1./512.);
    outlet_new(&x->x_obj, gensym("signal"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_float, gensym("ft1"));
    x->x_f = 0;
    return (x);
}

static t_int *tabosc4_tilde_perform(t_int *w)
{
    t_tabosc4_tilde *x = (t_tabosc4_tilde *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);
    t_sample *out = (t_sample *)(w[3]);
    int n = (int)(w[4]);
    int normhipart;
    union tabfudge tf;
    t_float fnpoints = x->x_fnpoints;
    int mask = fnpoints - 1;
    t_float conv = fnpoints * x->x_conv;
    t_word *tab = x->x_vec, *addr;
    double dphase = fnpoints * x->x_phase + UNITBIT32;

    if (!tab) goto zero;
    tf.tf_d = UNITBIT32;
    normhipart = tf.tf_i[HIOFFSET];

#if 1
    while (n--)
    {
        t_sample frac,  a,  b,  c,  d, cminusb;
        tf.tf_d = dphase;
        dphase += *in++ * conv;
        addr = tab + (tf.tf_i[HIOFFSET] & mask);
        tf.tf_i[HIOFFSET] = normhipart;
        frac = tf.tf_d - UNITBIT32;
        a = addr[0].w_float;
        b = addr[1].w_float;
        c = addr[2].w_float;
        d = addr[3].w_float;
        cminusb = c-b;
        *out++ = b + frac * (
            cminusb - 0.1666667f * (1.-frac) * (
                (d - a - 3.0f * cminusb) * frac + (d + 2.0f*a - 3.0f*b)
            )
        );
    }
#endif

    tf.tf_d = UNITBIT32 * fnpoints;
    normhipart = tf.tf_i[HIOFFSET];
    tf.tf_d = dphase + (UNITBIT32 * fnpoints - UNITBIT32);
    tf.tf_i[HIOFFSET] = normhipart;
    x->x_phase = (tf.tf_d - UNITBIT32 * fnpoints)  * x->x_finvnpoints;
    return (w+5);
 zero:
    while (n--) *out++ = 0;

    return (w+5);
}

static void tabosc4_tilde_set(t_tabosc4_tilde *x, t_symbol *s)
{
    t_garray *a;
    int npoints, pointsinarray;

    x->x_arrayname = s;
    if (!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
    {
        if (*s->s_name)
            pd_error(x, "tabosc4~: %s: no such array", x->x_arrayname->s_name);
        x->x_vec = 0;
    }
    else if (!garray_getfloatwords(a, &pointsinarray, &x->x_vec))
    {
        pd_error(x, "%s: bad template for tabosc4~", x->x_arrayname->s_name);
        x->x_vec = 0;
    }
    else if ((npoints = pointsinarray - 3) != (1 << ilog2(pointsinarray - 3)))
    {
        pd_error(x, "%s: number of points (%d) not a power of 2 plus three",
            x->x_arrayname->s_name, pointsinarray);
        x->x_vec = 0;
    }
    else
    {
        x->x_fnpoints = npoints;
        x->x_finvnpoints = 1./npoints;
        garray_usedindsp(a);
    }
}

static void tabosc4_tilde_ft1(t_tabosc4_tilde *x, t_float f)
{
    x->x_phase = f;
}

static void tabosc4_tilde_dsp(t_tabosc4_tilde *x, t_signal **sp)
{
    x->x_conv = 1. / sp[0]->s_sr;
    tabosc4_tilde_set(x, x->x_arrayname);

    dsp_add(tabosc4_tilde_perform, 4, x,
        sp[0]->s_vec, sp[1]->s_vec, (t_int)sp[0]->s_n);
}

static void tabosc4_tilde_setup(void)
{
    tabosc4_tilde_class = class_new(gensym("tabosc4~"),
        (t_newmethod)tabosc4_tilde_new, 0,
        sizeof(t_tabosc4_tilde), 0, A_DEFSYM, 0);
    CLASS_MAINSIGNALIN(tabosc4_tilde_class, t_tabosc4_tilde, x_f);
    class_addmethod(tabosc4_tilde_class, (t_method)tabosc4_tilde_dsp,
        gensym("dsp"), A_CANT, 0);
    class_addmethod(tabosc4_tilde_class, (t_method)tabosc4_tilde_set,
        gensym("set"), A_SYMBOL, 0);
    class_addmethod(tabosc4_tilde_class, (t_method)tabosc4_tilde_ft1,
        gensym("ft1"), A_FLOAT, 0);
}

#define TABSIZE COSTABSIZE
#define COSPERF cos_perform
#define OSCPERF osc_perform
#define SIGVCFPERF sigvcf_perform

#include "d_osc.h"    /* include normal perf routines */

#ifdef OLDTABSIZE
#undef COSTABLESIZE
#define COSTABLESIZE OLDTABSIZE
#undef COSTABLENAME
#define COSTABLENAME cos_table
#undef COSPERF
#define COSPERF cos_perform_old
#undef OSCPERF
#define OSCPERF osc_perform_old
#undef SIGVCFPERF
#define SIGVCFPERF sigvcf_perform_old
#include "d_osc.h"    /* include 512-point compatibility perf routines */
#endif

/* ----------------------- global setup routine ---------------- */
void d_osc_setup(void)
{
    phasor_setup();
    cos_setup();
    osc_setup();
    sigvcf_setup();
    noise_setup();
    tabosc4_tilde_setup();
}
