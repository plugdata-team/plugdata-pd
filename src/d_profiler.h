/* Copyright (c) 2026 Timothy Schoen
* For information on usage and redistribution, and for a DISCLAIMER OF ALL
* WARRANTIES, see the file, "LICENSE.txt," in this distribution.  */

/* plugdata: DSP profiler, times the perform routines of each DSP object */

#ifndef __d_profiler_h_
#define __d_profiler_h_

#include "m_pd.h"

#if defined(_LANGUAGE_C_PLUS_PLUS) || defined(__cplusplus)
extern "C" {
#endif

typedef struct _dspprofilerentry
{
    t_object *e_object;     /* 0 for time spent between objects */
    int e_parent;           /* entry of the subpatch, abstraction or clone it's in, -1 at toplevel */
    double e_time;          /* seconds spent in its perform routines since the last read */
} t_dspprofilerentry;

    /* enabling it rebuilds the DSP chain */
EXTERN void dsp_profiler_setenabled(int enabled);
    /* number of entries the next read returns, 0 while DSP is off */
EXTERN int dsp_profiler_getsize(void);
    /* fills in at most n entries and starts counting from 0 again. audiotime is the duration of the
    audio computed since the last read, generation changes whenever the DSP chain is rebuilt */
EXTERN int dsp_profiler_read(t_dspprofilerentry *entries, int n, double *audiotime, unsigned int *generation);

#if defined(_LANGUAGE_C_PLUS_PLUS) || defined(__cplusplus)
}
#endif

/* The rest is only compiled into d_ugen.c, which keeps the profiler in THIS */
#ifdef DSP_PROFILER_IMPLEMENTATION

#include <string.h>
#ifdef _WIN32
#include <windows.h>
#elif defined(__APPLE__)
#include <mach/mach_time.h>
#else
#include <time.h>
#endif

/* While profiling, dsp_tick() times every perform routine. Each routine is
counted for the innermost object whose "dsp" method added it, so a subpatch,
abstraction or clone only gets what its objects don't add themselves, like
reblocking and summing signals inside it. Slot 0 gets the same outside of any
subpatch. Everything here runs with Pd's lock held. */

typedef struct _profilerslot
{
    t_object *s_object;
    int s_parent;           /* slot of the object it's in, -1 at toplevel */
    int s_onset;            /* where its perform routines start in the chain */
    uint64_t s_time;        /* clock ticks since the last read */
} t_profilerslot;

typedef struct _dspprofiler
{
    int *p_owner;           /* slot of each element of the DSP chain */
    int p_ownersize;
    int p_ownercapacity;
    t_profilerslot *p_slots;
    int p_nslots;
    int p_slotcapacity;
    uint64_t p_nticks;      /* DSP ticks since the last read */
    unsigned int p_generation;
} t_dspprofiler;

static uint64_t dsp_profiler_now(void)
{
#ifdef _WIN32
    LARGE_INTEGER now;
    QueryPerformanceCounter(&now);
    return ((uint64_t)now.QuadPart);
#elif defined(__APPLE__)
    return (mach_absolute_time());
#else
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return ((uint64_t)now.tv_sec * 1000000000u + (uint64_t)now.tv_nsec);
#endif
}

static double dsp_profiler_toseconds(uint64_t ticks)
{
    static double scale;
    if (!scale)
    {
#ifdef _WIN32
        LARGE_INTEGER frequency;
        QueryPerformanceFrequency(&frequency);
        scale = 1. / (double)frequency.QuadPart;
#elif defined(__APPLE__)
        mach_timebase_info_data_t timebase;
        mach_timebase_info(&timebase);
        scale = 1e-9 * (double)timebase.numer / (double)timebase.denom;
#else
        scale = 1e-9;
#endif
    }
    return (ticks * scale);
}

    /* forget the slots of the old DSP chain */
static void dsp_profiler_clear(void)
{
    t_dspprofiler *p = THIS->u_profiler;
    p->p_ownersize = 0;
    p->p_nslots = 1;
    p->p_slots[0].s_object = 0;
    p->p_slots[0].s_parent = -1;
    p->p_slots[0].s_onset = 0;
    p->p_slots[0].s_time = 0;
    p->p_nticks = 0;
    p->p_generation++;
}

static void dsp_profiler_free(void)
{
    t_dspprofiler *p = THIS->u_profiler;
    if (!p)
        return;
    freebytes(p->p_owner, p->p_ownercapacity * sizeof(*p->p_owner));
    freebytes(p->p_slots, p->p_slotcapacity * sizeof(*p->p_slots));
    freebytes(p, sizeof(*p));
    THIS->u_profiler = 0;
}

    /* called after an object's "dsp" method added the chain elements from
    onset to end */
static void dsp_profiler_record(t_object *obj, int onset, int end)
{
    t_dspprofiler *p = THIS->u_profiler;
    int i, slot;
    if (end <= onset)
        return;

    if (end > p->p_ownercapacity)
    {
        int capacity = (end > 2 * p->p_ownercapacity ?
            end : 2 * p->p_ownercapacity);
        p->p_owner = (int *)resizebytes(p->p_owner,
            p->p_ownercapacity * sizeof(*p->p_owner),
                capacity * sizeof(*p->p_owner));
        p->p_ownercapacity = capacity;
    }
    if (end > p->p_ownersize)
    {
        memset(p->p_owner + p->p_ownersize, 0,
            (end - p->p_ownersize) * sizeof(*p->p_owner));
        p->p_ownersize = end;
    }
    if (p->p_nslots == p->p_slotcapacity)
    {
        p->p_slots = (t_profilerslot *)resizebytes(p->p_slots,
            p->p_slotcapacity * sizeof(*p->p_slots),
                2 * p->p_slotcapacity * sizeof(*p->p_slots));
        p->p_slotcapacity *= 2;
    }

    slot = p->p_nslots++;
    p->p_slots[slot].s_object = obj;
    p->p_slots[slot].s_parent = -1;
    p->p_slots[slot].s_onset = onset;
    p->p_slots[slot].s_time = 0;

        /* objects inside this one (a subpatch's contents) were recorded
        just before it */
    for (i = slot - 1; i > 0 && p->p_slots[i].s_onset >= onset; i--)
        if (p->p_slots[i].s_parent < 0)
            p->p_slots[i].s_parent = slot;
    for (i = onset; i < end; i++)
        if (!p->p_owner[i])
            p->p_owner[i] = slot;
}

static void dsp_tick_profiled(void)
{
    t_dspprofiler *p = THIS->u_profiler;
    t_int *chain = THIS->u_dspchain, *ip = chain;
    uint64_t then = dsp_profiler_now(), now;
    while (ip)
    {
        int index = (int)(ip - chain);
        t_profilerslot *slot = p->p_slots +
            (index < p->p_ownersize ? p->p_owner[index] : 0);
        ip = (*(t_perfroutine)(*ip))(ip);
        now = dsp_profiler_now();
        slot->s_time += now - then;
        then = now;
    }
    p->p_nticks++;
}

void dsp_profiler_setenabled(int enabled)
{
    if (enabled && !THIS->u_profiler)
    {
        t_dspprofiler *p = (t_dspprofiler *)getbytes(sizeof(*p));
        p->p_slotcapacity = 64;
        p->p_slots = (t_profilerslot *)getbytes(
            p->p_slotcapacity * sizeof(*p->p_slots));
        THIS->u_profiler = p;
        dsp_profiler_clear();
            /* rebuild the DSP chain to find out which object added what */
        canvas_update_dsp();
    }
    else if (!enabled)
        dsp_profiler_free();
}

int dsp_profiler_getsize(void)
{
    return (THIS->u_profiler && THIS->u_dspchain ?
        THIS->u_profiler->p_nslots : 0);
}

int dsp_profiler_read(t_dspprofilerentry *entries, int n, double *audiotime,
    unsigned int *generation)
{
    t_dspprofiler *p = THIS->u_profiler;
    t_float sr = sys_getsr();
    int i, count = dsp_profiler_getsize();
    if (!p)
    {
        *audiotime = 0;
        *generation = 0;
        return (0);
    }
    if (count > n)
        count = n;
    for (i = 0; i < p->p_nslots; i++)
    {
        if (i < count)
        {
            entries[i].e_object = p->p_slots[i].s_object;
            entries[i].e_parent = p->p_slots[i].s_parent;
            entries[i].e_time = dsp_profiler_toseconds(p->p_slots[i].s_time);
        }
        p->p_slots[i].s_time = 0;
    }
    *audiotime = (sr > 0 ? (double)p->p_nticks * DEFDACBLKSIZE / sr : 0);
    *generation = p->p_generation;
    p->p_nticks = 0;
    return (count);
}

#endif /* DSP_PROFILER_IMPLEMENTATION */
#endif /* __d_profiler_h_ */
