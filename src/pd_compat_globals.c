/* pd_compat_globals.c
 *
 * Compatibility layer for third-party Pd externals compiled against vanilla Pd
 * (non-PDINSTANCE mode). In PDINSTANCE mode, symbols like s_bang, s_float, etc.
 * become per-instance macros and are NOT exported as global variables. This file
 * provides real global t_symbol variables that satisfy the linker, and a
 * vst_cleanser() function that translates these "stale" global pointers to the
 * correct per-instance symbols at Pd API entry points.
 *
 * See m_pd.h lines 1127-1137 for the VST_CLEANSER design by Miller Puckette.
 */

#include "m_pd.h"

#if defined(PDINSTANCE) && defined(VST_CLEANSER)

/* Undef the per-instance macros so we can define real global variables.
 * After this point, s_bang etc. refer to our globals, NOT pd_this->pd_s_bang.
 * We access the per-instance versions via pd_this->pd_s_xxx directly. */
#undef s_pointer
#undef s_float
#undef s_symbol
#undef s_bang
#undef s_list
#undef s_anything
#undef s_signal
#undef s__N
#undef s__X
#undef s_x
#undef s_y
#undef s_

/* Global t_symbol variables exported from pd.dll.
 * Initialized with correct s_name so that direct dereferences (e.g.
 * s_bang.s_name) return the expected string. The addresses of these
 * globals serve as identity markers for vst_cleanser(). */
t_symbol s_pointer  = {"pointer",  0, 0};
t_symbol s_float    = {"float",    0, 0};
t_symbol s_symbol   = {"symbol",   0, 0};
t_symbol s_bang     = {"bang",     0, 0};
t_symbol s_list     = {"list",     0, 0};
t_symbol s_anything = {"anything", 0, 0};
t_symbol s_signal   = {"signal",   0, 0};
t_symbol s__N       = {"N",        0, 0};
t_symbol s__X       = {"X",        0, 0};
t_symbol s_x        = {"x",        0, 0};
t_symbol s_y        = {"y",        0, 0};
t_symbol s_         = {"",         0, 0};

/* vst_cleanser - translate global symbol pointers to per-instance equivalents.
 *
 * Called at Pd API entry points (pd_bind, pd_unbind, pd_symbol, pd_list,
 * pd_anything, class_doaddmethod) to intercept "stale" global symbol pointers
 * from vanilla-Pd-compiled externals and replace them with the correct
 * per-instance symbol from the current pd_this. */
void vst_cleanser(t_symbol **s)
{
    t_symbol *sym = *s;
    if      (sym == &s_pointer)  *s = &(pd_this->pd_s_pointer);
    else if (sym == &s_float)    *s = &(pd_this->pd_s_float);
    else if (sym == &s_symbol)   *s = &(pd_this->pd_s_symbol);
    else if (sym == &s_bang)     *s = &(pd_this->pd_s_bang);
    else if (sym == &s_list)     *s = &(pd_this->pd_s_list);
    else if (sym == &s_anything) *s = &(pd_this->pd_s_anything);
    else if (sym == &s_signal)   *s = &(pd_this->pd_s_signal);
    else if (sym == &s__N)       *s = &(pd_this->pd_s__N);
    else if (sym == &s__X)       *s = &(pd_this->pd_s__X);
    else if (sym == &s_x)        *s = &(pd_this->pd_s_x);
    else if (sym == &s_y)        *s = &(pd_this->pd_s_y);
    else if (sym == &s_)         *s = &(pd_this->pd_s_);
}

#endif /* PDINSTANCE && VST_CLEANSER */
