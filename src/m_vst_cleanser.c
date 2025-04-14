#include "m_pd.h"

#ifdef VST_CLEANSER
#    undef s_pointer
#    undef s_float
#    undef s_symbol
#    undef s_bang
#    undef s_list
#    undef s_anything
#    undef s_signal
#    undef s__N
#    undef s__X
#    undef s_x
#    undef s_y
#    undef s_

#    if _MSC_VER
__declspec(dllexport)
#    endif
t_symbol s_pointer,
    s_float, s_symbol, s_bang, s_list, s_anything,
    s_signal, s__N, s__X, s_x, s_y, s_;

void vst_cleanser(t_symbol** s)
{
    if (*s >= &s_pointer && *s <= &s_) {
        *s = gensym((*s)->s_name);
    }
}
#endif
