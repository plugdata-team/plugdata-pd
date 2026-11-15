#pragma once
#include <m_pd.h>
#include <stdarg.h>

typedef void (*pd_gui_callback)(void*, char const*, int, t_atom*);
typedef void (*pd_message_callback)(void*, int, void*, t_symbol*, int, t_atom*);
/* Same arguments as pd_message_callback: target, type, object or connection, selector, argc, argv */
typedef void (*pd_trace_callback)(void*, int, void*, t_symbol*, int, t_atom*);

void register_gui_triggers(t_pdinstance* instance, void* target, pd_gui_callback gui_callback, pd_message_callback message_callback);
/* Sets a second callback that gets every message plugdata_forward_message() forwards while debugging is on,
   in the order Pd sends them. A NULL target removes it. Call it with Pd's lock held, like Pd calls the callback */
void register_trace_trigger(t_pdinstance* instance, void* target, pd_trace_callback trace_callback);

void setup_lock(const void* lock, void(*lock_func)(void*), void(*unlock_func)(void*));

void setup_weakreferences(void(*clear_references_func)(void*, void*), void(*register_reference_func)(void*, void*, void*), void(*unregister_reference_func)(void*, void*, void*), int(*is_reference_valid_func)(void*));

void plugdata_gui_message(const char* message, va_list args);

void plugdata_forward_message(int type, void *x, t_symbol *s, int argc, t_atom *argv);

int plugdata_debugging_enabled();
int plugdata_activity_enabled();
int plugdata_debugging_or_activity_enabled();

void set_plugdata_debugging_enabled(int enabled_debugging);
void set_plugdata_activity_enabled(int enabled_activity);

#ifndef __cplusplus
#include <stdatomic.h>

/* The flags behind plugdata_activity_enabled() and plugdata_debugging_or_activity_enabled(), for code that
   checks them on every message (outlets, inlets). Reading them inline avoids a function call per check,
   which cost ~25% of message throughput when called for every connection */
extern atomic_bool plugdata_activity_flag;
extern atomic_bool plugdata_debug_or_activity_flag;

static inline int plugdata_activity_enabled_inline(void)
{
    return atomic_load_explicit(&plugdata_activity_flag, memory_order_relaxed);
}

static inline int plugdata_debugging_or_activity_enabled_inline(void)
{
    return atomic_load_explicit(&plugdata_debug_or_activity_flag, memory_order_relaxed);
}
#endif
