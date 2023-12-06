#pragma once
#include <m_pd.h>

typedef void (*pd_gui_callback)(void*, char const*, int, t_atom*);
typedef void (*pd_message_callback)(void*, void*, t_symbol*, int, t_atom*);

void register_gui_triggers(t_pdinstance* instance, void* target, pd_gui_callback gui_callback, pd_message_callback message_callback);

void setup_lock(const void* lock, void(*lock_func)(void*), void(*unlock_func)(void*));

void setup_weakreferences(void(*clear_references_func)(void*, void*), void(*register_reference_func)(void*, void*, void*), void(*unregister_reference_func)(void*, void*, void*), int(*is_reference_valid_func)(void*));

void plugdata_gui_message(const char* message, va_list args);

void plugdata_forward_message(void *x, t_symbol *s, int argc, t_atom *argv);
