#pragma once
#include <m_pd.h>

typedef void (*pd_gui_callback)(void*, char const*, int, t_atom*);
typedef void (*pd_message_callback)(void*, void*, t_symbol*, int, t_atom*);

void register_gui_triggers(t_pdinstance* instance, void* target, pd_gui_callback gui_callback, pd_message_callback message_callback);

void set_instance_lock(const void* lock, void(*lock_func)(void*), void(*unlock_func)(void*), void(*clear_references_func)(void*, t_pd*));

void update_gui();

void plugdata_gui_message(const char* message, va_list args);

void plugdata_forward_message(void *x, t_symbol *s, int argc, t_atom *argv);
