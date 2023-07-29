#pragma once
#include <m_pd.h>

typedef void (*pd_gui_callback)(void*, char const*, t_atom*, t_atom*, t_atom*);
typedef void (*pd_message_callback)(void*, void*, t_symbol*, int, t_atom*);

void register_gui_triggers(t_pdinstance* instance, void* target, pd_gui_callback gui_callback, pd_message_callback message_callback);

void set_lock(const void* lock, void(*lock_func)(void*), void(*unlock_func)(void*));

void set_global_lock(void(*instance_read_lock_func)(), void(*instance_read_unlock_func)(), void(*instance_write_lock_func)(), void(*instance_write_unlock_func)());

void set_free_callback(void(*free_callback)(void*, t_pd*));

void create_panel(int openpanel, char const* path, char const* snd);
void trigger_open_file(const char* file);
void update_gui();

void plugdata_forward_message(void *x, t_symbol *s, int argc, t_atom *argv);
