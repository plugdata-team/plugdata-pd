/* Copyright (c) 2026 Timothy Schoen
* For information on usage and redistribution, and for a DISCLAIMER OF ALL
* WARRANTIES, see the file, "LICENSE.txt," in this distribution.  */

EXTERN int plugdata_object_probe_enabled();
EXTERN void set_plugdata_object_probe_enabled(int);

#define PROBE_IDX_BANG    0
#define PROBE_IDX_FLOAT   1
#define PROBE_IDX_SYMBOL  2
#define PROBE_IDX_POINTER 3
#define PROBE_IDX_LIST    4
#define PROBE_IDX_ANY     5

static void probe_bang(t_pd *x)
{
    t_class *c = *((t_class **)x);
    ((t_bangmethod)c->c_methods[PROBE_IDX_BANG].me_fun)(x);
    plugdata_forward_message(x, &s_bang, 0, NULL);
}

static void probe_float(t_pd *x, t_float f)
{
    t_class *c = *((t_class **)x);
    ((t_floatmethod)c->c_methods[PROBE_IDX_FLOAT].me_fun)(x, f);
    t_atom a; SETFLOAT(&a, f);
    plugdata_forward_message(x, &s_float, 1, &a);
}

static void probe_symbol(t_pd *x, t_symbol *s)
{
    t_class *c = *((t_class **)x);
    ((t_symbolmethod)c->c_methods[PROBE_IDX_SYMBOL].me_fun)(x, s);
    t_atom a; SETSYMBOL(&a, s);
    plugdata_forward_message(x, &s_symbol, 1, &a);
}

static void probe_pointer(t_pd *x, t_gpointer *gp)
{
    t_class *c = *((t_class **)x);
    ((t_pointermethod)c->c_methods[PROBE_IDX_POINTER].me_fun)(x, gp);
    t_atom a; SETPOINTER(&a, gp);
    plugdata_forward_message(x, &s_pointer, 1, &a);
}

static void probe_list(t_pd *x, t_symbol *s, int argc, t_atom *argv)
{
    t_class *c = *((t_class **)x);
    ((t_listmethod)c->c_methods[PROBE_IDX_LIST].me_fun)(x, s, argc, argv);
    plugdata_forward_message(x, &s_list, argc, argv);
}

static void probe_any(t_pd *x, t_symbol *s, int argc, t_atom *argv)
{
    t_class *c = *((t_class **)x);
    ((t_anymethod)c->c_methods[PROBE_IDX_ANY].me_fun)(x, s, argc, argv);
    plugdata_forward_message(x, s, argc, argv);
}

static inline void probe_named_wrapper(t_pd *x, t_symbol *s, int argc, t_atom *argv)
{
    /* Call the original, stored under "$sel" */
    char mangled[MAXPDSTRING];
    mangled[0] = '$';
    strncpy(mangled + 1, s->s_name, MAXPDSTRING - 2);
    mangled[MAXPDSTRING - 1] = '\0';
    pd_typedmess(x, gensym(mangled), argc, argv);
    plugdata_forward_message(x, s, argc, argv);
}

static inline void plugdata_fwd_bang(t_class *c, t_method fn)
{
    c->c_methods[PROBE_IDX_BANG].me_fun = (t_gotfn)fn;
    c->c_bangmethod = probe_bang;
}

static inline void plugdata_fwd_float(t_class *c, t_method fn)
{
    c->c_methods[PROBE_IDX_FLOAT].me_fun = (t_gotfn)fn;
    c->c_floatmethod = probe_float;
}

static inline void plugdata_fwd_symbol(t_class *c, t_method fn)
{
    c->c_methods[PROBE_IDX_SYMBOL].me_fun = (t_gotfn)fn;
    c->c_symbolmethod = probe_symbol;
}

static inline void plugdata_fwd_pointer(t_class *c, t_method fn)
{
    c->c_methods[PROBE_IDX_POINTER].me_fun = (t_gotfn)fn;
    c->c_pointermethod = probe_pointer;
}

static inline void plugdata_fwd_list(t_class *c, t_method fn)
{
    c->c_methods[PROBE_IDX_LIST].me_fun = (t_gotfn)fn;
    c->c_listmethod = probe_list;
}

static inline void plugdata_fwd_anything(t_class *c, t_method fn)
{
    c->c_methods[PROBE_IDX_ANY].me_fun = (t_gotfn)fn;
    c->c_anymethod = probe_any;
}
