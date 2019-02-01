/* godot_gsl.cpp */

#include "godot_gsl.h"

void GodotGSL::_bind_methods()
{
    ClassDB::bind_method(D_METHOD("matrix_new", "vn", "row_count", "col_count"), &GodotGSL::matrix_new);
    ClassDB::bind_method(D_METHOD("matrix_new_from_array", "vn", "a"), &GodotGSL::matrix_new_from_array);
    ClassDB::bind_method(D_METHOD("add", "vn1", "vn2"), &GodotGSL::add);
    ClassDB::bind_method(D_METHOD("to_array", "vn"), &GodotGSL::to_array);
    ClassDB::bind_method(D_METHOD("set_zero", "vn"), &GodotGSL::set_zero);
    ClassDB::bind_method(D_METHOD("prod", "vn1", "vn2", "to"), &GodotGSL::prod);
    ClassDB::bind_method(D_METHOD("vector_new", "vn", "size", "is_row_vector"), &GodotGSL::vector_new, NULL, NULL, false);
    ClassDB::bind_method(D_METHOD("vector_new_from_array", "vn", "a", "is_row_vector"), &GodotGSL::vector_new_from_array, NULL, NULL, false);
    ClassDB::bind_method(D_METHOD("fx", "fn", "a", "to"), &GodotGSL::fx, "=", "", "");
    ClassDB::bind_method(D_METHOD("function", "fn", "args"), &GodotGSL::function, "", Array());
    ClassDB::bind_method(D_METHOD("instruction", "ins_name", "args"), &GodotGSL::instruction, "", Array());
    ClassDB::bind_method(D_METHOD("callf", "fn"), &GodotGSL::callf);
    ClassDB::bind_method(D_METHOD("ode", "on", "dim"), &GodotGSL::ode);
    ClassDB::bind_method(D_METHOD("ode_set_fx", "on", "fn"), &GodotGSL::ode_set_fx);
    ClassDB::bind_method(D_METHOD("ode_run", "on", "x0", "time_interval", "dt"), &GodotGSL::ode_run);
    ClassDB::bind_method(D_METHOD("ode_set_node_path", "on", "object", "property", "index"), &GodotGSL::ode_set_node_path);
    ClassDB::bind_method(D_METHOD("ode_set_init_cond", "on", "x0", "t0"), &GodotGSL::ode_set_init_cond);
    ClassDB::bind_method(D_METHOD("ode_run_delta", "on", "dt"), &GodotGSL::ode_run_delta);
}

GodotGSL::GodotGSL()
{
    return;
}

GodotGSL::~GodotGSL()
{
    variables.clear();
    functions.clear();
    odes.clear();
}

void GodotGSL::_add_variable(String vn, GodotGSLMatrix* mtx)
{
    if (variables.has(vn))
    {
        variables.erase(vn);
        variables[vn] = mtx;
    }
    else
    {
        variables.insert(vn, mtx);
    }
}

void GodotGSL::matrix_new(String vn, const size_t row_count, const size_t col_count)
{
    GodotGSLMatrix *mtx = memnew(GodotGSLMatrix);
    mtx->init(row_count, col_count);

    _add_variable(vn, mtx);
}

void GodotGSL::matrix_new_from_array(String vn, const Array a)
{
    GodotGSLMatrix *mtx = memnew(GodotGSLMatrix);
    mtx->init(a);

    _add_variable(vn, mtx);
}

void GodotGSL::add(String vn1, String vn2)
{
    GodotGSLMatrix *mtx1 = variables[vn1];
    GodotGSLMatrix *mtx2 = variables[vn2];

    mtx1->add(*mtx2, NULL);
}

void GodotGSL::set_zero(String vn)
{
    GodotGSLMatrix *mtx = variables[vn];

    mtx->set_zero();
}

Array GodotGSL::to_array(String vn)
{
    GodotGSLMatrix *mtx = variables[vn];

    size_t row_count = mtx->size[0];
    size_t col_count = mtx->size[1];

    Array rows;
    for (size_t k = 0; k < row_count; k++)
    {
        Array row;
        for (size_t l = 0; l < col_count; l++)
        {
            row.append(mtx->get(k, l));
        }

        rows.append(row);
    }

    return rows;
}

void GodotGSL::prod(String vn1, String vn2, String to)
{
    GodotGSLMatrix *mtx1 = variables[vn1];
    GodotGSLMatrix *mtx2 = variables[vn2];

    GodotGSLMatrix *result = NULL;
    mtx2->prod(*mtx1, result, NULL);

    _add_variable(to, result);
}


void GodotGSL::vector_new(String vn, const size_t size, bool is_row_vector)
{
    GodotGSLMatrix *mtx = memnew(GodotGSLMatrix);

    size_t row_count = size;
    size_t col_count = 1;

    if (is_row_vector)
    {
        row_count = 1;
        col_count = size;
    }

    mtx->init(row_count, col_count);

    _add_variable(vn, mtx);
}

void GodotGSL::vector_new_from_array(String vn, const Array a, bool is_row_vector)
{
    GodotGSLMatrix *mtx = memnew(GodotGSLMatrix);

    Array b;

    if (is_row_vector)
    {
        b.append(a);
    }
    else
    {
        size_t size = a.size();
        for (size_t k = 0; k < size; k++)
        {
            STYPE value = a.get(k);
            Array c;
            c.append(value);
            b.append(c);
        }
    }

    mtx->init(b);

    _add_variable(vn, mtx);
}

void GodotGSL::fx(const String fn, const String a, String to)
{
    GodotGSLMatrix *mtx_a;
    GodotGSLMatrix *mtx_to;

    if (fn == FN_EQ)
    {
        if (to.empty())
        {
            return;
        }
    }
    else if (fn == FN_SIN)
    {
        if (!variables.has(a))
        {
            GGSL_MESSAGE("GodotGSL::fx: !variables.has(a)");
            return;
        }

        if (to.empty())
        {
            to = a;
        }
        else
        {
            if (!variables.has(to))
            {
                GGSL_MESSAGE("GodotGSL::fx: !variables.has(to)");
                return;
            }
        }

        mtx_a = variables[a];
        mtx_to = variables[to];
        /* TODO: find a way to supply default bounds
         * Possibly you are going to put this in fx
         * Cause you have all bound informations of matrices
         */
        mtx_a->fx(fn, mtx_to, NULL);
    }

    /*
     * TODO: If there is no 'to' variable
     * It may be added automaticall
     * SOON

    _add_variable(to, mtx_to);

    */
}

void GodotGSL::_add_function(const String fn, GodotGSLFunction *fnc)
{
    if (functions.has(fn))
    {
        functions.erase(fn);
        functions[fn] = fnc;
    }
    else
    {
        functions.insert(fn, fnc);
    }
}

void GodotGSL::function(const String fn, const Array args)
{
    GodotGSLFunction *fnc = memnew(GodotGSLFunction(fn));
    _add_function(fn, fnc);

    int size = args.size();
    if (size > 0)
    {
        GodotGSLMatrix **argv;
        GGSL_ALLOC(argv, size);

        for (int k = 0; k < size; k++)
        {
            String vn = args[k];
            // int sindex = vn.find_char('[');
            // if (sindex < 0)
            // {
            //     sindex = vn.size();
            // }

            // vn = vn.substr(0, sindex);
            if (variables.has(vn))
            {
                argv[k] = variables[vn];
                GGSL_DEBUG_MSG(vn.utf8().get_data(), 1);
            }
            else
            {
                GGSL_MESSAGE("GodotGSL::function: !variables.has(vn)");
                GGSL_MESSAGE_I(vn.utf8().get_data(), 1);
                return;
            }
        }

        fnc->add_arguments(args, argv);
        GGSL_DEBUG_MSG("GodotGSL::function: args are added", 0);
    }

    current = fnc;
    GGSL_DEBUG_MSG("GodotGSL::function: procedure ends", 0);
}

void GodotGSL::instruction(const String in, const Array args)
{
    if (current == NULL)
    {
        GGSL_MESSAGE("GodotGSL::instruction: current == NULL");
        return;
    }

    current->add_instruction(in, args);
    GGSL_DEBUG_MSG("GodotGSL::instruction: procedure ends", 0);
}

void GodotGSL::callf(const String fn)
{
    GGSL_DEBUG_MSG("GodotGSL::callf: procedure starts", 0);
    if (!functions.has(fn))
    {
        GGSL_MESSAGE("GodotGSL::call: !functions.has(fn)");
        return;
    }

    GodotGSLFunction *fnc = functions[fn];
    fnc->execute();
    GGSL_DEBUG_MSG("GodotGSL::callf: procedure ends", 0);
}

void GodotGSL::ode(const String on, const size_t dim)
{
    if (odes.has(on))
    {
        GodotGSLODE *to_rm = odes[on];
        odes[on] = memnew(GodotGSLODE(dim));
        delete to_rm;
    }
    else
    {
        odes[on] = memnew(GodotGSLODE(dim));
    }
}

void GodotGSL::ode_set_fx(const String on, const String fn)
{
    GodotGSLODE *ode;
    GodotGSLFunction *fnc;

    if (odes.has(on))
    {
        ode = odes[on];
    }
    else
    {
        GGSL_MESSAGE("GodotGSL::ode_set_fx: !odes.has(on)");
        return;
    }

    if (functions.has(fn))
    {
        fnc = functions[fn];
    }
    else
    {
        GGSL_MESSAGE("GodotGSL::ode_set_fx: !functions.has(on)");
        return;
    }

    ode->set_function(fnc);
}

void GodotGSL::ode_run(const String on, const Array x0, const Array time_interval, const double dt)
{
    if (!odes.has(on))
    {
        GGSL_MESSAGE("GodotGSL::ode_run: !odes.has(on)");
        return;
    }

    GodotGSLODE *ode = odes[on];
    ode->x->copy_vector_from_array(x0);

    double start_time, end_time;
    if (time_interval.size() == 1)
    {
        start_time = 0.0;
        end_time = time_interval[0];
    }
    else if (time_interval.size() == 2)
    {
        start_time = time_interval[0];
        end_time = time_interval[1];
    }

    ode->run(start_time, end_time, dt);
}

void GodotGSL::ode_set_init_cond(const String on, const Array x0, const double t0)
{
    if (!odes.has(on))
    {
        GGSL_MESSAGE("GodotGSL::ode_run: !odes.has(on)");
        return;
    }

    GodotGSLODE *ode = odes[on];
    ode->set_initial_conditions(x0, t0);
}

void GodotGSL::ode_run_delta(const String on, const double delta)
{
    if (!odes.has(on))
    {
        GGSL_MESSAGE("GodotGSL::ode_run: !odes.has(on)");
        return;
    }

    GodotGSLODE *ode = odes[on];

    ode->step(delta);
}

void GodotGSL::ode_set_node_path(const String on, Variant obj_var, const String subpath, const int index)
{
    if (!odes.has(on))
    {
        GGSL_MESSAGE("GodotGSL::ode_set_node_path: !odes.has(on)");
        return;
    }

    Object *obj = (Object*) obj_var;

    GodotGSLODE *ode = odes[on];
    /* TODO: Use weak ptr */
    ode->set_node_path(obj, subpath, index);
}
