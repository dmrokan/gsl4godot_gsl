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
}

GodotGSL::GodotGSL()
{
    return;
}

GodotGSL::~GodotGSL()
{
    variables.clear();
    functions.clear();
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

    mtx1->add(*mtx2);
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

    GodotGSLMatrix *result = mtx2->prod(*mtx1);
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
        mtx_a->fx(fn, mtx_to);
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
            if (variables.has(vn))
            {
                argv[k] = variables[vn];
                GGSL_DEBUG_MSG(vn.ascii(), 1);
            }
            else
            {
                GGSL_MESSAGE("GodotGSL::function: !variables.has(vn)");
                GGSL_MESSAGE_I(vn.ascii(), 1);
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
