#ifndef GODOT_GSL_FUNCTION_H
#define GODOT_GSL_FUNCTION_H

#include "godot_gsl_matrix.h"

class GodotGSLFunction
{
public:
    class GodotGSLInstruction
    {
    public:
        GodotGSLInstruction() { }

        GodotGSLInstruction(const String in)
        {
            name = in;
        }

        ~GodotGSLInstruction()
        {
            if (argv != NULL)
            {
                GGSL_FREE(argv);
            }
        }

        GodotGSLInstruction *next()
        {
            return nxt;
        }

        GodotGSLInstruction *next(GodotGSLInstruction *ins)
        {
            if (nxt != NULL)
            {
                GodotGSLInstruction *backup = nxt;
                nxt = ins;
                ins->next(backup);
            }
            else
            {
                nxt = ins;
            }

            return ins;
        }

        void execute()
        {
            if (argc == 2)
            {
                GodotGSLMatrix *mtx_a = argv[0];
                GodotGSLMatrix *mtx_to = argv[1];
                mtx_a->fx(name, mtx_to, argv_bounds);
            }
            else if (argc == 3)
            {
                GodotGSLMatrix *mtx_a = argv[0];
                GodotGSLMatrix *mtx_b = argv[1];
                GodotGSLMatrix *mtx_to = argv[2];
                mtx_a->fx(name, mtx_b, mtx_to, argv_bounds);
            }
            else
            {
                GGSL_MESSAGE("GodotGSLInstruction::execute: Wrong number of arguments");
                return;
            }
        }

    private:
        friend class GodotGSLFunction;
        GodotGSLInstruction *nxt = NULL;
        String name;
        size_t argc;
        GodotGSLMatrix **argv = NULL;
        GGSL_BOUNDS *argv_bounds;
    };

    GodotGSLFunction() { }
    GodotGSLFunction(const String fn);
    ~GodotGSLFunction();
    void execute();
    void add_arguments(const Array args, GodotGSLMatrix **a);
    void add_argument(const String vn, GodotGSLMatrix *a);
    void add_instruction(const String in);
    void add_instruction(const String in, const Array args);
    void set_instruction_arguments(const Array args);
    void set_instruction_arguments(GodotGSLInstruction *ins, const Array args);
    GodotGSLMatrix *get_arg(const String vn);
    GodotGSLMatrix *get_arg(const size_t index);
    size_t get_argc() { return argc; };

private:
    Array arg_names;
    String name;
    GodotGSLMatrix **argv = NULL;
    GodotGSLInstruction **instructions = NULL;
    size_t instruction_count = 0;
    size_t argc = 0;
    GodotGSLInstruction *current = NULL;
};

#endif
