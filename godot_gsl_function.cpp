#include "godot_gsl_function.h"

/***** GodotGSLFunction methods *****/

GodotGSLFunction::GodotGSLFunction(const String fn)
{
    name = fn;
}

GodotGSLFunction::~GodotGSLFunction()
{
    for (int k = 0; k < instruction_count; k++)
    {
        delete instructions[k];
    }

    if (argv != NULL)
    {
        GGSL_FREE(argv);
    }

    if (instructions != NULL)
    {
        GGSL_FREE(instructions);
    }
}

void GodotGSLFunction::add_arguments(const Array args, GodotGSLMatrix **a)
{
    GGSL_DEBUG_MSG("GodotGSLFunction::add_arguments: procedure starts", 0);
    arg_names = args;

    GodotGSLMatrix **backup = NULL;
    if (argv != NULL)
    {
        backup = argv;
    }

    argv = a;
    GGSL_DEBUG_MSG("GodotGSLFunction::add_arguments: procedure contd", 0);

    if (backup != NULL)
    {
        GGSL_FREE(backup);
    }
    GGSL_DEBUG_MSG("GodotGSLFunction::add_arguments: procedure ends", 0);

    argc = args.size();
}

void GodotGSLFunction::add_instruction(const String in, const Array args)
{
    GodotGSLInstruction **backup = NULL;
    if (instructions != NULL)
    {
        backup = instructions;
    }

    size_t base_size = sizeof(GodotGSLInstruction*);
    size_t mem_size = base_size * (++instruction_count);
    uint8_t* ins_tmp = memnew_arr(uint8_t, mem_size);

    if (backup != NULL)
    {
        memcpy(ins_tmp, backup, mem_size - base_size);
        GGSL_FREE(backup);
    }
    instructions = (GodotGSLInstruction**) ins_tmp;

    GodotGSLInstruction *ins = memnew(GodotGSLInstruction);
    ins->name = in;

    int size = args.size();
    if (size > 0)
    {
        set_instruction_arguments(ins, args);
    }

    GGSL_MESSAGE_I(String("GodotGSLFunction::add_instruction: Count is {_}").format(instruction_count).ascii(), 1);
    printf("aaaaaaaaaaaaaaaaa %d\n", instruction_count);
    instructions[instruction_count - 1] = ins;
    current = ins;
}

void GodotGSLFunction::set_instruction_arguments(Array args)
{
    if (current == NULL)
    {
        GGSL_MESSAGE("GodotGSLFunction::set_instruction_arguments: current == NULL");
        return;
    }

    set_instruction_arguments(current, args);
}

void GodotGSLFunction::set_instruction_arguments(GodotGSLInstruction *ins, Array args)
{
    if (ins == NULL)
    {
        GGSL_MESSAGE("GodotGSLFunction::set_instruction_arguments: ins == NULL");
        return;
    }

    int size = args.size();
    if (size == 0 || size > argc)
    {
        GGSL_MESSAGE("GodotGSLFunction::set_instruction_arguments: size == 0 || size > argc");
        return;
    }

    if (ins->argv != NULL)
    {
        GGSL_FREE(ins->argv);
    }

    GodotGSLMatrix **argv;
    GGSL_ALLOC(ins->argv, size);

    for (int k = 0; k < size; k++)
    {
        String vn = args[k];
        int index = arg_names.find(vn);

        if (index > -1)
        {
            ins->argv[k] = argv[index];
        }
        else
        {
            GGSL_MESSAGE("GodotGSLFunction::set_instruction_arguments: !arg_names.has(vn)");
            return;
        }
    }

    ins->argc = size;
}

void GodotGSLFunction::execute()
{
    if (instructions == NULL)
    {
        GGSL_MESSAGE("GodotGSLFunction::execute: instructions == NULL");
        return;
    }

    GodotGSLInstruction *first = instructions[0];

    do {
        first->execute();
    } while ((first = first->next()));
}
