#include "godot_gsl_function.h"

#define ARGV_BOUND_DATA_COUNT 4
#define GGSL_ARGV_BOUNDS_ALLOC(__var, size)                           \
    {                                                                 \
    size_t byte_size = size * sizeof(size_t) * ARGV_BOUND_DATA_COUNT; \
    uint8_t *mem_tmp = memnew_arr(uint8_t, byte_size);                \
    __var = (size_t**) mem_tmp;                                       \
    }                                                                 \


#define GGSL_ARGV_BOUNDS_SET(__var, __bounds, k)                        \
    {                                                                   \
        memcpy(&__var[k], __bounds,                                     \
               sizeof(size_t) * ARGV_BOUND_DATA_COUNT);                 \
    }                                                                   \


typedef enum GGSL_PARSER_STATE {
    S_IDLE,
    S_NUMBER,
    S_SEP,
} GGSL_PARSER_STATE;

size_t *argv_bounds_parse(const String vn, size_t *size)
{
    static GGSL_PARSER_STATE state = S_IDLE;
    static size_t bounds[ARGV_BOUND_DATA_COUNT];
    bounds[0] = bounds[2] = 0;
    bounds[1] = size[0];
    bounds[3] = size[1];

    int l = 0;
    String sub = "";
    for (int k = 0; k < vn.size(); k++)
    {
        char c = (vn[k]);

        if (state == S_IDLE)
        {
            if (c == '[')
            {
                state = S_NUMBER;
                continue;
            }
        }
        else if (state == S_NUMBER)
        {
            if (c >= '0' && c <= '9')
            {
                sub += c;
                continue;
            }
            else if (c == ':' || c == ',')
            {
                bounds[l] = sub.to_int();
                l++;
                state = S_SEP;
                continue;
            }
            else if (c == ']')
            {
                bounds[l] = sub.to_int();
                state = S_IDLE;
                break;
            }
        }
        else if (state == S_SEP)
        {
            sub.clear();
            k--;
            state = S_NUMBER;
            continue;
        }
    }

    return &bounds[0];
}

String remove_whitespace(const String str)
{
    String result = str;
    result.replace(" ", "");
    result.replace("\t", "");
    result.replace("\n", "");
    result.replace("\r", "");

    return result;
}


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

void GodotGSLFunction::add_argument(const String vn, GodotGSLMatrix *a)
{
    GodotGSLMatrix **_argv;
    GGSL_ALLOC(_argv, argc + 1);

    memcpy(_argv, argv, argc * sizeof(GodotGSLMatrix*));

    argv[argc] = a;
    arg_names.append(vn);

    argc++;
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
    
    instructions[instruction_count - 1] = ins;

    if (current == NULL)
    {
        current = ins;
    }
    else
    {
        current->nxt = ins;
        current = ins;
    }
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
        memdelete_arr(ins->argv_bounds);
    }

    GGSL_ALLOC(ins->argv, size);
    ins->argv_bounds = memnew_arr(GGSL_BOUNDS, size);

    for (int k = 0; k < size; k++)
    {
        String arg = remove_whitespace(args[k]);
        int sindex = arg.find_char('[');
        if (sindex < 0)
        {
            sindex = arg.size();
        }

        String vn = arg.substr(0, sindex);
        int index = arg_names.find(vn);

        if (index > -1)
        {
            ins->argv[k] = argv[index];
            GodotGSLMatrix *argv_mtx = (GodotGSLMatrix*) argv[index];
            size_t *bounds = argv_bounds_parse(arg, argv_mtx->size);
            /* TODO: memcpy may no be the right wway to do this */
            GGSL_ARGV_BOUNDS_SET(ins->argv_bounds, bounds, k);
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

GodotGSLMatrix *GodotGSLFunction::get_arg(const String vn)
{
    int index = arg_names.find(vn);

    if (index < 0)
    {
        GGSL_ERR_MESSAGE("GodotGSLFunction::get_arg: index < 0");
        return NULL;
    }

    return argv[index];
}

GodotGSLMatrix *GodotGSLFunction::get_arg(const size_t index)
{
    if (index >= argc)
    {
        GGSL_ERR_MESSAGE("GodotGSLFunction::get_arg: index >= argc");
        return NULL;
    }

    return argv[index];
}
