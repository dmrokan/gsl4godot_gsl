#ifndef GODOT_GSL_MACROS_H
#define GODOT_GSL_MACROS_H

#include "core/error_macros.h"

#define GGSL_DEBUG_ON

#define GGSL_MESSAGE(msg) \
    {                     \
        WARN_PRINT(msg)   \
    }                     \

#define GGSL_MESSAGE_I(msg, lvl) \
    {                            \
        int k = 0;               \
        do {                     \
            printf(" ");         \
        } while (lvl > k++);     \
        WARN_PRINT(msg);         \
    }                            \

#define GGSL_ALLOC(__var, size)                                 \
    {                                                           \
        size_t size_in_bytes = size * sizeof(GodotGSLMatrix*);            \
        uint8_t *argv_tmp = memnew_arr(uint8_t, size_in_bytes); \
        __var = (GodotGSLMatrix**) argv_tmp;                              \
    }                                                           \

#define GGSL_FREE(__var)                        \
    {                                           \
        memdelete_arr((uint8_t*) __var);        \
    }                                           \

#ifdef GGSL_DEBUG_ON
#define GGSL_DEBUG_MSG(msg, lvl) \
    {                            \
        int k = 0;               \
        do {                     \
            printf(" DBG:");     \
        } while(lvl > k++);      \
        WARN_PRINT(msg);         \
    }                            \

#else
#define GGSL_DEBUG_MSG(msg, lvl) 0
#endif

#endif
