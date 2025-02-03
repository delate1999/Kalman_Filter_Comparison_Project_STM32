/*yafl_config.h*/

#ifndef YAFL_CONFIG_H
#define YAFL_CONFIG_H

#include <math.h>
#include <stdint.h>

#ifdef DEBUG
    /*
    In this example we will use standard output.
    You can actually use any printf implementation you want.
    */
#   include <stdio.h>
#   define YAFL_DBG(...) fprintf(stderr, __VA_ARGS__)

    /*
    Using branch speculation may save some clocks...
    */
#   ifdef __GNUC__
#       define YAFL_UNLIKELY(x) __builtin_expect((x), 0)
#   else /*__GNUC__*/
#       define YAFL_UNLIKELY(x) (x)
#   endif/*__GNUC__*/
#else /*DEBUG*/
#   define YAFL_DBG(...) /*Do nothing here*/
    /*
    Here we have "Never" actually, but you can use some of above definitions if you want.
    */
#   define YAFL_UNLIKELY(x) (0)
#endif/*DEBUG*/

#define YAFL_EPS  (1.0e-6)
#define YAFL_SQRT sqrtf
#define YAFL_ABS  fabsf
#define YAFL_EXP  expf
#define YAFL_LOG  logf


typedef double   yaflFloat;
typedef int32_t   yaflInt;

#endif/*YAFL_CONFIG_H*/