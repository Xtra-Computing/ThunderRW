//
// Created by Shixuan Sun on 2020/10/17.
//

#ifndef TYPE_CONFIG_H
#define TYPE_CONFIG_H

#include <limits.h>
#include <stdint.h>

#if defined(LONG)

typedef long intT;
typedef unsigned long uintT;
#define INT_T_MAX LONG_MAX
#define UINT_T_MAX ULONG_MAX

#else

typedef int intT;
typedef unsigned int uintT;
typedef struct {
    double alias_value_;
    intT first_;
    intT second_;
} AliasSlot;

#define INT_T_MAX INT_MAX
#define UINT_T_MAX UINT_MAX

#endif


#if defined(EDGE_LONG)

typedef long intE;
typedef unsigned long uintE;
#define INT_E_MAX LONG_MAX
#define UINT_E_MAX ULONG_MAX

#else

typedef int intE;
typedef unsigned int uintE;
#define INT_E_MAX INT_MAX
#define UINT_E_MAX UINT_MAX

#endif


#endif