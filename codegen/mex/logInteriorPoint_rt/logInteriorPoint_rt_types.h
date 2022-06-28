/*
 * logInteriorPoint_rt_types.h
 *
 * Code generation for function 'logInteriorPoint_rt'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

#ifndef typedef_logInteriorPoint_rtStackData
#define typedef_logInteriorPoint_rtStackData

typedef struct {
  struct {
    real_T invW[10000];
    real_T G[10000];
    real_T b_G[10000];
  } f0;
} logInteriorPoint_rtStackData;

#endif                                 /*typedef_logInteriorPoint_rtStackData*/

/* End of code generation (logInteriorPoint_rt_types.h) */
