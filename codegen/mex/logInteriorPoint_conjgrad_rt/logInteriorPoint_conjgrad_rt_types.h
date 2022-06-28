/*
 * logInteriorPoint_conjgrad_rt_types.h
 *
 * Code generation for function 'logInteriorPoint_conjgrad_rt'
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

#ifndef typedef_c_logInteriorPoint_conjgrad_rtS
#define typedef_c_logInteriorPoint_conjgrad_rtS

typedef struct {
  struct {
    real_T invW[10000];
  } f0;
} c_logInteriorPoint_conjgrad_rtS;

#endif                                 /*typedef_c_logInteriorPoint_conjgrad_rtS*/

/* End of code generation (logInteriorPoint_conjgrad_rt_types.h) */
