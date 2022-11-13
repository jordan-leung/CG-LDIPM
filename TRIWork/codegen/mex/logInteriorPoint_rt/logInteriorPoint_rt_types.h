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
#ifndef typedef_logInteriorPoint_rtStackData
#define typedef_logInteriorPoint_rtStackData

typedef struct {
  union
  {
    struct {
      real_T Q[14400];
      real_T y[14400];
    } f0;

    struct {
      real_T Q[14400];
      real_T y[14400];
      real_T b_y[14400];
    } f1;
  } u1;

  struct {
    real_T Q[14400];
    real_T y[14400];
    real_T b_y[14400];
  } f2;
} logInteriorPoint_rtStackData;

#endif                                 /*typedef_logInteriorPoint_rtStackData*/

/* End of code generation (logInteriorPoint_rt_types.h) */
