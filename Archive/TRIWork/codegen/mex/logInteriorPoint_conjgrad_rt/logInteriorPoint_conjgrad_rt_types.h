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
#ifndef typedef_c_logInteriorPoint_conjgrad_rtS
#define typedef_c_logInteriorPoint_conjgrad_rtS

typedef struct {
  union
  {
    struct {
      real_T dv[14400];
    } f0;

    struct {
      real_T dv[14400];
    } f1;
  } u1;

  struct {
    real_T invW[14400];
    real_T y[14400];
  } f2;
} c_logInteriorPoint_conjgrad_rtS;

#endif                                 /*typedef_c_logInteriorPoint_conjgrad_rtS*/

/* End of code generation (logInteriorPoint_conjgrad_rt_types.h) */
