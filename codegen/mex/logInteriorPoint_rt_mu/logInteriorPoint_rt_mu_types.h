/*
 * logInteriorPoint_rt_mu_types.h
 *
 * Code generation for function 'logInteriorPoint_rt_mu'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"

/* Type Definitions */
#ifndef typedef_logInteriorPoint_rt_muStackData
#define typedef_logInteriorPoint_rt_muStackData

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
} logInteriorPoint_rt_muStackData;

#endif                                 /*typedef_logInteriorPoint_rt_muStackData*/

/* End of code generation (logInteriorPoint_rt_mu_types.h) */
