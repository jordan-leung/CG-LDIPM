/*
 * logInteriorPoint_conjgrad_rt_mu_types.h
 *
 * Code generation for function 'logInteriorPoint_conjgrad_rt_mu'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"

/* Type Definitions */
#ifndef typedef_e_logInteriorPoint_conjgrad_rt_
#define typedef_e_logInteriorPoint_conjgrad_rt_

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
} e_logInteriorPoint_conjgrad_rt_;

#endif                                 /*typedef_e_logInteriorPoint_conjgrad_rt_*/

/* End of code generation (logInteriorPoint_conjgrad_rt_mu_types.h) */
