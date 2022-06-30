/*
 * logInteriorPoint_rt.h
 *
 * Code generation for function 'logInteriorPoint_rt'
 *
 */

#pragma once

/* Include files */
#include "logInteriorPoint_rt_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void logInteriorPoint_rt(logInteriorPoint_rtStackData *SD, const emlrtStack *sp,
  const real_T W[14400], const real_T c[120], real_T A[14400], const real_T
  bineq[120], real_T mu_f, real_T mu_0, const real_T v0[120], real_T maxIter,
  real_T x[120], real_T *mu, real_T *execTime, real_T *numIter);

/* End of code generation (logInteriorPoint_rt.h) */
