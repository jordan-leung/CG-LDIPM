/*
 * logInteriorPoint_conjgrad_rt.h
 *
 * Code generation for function 'logInteriorPoint_conjgrad_rt'
 *
 */

#pragma once

/* Include files */
#include "logInteriorPoint_conjgrad_rt_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void logInteriorPoint_conjgrad_rt(const emlrtStack *sp, const real_T W[100],
  const real_T c[10], real_T Aineq[200], const real_T bineq[20], real_T mu_f,
  real_T mu_0, const real_T v0[20], real_T maxIter, real_T maxCGIter, real_T
  preCondFlag, const real_T xStar[10], real_T xTol, real_T x[10],
  emxArray_real_T *xError_vec, real_T *execTime, real_T *numIter);

/* End of code generation (logInteriorPoint_conjgrad_rt.h) */
