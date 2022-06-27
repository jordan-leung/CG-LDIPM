/*
 * projGradSolver_rt.h
 *
 * Code generation for function 'projGradSolver_rt'
 *
 */

#pragma once

/* Include files */
#include "projGradSolver_rt_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void projGradSolver_rt(const emlrtStack *sp, const real_T H[100], const real_T
  f[10], const real_T x0[10], const real_T xl[10], const real_T xu[10], real_T
  MaxIter, const real_T xOpt[10], real_T xTol, real_T x[10], real_T *iterCount,
  emxArray_real_T *xError_vec, real_T *execTime);

/* End of code generation (projGradSolver_rt.h) */
