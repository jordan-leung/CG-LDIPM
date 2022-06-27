/*
 * projGradSolver_rt_emxutil.h
 *
 * Code generation for function 'projGradSolver_rt_emxutil'
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
void emxEnsureCapacity_real_T(const emlrtStack *sp, emxArray_real_T *emxArray,
  int32_T oldNumel, const emlrtRTEInfo *srcLocation);
void emxFree_real_T(emxArray_real_T **pEmxArray);
void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray, int32_T
                    numDimensions, const emlrtRTEInfo *srcLocation, boolean_T
                    doPush);

/* End of code generation (projGradSolver_rt_emxutil.h) */
