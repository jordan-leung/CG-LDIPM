/*
 * _coder_logInteriorPoint_conjgrad_rt_mu_mex.h
 *
 * Code generation for function '_coder_logInteriorPoint_conjgrad_rt_mu_mex'
 *
 */

#pragma once

/* Include files */
#include "logInteriorPoint_conjgrad_rt_mu_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void logInteriorPoint_conjgrad_rt_mu_mexFunction(e_logInteriorPoint_conjgrad_rt_
  *SD, int32_T nlhs, mxArray *plhs[5], int32_T nrhs, const mxArray *prhs[10]);
MEXFUNCTION_LINKAGE void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
  const mxArray *prhs[]);
emlrtCTX mexFunctionCreateRootTLS(void);

/* End of code generation (_coder_logInteriorPoint_conjgrad_rt_mu_mex.h) */
