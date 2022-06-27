/*
 * _coder_projGradSolver_rt_mex.c
 *
 * Code generation for function '_coder_projGradSolver_rt_mex'
 *
 */

/* Include files */
#include "_coder_projGradSolver_rt_mex.h"
#include "_coder_projGradSolver_rt_api.h"
#include "projGradSolver_rt_data.h"
#include "projGradSolver_rt_initialize.h"
#include "projGradSolver_rt_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(&projGradSolver_rt_atexit);

  /* Module initialization. */
  projGradSolver_rt_initialize();

  /* Dispatch the entry-point. */
  projGradSolver_rt_mexFunction(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  projGradSolver_rt_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

void projGradSolver_rt_mexFunction(int32_T nlhs, mxArray *plhs[4], int32_T nrhs,
  const mxArray *prhs[8])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  const mxArray *outputs[4];
  int32_T b_nlhs;
  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 8) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 8, 4,
                        17, "projGradSolver_rt");
  }

  if (nlhs > 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 17,
                        "projGradSolver_rt");
  }

  /* Call the function. */
  projGradSolver_rt_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);
}

/* End of code generation (_coder_projGradSolver_rt_mex.c) */
