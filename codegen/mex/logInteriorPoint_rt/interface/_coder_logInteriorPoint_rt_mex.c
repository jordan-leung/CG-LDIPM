/*
 * _coder_logInteriorPoint_rt_mex.c
 *
 * Code generation for function '_coder_logInteriorPoint_rt_mex'
 *
 */

/* Include files */
#include "_coder_logInteriorPoint_rt_mex.h"
#include "_coder_logInteriorPoint_rt_api.h"
#include "logInteriorPoint_rt_data.h"
#include "logInteriorPoint_rt_initialize.h"
#include "logInteriorPoint_rt_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void logInteriorPoint_rt_mexFunction(int32_T nlhs, mxArray *plhs[4], int32_T
  nrhs, const mxArray *prhs[10])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  const mxArray *outputs[4];
  int32_T b_nlhs;
  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 10) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 10, 4,
                        19, "logInteriorPoint_rt");
  }

  if (nlhs > 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 19,
                        "logInteriorPoint_rt");
  }

  /* Call the function. */
  logInteriorPoint_rt_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(&logInteriorPoint_rt_atexit);

  /* Module initialization. */
  logInteriorPoint_rt_initialize();

  /* Dispatch the entry-point. */
  logInteriorPoint_rt_mexFunction(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  logInteriorPoint_rt_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_logInteriorPoint_rt_mex.c) */
