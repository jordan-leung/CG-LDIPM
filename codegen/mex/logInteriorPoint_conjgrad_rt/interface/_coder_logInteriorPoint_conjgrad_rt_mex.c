/*
 * _coder_logInteriorPoint_conjgrad_rt_mex.c
 *
 * Code generation for function '_coder_logInteriorPoint_conjgrad_rt_mex'
 *
 */

/* Include files */
#include "_coder_logInteriorPoint_conjgrad_rt_mex.h"
#include "_coder_logInteriorPoint_conjgrad_rt_api.h"
#include "logInteriorPoint_conjgrad_rt_data.h"
#include "logInteriorPoint_conjgrad_rt_initialize.h"
#include "logInteriorPoint_conjgrad_rt_terminate.h"
#include "logInteriorPoint_conjgrad_rt_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void logInteriorPoint_conjgrad_rt_mexFunction(c_logInteriorPoint_conjgrad_rtS
  *SD, int32_T nlhs, mxArray *plhs[5], int32_T nrhs, const mxArray *prhs[11])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  const mxArray *outputs[5];
  int32_T b_nlhs;
  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 11) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 11, 4,
                        28, "logInteriorPoint_conjgrad_rt");
  }

  if (nlhs > 5) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 28,
                        "logInteriorPoint_conjgrad_rt");
  }

  /* Call the function. */
  d_logInteriorPoint_conjgrad_rt_(SD, prhs, nlhs, outputs);

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
  c_logInteriorPoint_conjgrad_rtS *d_logInteriorPoint_conjgrad_rtS = NULL;
  d_logInteriorPoint_conjgrad_rtS = (c_logInteriorPoint_conjgrad_rtS *)
    emlrtMxCalloc(1, (size_t)1U * sizeof(c_logInteriorPoint_conjgrad_rtS));
  mexAtExit(&logInteriorPoint_conjgrad_rt_atexit);

  /* Module initialization. */
  logInteriorPoint_conjgrad_rt_initialize();

  /* Dispatch the entry-point. */
  logInteriorPoint_conjgrad_rt_mexFunction(d_logInteriorPoint_conjgrad_rtS, nlhs,
    plhs, nrhs, prhs);

  /* Module termination. */
  logInteriorPoint_conjgrad_rt_terminate();
  emlrtMxFree(d_logInteriorPoint_conjgrad_rtS);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_logInteriorPoint_conjgrad_rt_mex.c) */
