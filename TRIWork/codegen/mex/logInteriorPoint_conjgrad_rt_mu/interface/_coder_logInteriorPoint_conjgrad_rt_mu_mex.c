/*
 * _coder_logInteriorPoint_conjgrad_rt_mu_mex.c
 *
 * Code generation for function '_coder_logInteriorPoint_conjgrad_rt_mu_mex'
 *
 */

/* Include files */
#include "_coder_logInteriorPoint_conjgrad_rt_mu_mex.h"
#include "_coder_logInteriorPoint_conjgrad_rt_mu_api.h"
#include "logInteriorPoint_conjgrad_rt_mu_data.h"
#include "logInteriorPoint_conjgrad_rt_mu_initialize.h"
#include "logInteriorPoint_conjgrad_rt_mu_terminate.h"
#include "logInteriorPoint_conjgrad_rt_mu_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void logInteriorPoint_conjgrad_rt_mu_mexFunction(e_logInteriorPoint_conjgrad_rt_
  *SD, int32_T nlhs, mxArray *plhs[5], int32_T nrhs, const mxArray *prhs[10])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  const mxArray *outputs[5];
  int32_T b_nlhs;
  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 10) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 10, 4,
                        31, "logInteriorPoint_conjgrad_rt_mu");
  }

  if (nlhs > 5) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 31,
                        "logInteriorPoint_conjgrad_rt_mu");
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
  e_logInteriorPoint_conjgrad_rt_ *f_logInteriorPoint_conjgrad_rt_ = NULL;
  f_logInteriorPoint_conjgrad_rt_ = (e_logInteriorPoint_conjgrad_rt_ *)
    emlrtMxCalloc(1, (size_t)1U * sizeof(e_logInteriorPoint_conjgrad_rt_));
  mexAtExit(&logInteriorPoint_conjgrad_rt_mu_atexit);

  /* Module initialization. */
  logInteriorPoint_conjgrad_rt_mu_initialize();

  /* Dispatch the entry-point. */
  logInteriorPoint_conjgrad_rt_mu_mexFunction(f_logInteriorPoint_conjgrad_rt_,
    nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  logInteriorPoint_conjgrad_rt_mu_terminate();
  emlrtMxFree(f_logInteriorPoint_conjgrad_rt_);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_logInteriorPoint_conjgrad_rt_mu_mex.c) */
