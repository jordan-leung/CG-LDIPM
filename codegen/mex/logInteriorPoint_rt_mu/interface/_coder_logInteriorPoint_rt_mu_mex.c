/*
 * _coder_logInteriorPoint_rt_mu_mex.c
 *
 * Code generation for function '_coder_logInteriorPoint_rt_mu_mex'
 *
 */

/* Include files */
#include "_coder_logInteriorPoint_rt_mu_mex.h"
#include "_coder_logInteriorPoint_rt_mu_api.h"
#include "logInteriorPoint_rt_mu_data.h"
#include "logInteriorPoint_rt_mu_initialize.h"
#include "logInteriorPoint_rt_mu_terminate.h"
#include "logInteriorPoint_rt_mu_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void logInteriorPoint_rt_mu_mexFunction(logInteriorPoint_rt_muStackData *SD,
  int32_T nlhs, mxArray *plhs[4], int32_T nrhs, const mxArray *prhs[8])
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
                        22, "logInteriorPoint_rt_mu");
  }

  if (nlhs > 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 22,
                        "logInteriorPoint_rt_mu");
  }

  /* Call the function. */
  logInteriorPoint_rt_mu_api(SD, prhs, nlhs, outputs);

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
  logInteriorPoint_rt_muStackData *c_logInteriorPoint_rt_muStackDa = NULL;
  c_logInteriorPoint_rt_muStackDa = (logInteriorPoint_rt_muStackData *)
    emlrtMxCalloc(1, (size_t)1U * sizeof(logInteriorPoint_rt_muStackData));
  mexAtExit(&logInteriorPoint_rt_mu_atexit);

  /* Module initialization. */
  logInteriorPoint_rt_mu_initialize();

  /* Dispatch the entry-point. */
  logInteriorPoint_rt_mu_mexFunction(c_logInteriorPoint_rt_muStackDa, nlhs, plhs,
    nrhs, prhs);

  /* Module termination. */
  logInteriorPoint_rt_mu_terminate();
  emlrtMxFree(c_logInteriorPoint_rt_muStackDa);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_logInteriorPoint_rt_mu_mex.c) */
