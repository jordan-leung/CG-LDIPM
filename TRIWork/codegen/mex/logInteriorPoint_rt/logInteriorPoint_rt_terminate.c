/*
 * logInteriorPoint_rt_terminate.c
 *
 * Code generation for function 'logInteriorPoint_rt_terminate'
 *
 */

/* Include files */
#include "logInteriorPoint_rt_terminate.h"
#include "_coder_logInteriorPoint_rt_mex.h"
#include "logInteriorPoint_rt_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void logInteriorPoint_rt_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void logInteriorPoint_rt_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (logInteriorPoint_rt_terminate.c) */
