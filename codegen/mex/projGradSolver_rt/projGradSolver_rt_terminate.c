/*
 * projGradSolver_rt_terminate.c
 *
 * Code generation for function 'projGradSolver_rt_terminate'
 *
 */

/* Include files */
#include "projGradSolver_rt_terminate.h"
#include "_coder_projGradSolver_rt_mex.h"
#include "projGradSolver_rt_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void projGradSolver_rt_atexit(void)
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

void projGradSolver_rt_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (projGradSolver_rt_terminate.c) */
