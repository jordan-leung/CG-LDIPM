/*
 * logInteriorPoint_conjgrad_rt_mu_initialize.c
 *
 * Code generation for function 'logInteriorPoint_conjgrad_rt_mu_initialize'
 *
 */

/* Include files */
#include "logInteriorPoint_conjgrad_rt_mu_initialize.h"
#include "_coder_logInteriorPoint_conjgrad_rt_mu_mex.h"
#include "logInteriorPoint_conjgrad_rt_mu_data.h"
#include "rt_nonfinite.h"
#include "timeKeeper.h"

/* Function Declarations */
static void c_logInteriorPoint_conjgrad_rt_(void);

/* Function Definitions */
static void c_logInteriorPoint_conjgrad_rt_(void)
{
  mex_InitInfAndNan();
  savedTime_not_empty_init();
}

void logInteriorPoint_conjgrad_rt_mu_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    c_logInteriorPoint_conjgrad_rt_();
  }
}

/* End of code generation (logInteriorPoint_conjgrad_rt_mu_initialize.c) */
