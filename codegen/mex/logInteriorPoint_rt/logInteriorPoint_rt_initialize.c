/*
 * logInteriorPoint_rt_initialize.c
 *
 * Code generation for function 'logInteriorPoint_rt_initialize'
 *
 */

/* Include files */
#include "logInteriorPoint_rt_initialize.h"
#include "_coder_logInteriorPoint_rt_mex.h"
#include "logInteriorPoint_rt_data.h"
#include "rt_nonfinite.h"
#include "timeKeeper.h"

/* Function Declarations */
static void logInteriorPoint_rt_once(void);

/* Function Definitions */
static void logInteriorPoint_rt_once(void)
{
  mex_InitInfAndNan();
  savedTime_not_empty_init();
}

void logInteriorPoint_rt_initialize(void)
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
    logInteriorPoint_rt_once();
  }
}

/* End of code generation (logInteriorPoint_rt_initialize.c) */
