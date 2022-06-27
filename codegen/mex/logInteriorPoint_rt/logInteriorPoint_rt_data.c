/*
 * logInteriorPoint_rt_data.c
 *
 * Code generation for function 'logInteriorPoint_rt_data'
 *
 */

/* Include files */
#include "logInteriorPoint_rt_data.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131595U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "logInteriorPoint_rt",               /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

emlrtRSInfo db_emlrtRSI = { 30,        /* lineNo */
  "xgetrf",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgetrf.m"/* pathName */
};

emlrtRSInfo qb_emlrtRSI = { 9,         /* lineNo */
  "getTime",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/shared/coder/coder/lib/+coder/+internal/+time/getTime.m"/* pathName */
};

emlrtRSInfo rb_emlrtRSI = { 19,        /* lineNo */
  "callEMLRTClockGettime",             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/shared/coder/coder/lib/+coder/+internal/+time/private/callEMLRTClockGettime.m"/* pathName */
};

emlrtRSInfo sb_emlrtRSI = { 29,        /* lineNo */
  "getTimeEMLRT",                      /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/shared/coder/coder/lib/+coder/+internal/+time/private/callEMLRTClockGettime.m"/* pathName */
};

emlrtRTEInfo b_emlrtRTEI = { 37,       /* lineNo */
  9,                                   /* colNo */
  "checkPOSIXStatus",                  /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/shared/coder/coder/lib/+coder/+internal/+time/private/callEMLRTClockGettime.m"/* pName */
};

/* End of code generation (logInteriorPoint_rt_data.c) */
