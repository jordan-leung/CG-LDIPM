/*
 * projGradSolver_rt_data.c
 *
 * Code generation for function 'projGradSolver_rt_data'
 *
 */

/* Include files */
#include "projGradSolver_rt_data.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131595U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "projGradSolver_rt",                 /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

emlrtRSInfo m_emlrtRSI = { 21,         /* lineNo */
  "eml_int_forloop_overflow_check",    /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"/* pathName */
};

emlrtRSInfo bb_emlrtRSI = { 71,        /* lineNo */
  "xzlarfg",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzlarfg.m"/* pathName */
};

emlrtRSInfo jb_emlrtRSI = { 64,        /* lineNo */
  "xgemv",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+blas/xgemv.m"/* pathName */
};

emlrtRSInfo lb_emlrtRSI = { 37,        /* lineNo */
  "xgemv",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+refblas/xgemv.m"/* pathName */
};

emlrtRSInfo mb_emlrtRSI = { 45,        /* lineNo */
  "xgerc",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+blas/xgerc.m"/* pathName */
};

emlrtRSInfo nb_emlrtRSI = { 45,        /* lineNo */
  "xger",                              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+blas/xger.m"/* pathName */
};

emlrtRSInfo ob_emlrtRSI = { 15,        /* lineNo */
  "xger",                              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+refblas/xger.m"/* pathName */
};

emlrtRSInfo pb_emlrtRSI = { 41,        /* lineNo */
  "xgerx",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+refblas/xgerx.m"/* pathName */
};

emlrtRSInfo qb_emlrtRSI = { 54,        /* lineNo */
  "xgerx",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+refblas/xgerx.m"/* pathName */
};

emlrtRSInfo pd_emlrtRSI = { 9,         /* lineNo */
  "getTime",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/shared/coder/coder/lib/+coder/+internal/+time/getTime.m"/* pathName */
};

emlrtRSInfo qd_emlrtRSI = { 19,        /* lineNo */
  "callEMLRTClockGettime",             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/shared/coder/coder/lib/+coder/+internal/+time/private/callEMLRTClockGettime.m"/* pathName */
};

emlrtRSInfo rd_emlrtRSI = { 29,        /* lineNo */
  "getTimeEMLRT",                      /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/shared/coder/coder/lib/+coder/+internal/+time/private/callEMLRTClockGettime.m"/* pathName */
};

emlrtRTEInfo d_emlrtRTEI = { 37,       /* lineNo */
  9,                                   /* colNo */
  "checkPOSIXStatus",                  /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/shared/coder/coder/lib/+coder/+internal/+time/private/callEMLRTClockGettime.m"/* pName */
};

/* End of code generation (projGradSolver_rt_data.c) */
