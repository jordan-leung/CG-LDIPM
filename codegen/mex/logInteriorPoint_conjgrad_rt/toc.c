/*
 * toc.c
 *
 * Code generation for function 'toc'
 *
 */

/* Include files */
#include "toc.h"
#include "logInteriorPoint_conjgrad_rt_data.h"
#include "rt_nonfinite.h"
#include "timeKeeper.h"
#include "emlrt.h"

/* Variable Definitions */
static emlrtRSInfo nc_emlrtRSI = { 31, /* lineNo */
  "toc",                               /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/timefun/toc.m"/* pathName */
};

static emlrtRSInfo oc_emlrtRSI = { 36, /* lineNo */
  "toc",                               /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/timefun/toc.m"/* pathName */
};

/* Function Definitions */
real_T toc(const emlrtStack *sp)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  emlrtTimespec tnow;
  real_T tstart_tv_nsec;
  real_T tstart_tv_sec;
  int32_T status;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &nc_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  b_timeKeeper(&st, &tstart_tv_sec, &tstart_tv_nsec);
  st.site = &oc_emlrtRSI;
  b_st.site = &jb_emlrtRSI;
  c_st.site = &kb_emlrtRSI;
  status = emlrtClockGettimeMonotonic(&tnow);
  d_st.site = &lb_emlrtRSI;
  if (status != 0) {
    emlrtErrorWithMessageIdR2018a(&d_st, &c_emlrtRTEI,
      "Coder:toolbox:POSIXCallFailed", "Coder:toolbox:POSIXCallFailed", 5, 4, 26,
      "emlrtClockGettimeMonotonic", 12, status);
  }

  return (tnow.tv_sec - tstart_tv_sec) + (tnow.tv_nsec - tstart_tv_nsec) /
    1.0E+9;
}

/* End of code generation (toc.c) */
