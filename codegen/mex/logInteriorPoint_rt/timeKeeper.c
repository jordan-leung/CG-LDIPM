/*
 * timeKeeper.c
 *
 * Code generation for function 'timeKeeper'
 *
 */

/* Include files */
#include "timeKeeper.h"
#include "logInteriorPoint_rt_data.h"
#include "rt_nonfinite.h"
#include "emlrt.h"

/* Variable Definitions */
static emlrtTimespec savedTime;
static boolean_T savedTime_not_empty;
static emlrtRSInfo tb_emlrtRSI = { 13, /* lineNo */
  "timeKeeper",                        /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/shared/coder/coder/lib/+coder/+internal/+time/+impl/timeKeeper.m"/* pathName */
};

static emlrtRTEInfo d_emlrtRTEI = { 11,/* lineNo */
  9,                                   /* colNo */
  "timeKeeper",                        /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/shared/coder/coder/lib/+coder/+internal/+time/+impl/timeKeeper.m"/* pName */
};

/* Function Definitions */
void b_timeKeeper(const emlrtStack *sp, real_T *outTime_tv_sec, real_T
                  *outTime_tv_nsec)
{
  if (!savedTime_not_empty) {
    emlrtErrorWithMessageIdR2018a(sp, &d_emlrtRTEI,
      "MATLAB:toc:callTicFirstNoInputs", "MATLAB:toc:callTicFirstNoInputs", 0);
  }

  *outTime_tv_sec = savedTime.tv_sec;
  *outTime_tv_nsec = savedTime.tv_nsec;
}

void savedTime_not_empty_init(void)
{
  savedTime_not_empty = false;
}

void timeKeeper(const emlrtStack *sp, const emlrtTimespec newTime)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  int32_T status;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  if (!savedTime_not_empty) {
    st.site = &tb_emlrtRSI;
    b_st.site = &qb_emlrtRSI;
    c_st.site = &rb_emlrtRSI;
    status = emlrtClockGettimeMonotonic(&savedTime);
    d_st.site = &sb_emlrtRSI;
    if (status != 0) {
      emlrtErrorWithMessageIdR2018a(&d_st, &b_emlrtRTEI,
        "Coder:toolbox:POSIXCallFailed", "Coder:toolbox:POSIXCallFailed", 5, 4,
        26, "emlrtClockGettimeMonotonic", 12, status);
    }

    savedTime_not_empty = true;
  }

  savedTime = newTime;
}

/* End of code generation (timeKeeper.c) */
