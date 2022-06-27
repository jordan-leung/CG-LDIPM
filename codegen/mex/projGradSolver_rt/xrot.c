/*
 * xrot.c
 *
 * Code generation for function 'xrot'
 *
 */

/* Include files */
#include "xrot.h"
#include "eml_int_forloop_overflow_check.h"
#include "projGradSolver_rt_data.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo oc_emlrtRSI = { 32, /* lineNo */
  "xrot",                              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+blas/xrot.m"/* pathName */
};

static emlrtRSInfo pc_emlrtRSI = { 24, /* lineNo */
  "xrot",                              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+refblas/xrot.m"/* pathName */
};

/* Function Definitions */
void b_xrot(const emlrtStack *sp, int32_T n, real_T x[100], int32_T ix0, int32_T
            iy0, real_T c, real_T s)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  real_T temp;
  int32_T ix;
  int32_T iy;
  int32_T k;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &oc_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  if (n >= 1) {
    ix = ix0 - 1;
    iy = iy0 - 1;
    b_st.site = &pc_emlrtRSI;
    if (n > 2147483646) {
      c_st.site = &m_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }

    for (k = 0; k < n; k++) {
      temp = c * x[ix] + s * x[iy];
      x[iy] = c * x[iy] - s * x[ix];
      x[ix] = temp;
      iy++;
      ix++;
    }
  }
}

void xrot(const emlrtStack *sp, int32_T n, real_T x[100], int32_T ix0, int32_T
          iy0, real_T c, real_T s)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  real_T temp;
  int32_T ix;
  int32_T iy;
  int32_T k;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &oc_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  ix = ix0 - 1;
  iy = iy0 - 1;
  b_st.site = &pc_emlrtRSI;
  if ((1 <= n) && (n > 2147483646)) {
    c_st.site = &m_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  for (k = 0; k < n; k++) {
    temp = c * x[ix] + s * x[iy];
    x[iy] = c * x[iy] - s * x[ix];
    x[ix] = temp;
    iy += 10;
    ix += 10;
  }
}

/* End of code generation (xrot.c) */
