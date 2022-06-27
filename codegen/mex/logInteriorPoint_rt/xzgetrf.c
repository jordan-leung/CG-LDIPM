/*
 * xzgetrf.c
 *
 * Code generation for function 'xzgetrf'
 *
 */

/* Include files */
#include "xzgetrf.h"
#include "eml_int_forloop_overflow_check.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo eb_emlrtRSI = { 50, /* lineNo */
  "xzgetrf",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzgetrf.m"/* pathName */
};

static emlrtRSInfo fb_emlrtRSI = { 58, /* lineNo */
  "xzgetrf",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzgetrf.m"/* pathName */
};

static emlrtRSInfo gb_emlrtRSI = { 21, /* lineNo */
  "eml_int_forloop_overflow_check",    /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"/* pathName */
};

static emlrtRSInfo hb_emlrtRSI = { 45, /* lineNo */
  "xgeru",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+blas/xgeru.m"/* pathName */
};

static emlrtRSInfo ib_emlrtRSI = { 45, /* lineNo */
  "xger",                              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+blas/xger.m"/* pathName */
};

static emlrtRSInfo jb_emlrtRSI = { 15, /* lineNo */
  "xger",                              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+refblas/xger.m"/* pathName */
};

static emlrtRSInfo kb_emlrtRSI = { 54, /* lineNo */
  "xgerx",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+refblas/xgerx.m"/* pathName */
};

/* Function Definitions */
void xzgetrf(const emlrtStack *sp, real_T A[100], int32_T ipiv[10], int32_T
             *info)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  real_T s;
  real_T smax;
  int32_T b_j;
  int32_T b_tmp;
  int32_T ijA;
  int32_T ix;
  int32_T iy;
  int32_T j;
  int32_T jp1j;
  int32_T k;
  int32_T n;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  for (iy = 0; iy < 10; iy++) {
    ipiv[iy] = iy + 1;
  }

  *info = 0;
  for (j = 0; j < 9; j++) {
    b_tmp = j * 11;
    jp1j = b_tmp + 2;
    n = 10 - j;
    iy = 0;
    ix = b_tmp;
    smax = muDoubleScalarAbs(A[b_tmp]);
    for (k = 2; k <= n; k++) {
      ix++;
      s = muDoubleScalarAbs(A[ix]);
      if (s > smax) {
        iy = k - 1;
        smax = s;
      }
    }

    if (A[b_tmp + iy] != 0.0) {
      if (iy != 0) {
        iy += j;
        ipiv[j] = iy + 1;
        ix = j;
        for (k = 0; k < 10; k++) {
          smax = A[ix];
          A[ix] = A[iy];
          A[iy] = smax;
          ix += 10;
          iy += 10;
        }
      }

      k = (b_tmp - j) + 10;
      st.site = &eb_emlrtRSI;
      for (iy = jp1j; iy <= k; iy++) {
        A[iy - 1] /= A[b_tmp];
      }
    } else {
      *info = j + 1;
    }

    n = 8 - j;
    iy = b_tmp + 10;
    st.site = &fb_emlrtRSI;
    b_st.site = &hb_emlrtRSI;
    c_st.site = &ib_emlrtRSI;
    d_st.site = &jb_emlrtRSI;
    jp1j = b_tmp + 12;
    for (b_j = 0; b_j <= n; b_j++) {
      smax = A[iy];
      if (A[iy] != 0.0) {
        ix = b_tmp + 1;
        k = (jp1j - j) + 8;
        e_st.site = &kb_emlrtRSI;
        if ((jp1j <= k) && (k > 2147483646)) {
          f_st.site = &gb_emlrtRSI;
          check_forloop_overflow_error(&f_st);
        }

        for (ijA = jp1j; ijA <= k; ijA++) {
          A[ijA - 1] += A[ix] * -smax;
          ix++;
        }
      }

      iy += 10;
      jp1j += 10;
    }
  }

  if ((*info == 0) && (!(A[99] != 0.0))) {
    *info = 10;
  }
}

/* End of code generation (xzgetrf.c) */
