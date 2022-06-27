/*
 * mldivide.c
 *
 * Code generation for function 'mldivide'
 *
 */

/* Include files */
#include "mldivide.h"
#include "logInteriorPoint_rt_data.h"
#include "rt_nonfinite.h"
#include "warning.h"
#include "xzgetrf.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo fc_emlrtRSI = { 20, /* lineNo */
  "mldivide",                          /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/ops/mldivide.m"/* pathName */
};

static emlrtRSInfo gc_emlrtRSI = { 42, /* lineNo */
  "mldiv",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/ops/mldivide.m"/* pathName */
};

static emlrtRSInfo hc_emlrtRSI = { 67, /* lineNo */
  "lusolve",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/lusolve.m"/* pathName */
};

static emlrtRSInfo ic_emlrtRSI = { 109,/* lineNo */
  "lusolveNxN",                        /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/lusolve.m"/* pathName */
};

static emlrtRSInfo jc_emlrtRSI = { 112,/* lineNo */
  "lusolveNxN",                        /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/lusolve.m"/* pathName */
};

static emlrtRSInfo kc_emlrtRSI = { 124,/* lineNo */
  "InvAtimesX",                        /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/lusolve.m"/* pathName */
};

static emlrtRSInfo lc_emlrtRSI = { 26, /* lineNo */
  "xgetrfs",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgetrfs.m"/* pathName */
};

static emlrtRSInfo mc_emlrtRSI = { 90, /* lineNo */
  "warn_singular",                     /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/lusolve.m"/* pathName */
};

/* Function Definitions */
void mldivide(const emlrtStack *sp, const real_T A[100], real_T B[10])
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack st;
  real_T b_A[100];
  real_T temp;
  int32_T ipiv[10];
  int32_T b_i;
  int32_T i;
  int32_T info;
  int32_T k;
  int32_T kAcol;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &fc_emlrtRSI;
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
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  b_st.site = &gc_emlrtRSI;
  c_st.site = &hc_emlrtRSI;
  d_st.site = &ic_emlrtRSI;
  e_st.site = &kc_emlrtRSI;
  f_st.site = &lc_emlrtRSI;
  memcpy(&b_A[0], &A[0], 100U * sizeof(real_T));
  g_st.site = &db_emlrtRSI;
  xzgetrf(&g_st, b_A, ipiv, &info);
  for (i = 0; i < 9; i++) {
    b_i = ipiv[i];
    if (b_i != i + 1) {
      temp = B[i];
      B[i] = B[b_i - 1];
      B[b_i - 1] = temp;
    }
  }

  for (k = 0; k < 10; k++) {
    kAcol = 10 * k;
    if (B[k] != 0.0) {
      b_i = k + 2;
      for (i = b_i; i < 11; i++) {
        B[i - 1] -= B[k] * b_A[(i + kAcol) - 1];
      }
    }
  }

  for (k = 9; k >= 0; k--) {
    kAcol = 10 * k;
    if (B[k] != 0.0) {
      B[k] /= b_A[k + kAcol];
      for (i = 0; i < k; i++) {
        B[i] -= B[k] * b_A[i + kAcol];
      }
    }
  }

  if (info > 0) {
    d_st.site = &jc_emlrtRSI;
    e_st.site = &mc_emlrtRSI;
    warning(&e_st);
  }
}

/* End of code generation (mldivide.c) */
