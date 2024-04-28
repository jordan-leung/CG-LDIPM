/*
 * mldivide.c
 *
 * Code generation for function 'mldivide'
 *
 */

/* Include files */
#include "mldivide.h"
#include "logInteriorPoint_rt_mu_data.h"
#include "rt_nonfinite.h"
#include "warning.h"
#include "blas.h"
#include "lapacke.h"
#include <stddef.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo ac_emlrtRSI = { 20, /* lineNo */
  "mldivide",                          /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/ops/mldivide.m"/* pathName */
};

static emlrtRSInfo bc_emlrtRSI = { 42, /* lineNo */
  "mldiv",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/ops/mldivide.m"/* pathName */
};

static emlrtRSInfo cc_emlrtRSI = { 67, /* lineNo */
  "lusolve",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/lusolve.m"/* pathName */
};

static emlrtRSInfo dc_emlrtRSI = { 109,/* lineNo */
  "lusolveNxN",                        /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/lusolve.m"/* pathName */
};

static emlrtRSInfo ec_emlrtRSI = { 112,/* lineNo */
  "lusolveNxN",                        /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/lusolve.m"/* pathName */
};

static emlrtRSInfo fc_emlrtRSI = { 124,/* lineNo */
  "InvAtimesX",                        /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/lusolve.m"/* pathName */
};

static emlrtRSInfo gc_emlrtRSI = { 19, /* lineNo */
  "xgetrfs",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgetrfs.m"/* pathName */
};

static emlrtRSInfo hc_emlrtRSI = { 108,/* lineNo */
  "cmldiv",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgetrfs.m"/* pathName */
};

static emlrtRSInfo ic_emlrtRSI = { 90, /* lineNo */
  "warn_singular",                     /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/lusolve.m"/* pathName */
};

/* Function Definitions */
void mldivide(const emlrtStack *sp, const real_T A[14400], real_T B[120])
{
  ptrdiff_t IPIV[120];
  ptrdiff_t INFO;
  ptrdiff_t LDA;
  ptrdiff_t N;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack st;
  real_T b_A[14400];
  int32_T info;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &ac_emlrtRSI;
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
  b_st.site = &bc_emlrtRSI;
  c_st.site = &cc_emlrtRSI;
  d_st.site = &dc_emlrtRSI;
  e_st.site = &fc_emlrtRSI;
  f_st.site = &gc_emlrtRSI;
  memcpy(&b_A[0], &A[0], 14400U * sizeof(real_T));
  N = (ptrdiff_t)120;
  LDA = (ptrdiff_t)120;
  INFO = LAPACKE_dgetrf_work(102, N, N, &b_A[0], LDA, &IPIV[0]);
  info = (int32_T)INFO;
  g_st.site = &hc_emlrtRSI;
  if (info < 0) {
    if (info == -1010) {
      emlrtErrorWithMessageIdR2018a(&g_st, &c_emlrtRTEI, "MATLAB:nomem",
        "MATLAB:nomem", 0);
    } else {
      emlrtErrorWithMessageIdR2018a(&g_st, &b_emlrtRTEI,
        "Coder:toolbox:LAPACKCallErrorInfo", "Coder:toolbox:LAPACKCallErrorInfo",
        5, 4, 19, cv, 12, info);
    }
  }

  LAPACKE_dgetrs_work(102, 'N', N, (ptrdiff_t)1, &b_A[0], LDA, &IPIV[0], &B[0],
                      (ptrdiff_t)120);
  if (info > 0) {
    d_st.site = &ec_emlrtRSI;
    e_st.site = &ic_emlrtRSI;
    warning(&e_st);
  }
}

/* End of code generation (mldivide.c) */
