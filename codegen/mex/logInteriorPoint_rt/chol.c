/*
 * chol.c
 *
 * Code generation for function 'chol'
 *
 */

/* Include files */
#include "chol.h"
#include "logInteriorPoint_rt_data.h"
#include "rt_nonfinite.h"
#include "lapacke.h"
#include <stddef.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo sb_emlrtRSI = { 74, /* lineNo */
  "cholesky",                          /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/chol.m"/* pathName */
};

static emlrtRSInfo tb_emlrtRSI = { 13, /* lineNo */
  "xpotrf",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/xpotrf.m"/* pathName */
};

static emlrtRSInfo ub_emlrtRSI = { 79, /* lineNo */
  "ceval_xpotrf",                      /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/xpotrf.m"/* pathName */
};

static emlrtRTEInfo g_emlrtRTEI = { 80,/* lineNo */
  23,                                  /* colNo */
  "cholesky",                          /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/chol.m"/* pName */
};

/* Function Definitions */
void cholesky(const emlrtStack *sp, real_T A[10000])
{
  static const char_T fname[19] = { 'L', 'A', 'P', 'A', 'C', 'K', 'E', '_', 'd',
    'p', 'o', 't', 'r', 'f', '_', 'w', 'o', 'r', 'k' };

  ptrdiff_t info_t;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  int32_T info;
  int32_T j;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &sb_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  b_st.site = &tb_emlrtRSI;
  info_t = LAPACKE_dpotrf_work(102, 'U', (ptrdiff_t)100, &A[0], (ptrdiff_t)100);
  info = (int32_T)info_t;
  c_st.site = &ub_emlrtRSI;
  if (info < 0) {
    if (info == -1010) {
      emlrtErrorWithMessageIdR2018a(&c_st, &emlrtRTEI, "MATLAB:nomem",
        "MATLAB:nomem", 0);
    } else {
      emlrtErrorWithMessageIdR2018a(&c_st, &b_emlrtRTEI,
        "Coder:toolbox:LAPACKCallErrorInfo", "Coder:toolbox:LAPACKCallErrorInfo",
        5, 4, 19, fname, 12, info);
    }
  }

  if (info != 0) {
    emlrtErrorWithMessageIdR2018a(sp, &g_emlrtRTEI, "Coder:MATLAB:posdef",
      "Coder:MATLAB:posdef", 0);
  }

  for (j = 0; j < 100; j++) {
    info = j + 2;
    if (info <= 100) {
      memset(&A[(j * 100 + info) + -1], 0, (101 - info) * sizeof(real_T));
    }
  }
}

/* End of code generation (chol.c) */
