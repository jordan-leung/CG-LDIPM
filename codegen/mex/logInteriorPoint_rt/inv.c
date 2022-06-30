/*
 * inv.c
 *
 * Code generation for function 'inv'
 *
 */

/* Include files */
#include "inv.h"
#include "logInteriorPoint_rt_data.h"
#include "logInteriorPoint_rt_mexutil.h"
#include "rt_nonfinite.h"
#include "warning.h"
#include "blas.h"
#include "lapacke.h"
#include "mwmathutil.h"
#include <stddef.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo w_emlrtRSI = { 21,  /* lineNo */
  "inv",                               /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo x_emlrtRSI = { 22,  /* lineNo */
  "inv",                               /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo y_emlrtRSI = { 173, /* lineNo */
  "invNxN",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo ab_emlrtRSI = { 180,/* lineNo */
  "invNxN",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo bb_emlrtRSI = { 183,/* lineNo */
  "invNxN",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo cb_emlrtRSI = { 190,/* lineNo */
  "invNxN",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo db_emlrtRSI = { 27, /* lineNo */
  "xgetrf",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgetrf.m"/* pathName */
};

static emlrtRSInfo eb_emlrtRSI = { 91, /* lineNo */
  "ceval_xgetrf",                      /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgetrf.m"/* pathName */
};

static emlrtRSInfo fb_emlrtRSI = { 67, /* lineNo */
  "xtrsm",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+blas/xtrsm.m"/* pathName */
};

static emlrtRSInfo hb_emlrtRSI = { 42, /* lineNo */
  "checkcond",                         /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo ib_emlrtRSI = { 46, /* lineNo */
  "checkcond",                         /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtMCInfo c_emlrtMCI = { 53,  /* lineNo */
  19,                                  /* colNo */
  "flt2str",                           /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/flt2str.m"/* pName */
};

static emlrtRSInfo sc_emlrtRSI = { 53, /* lineNo */
  "flt2str",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/flt2str.m"/* pathName */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[14]);
static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, emlrtMCInfo *location);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *a__output_of_sprintf_, const char_T *identifier, char_T y[14]);
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[14]);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[14])
{
  i_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  const mxArray *m;
  pArrays[0] = b;
  pArrays[1] = c;
  return emlrtCallMATLABR2012b(sp, 1, &m, 2, pArrays, "sprintf", true, location);
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *a__output_of_sprintf_, const char_T *identifier, char_T y[14])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(a__output_of_sprintf_), &thisId, y);
  emlrtDestroyArray(&a__output_of_sprintf_);
}

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[14])
{
  static const int32_T dims[2] = { 1, 14 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", false, 2U, dims);
  emlrtImportCharArrayR2015b(sp, src, &ret[0], 14);
  emlrtDestroyArray(&src);
}

void inv(const emlrtStack *sp, const real_T x[14400], real_T y[14400])
{
  static const int32_T iv[2] = { 1, 6 };

  static const char_T rfmt[6] = { '%', '1', '4', '.', '6', 'e' };

  ptrdiff_t ipiv_t[120];
  ptrdiff_t info_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t n_t;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *m;
  real_T b_x[14400];
  real_T n1x;
  real_T n1xinv;
  real_T s;
  int32_T ipiv[120];
  int32_T a;
  int32_T b_i;
  int32_T i;
  int32_T j;
  int32_T k;
  int32_T pipk;
  int32_T y_tmp;
  char_T str[14];
  char_T DIAGA1;
  char_T SIDE1;
  char_T TRANSA1;
  char_T UPLO1;
  int8_T p[120];
  boolean_T exitg1;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &w_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  memset(&y[0], 0, 14400U * sizeof(real_T));
  b_st.site = &y_emlrtRSI;
  memcpy(&b_x[0], &x[0], 14400U * sizeof(real_T));
  c_st.site = &db_emlrtRSI;
  info_t = LAPACKE_dgetrf_work(102, (ptrdiff_t)120, (ptrdiff_t)120, &b_x[0],
    (ptrdiff_t)120, &ipiv_t[0]);
  d_st.site = &eb_emlrtRSI;
  pipk = (int32_T)info_t;
  if (pipk < 0) {
    if (pipk == -1010) {
      emlrtErrorWithMessageIdR2018a(&d_st, &c_emlrtRTEI, "MATLAB:nomem",
        "MATLAB:nomem", 0);
    } else {
      emlrtErrorWithMessageIdR2018a(&d_st, &b_emlrtRTEI,
        "Coder:toolbox:LAPACKCallErrorInfo", "Coder:toolbox:LAPACKCallErrorInfo",
        5, 4, 19, cv, 12, pipk);
    }
  }

  for (k = 0; k < 120; k++) {
    ipiv[k] = (int32_T)ipiv_t[k];
  }

  for (i = 0; i < 120; i++) {
    p[i] = (int8_T)(i + 1);
  }

  for (k = 0; k < 119; k++) {
    i = ipiv[k];
    if (i > k + 1) {
      pipk = p[i - 1];
      p[i - 1] = p[k];
      p[k] = (int8_T)pipk;
    }
  }

  for (k = 0; k < 120; k++) {
    pipk = 120 * (p[k] - 1);
    y[k + pipk] = 1.0;
    b_st.site = &ab_emlrtRSI;
    for (j = k + 1; j < 121; j++) {
      i = (j + pipk) - 1;
      if (y[i] != 0.0) {
        a = j + 1;
        b_st.site = &bb_emlrtRSI;
        for (b_i = a; b_i < 121; b_i++) {
          y_tmp = (b_i + pipk) - 1;
          y[y_tmp] -= y[i] * b_x[(b_i + 120 * (j - 1)) - 1];
        }
      }
    }
  }

  b_st.site = &cb_emlrtRSI;
  c_st.site = &fb_emlrtRSI;
  s = 1.0;
  DIAGA1 = 'N';
  TRANSA1 = 'N';
  UPLO1 = 'U';
  SIDE1 = 'L';
  info_t = (ptrdiff_t)120;
  n_t = (ptrdiff_t)120;
  lda_t = (ptrdiff_t)120;
  ldb_t = (ptrdiff_t)120;
  dtrsm(&SIDE1, &UPLO1, &TRANSA1, &DIAGA1, &info_t, &n_t, &s, &b_x[0], &lda_t,
        &y[0], &ldb_t);
  st.site = &x_emlrtRSI;
  n1x = 0.0;
  j = 0;
  exitg1 = false;
  while ((!exitg1) && (j < 120)) {
    s = 0.0;
    for (b_i = 0; b_i < 120; b_i++) {
      s += muDoubleScalarAbs(x[b_i + 120 * j]);
    }

    if (muDoubleScalarIsNaN(s)) {
      n1x = rtNaN;
      exitg1 = true;
    } else {
      if (s > n1x) {
        n1x = s;
      }

      j++;
    }
  }

  n1xinv = 0.0;
  j = 0;
  exitg1 = false;
  while ((!exitg1) && (j < 120)) {
    s = 0.0;
    for (b_i = 0; b_i < 120; b_i++) {
      s += muDoubleScalarAbs(y[b_i + 120 * j]);
    }

    if (muDoubleScalarIsNaN(s)) {
      n1xinv = rtNaN;
      exitg1 = true;
    } else {
      if (s > n1xinv) {
        n1xinv = s;
      }

      j++;
    }
  }

  s = 1.0 / (n1x * n1xinv);
  if ((n1x == 0.0) || (n1xinv == 0.0) || (s == 0.0)) {
    b_st.site = &hb_emlrtRSI;
    warning(&b_st);
  } else {
    if (muDoubleScalarIsNaN(s) || (s < 2.2204460492503131E-16)) {
      b_st.site = &ib_emlrtRSI;
      b_y = NULL;
      m = emlrtCreateCharArray(2, &iv[0]);
      emlrtInitCharArrayR2013a(&b_st, 6, m, &rfmt[0]);
      emlrtAssign(&b_y, m);
      c_st.site = &sc_emlrtRSI;
      emlrt_marshallIn(&c_st, b_sprintf(&c_st, b_y, emlrt_marshallOut(s),
        &c_emlrtMCI), "<output of sprintf>", str);
      b_st.site = &ib_emlrtRSI;
      b_warning(&b_st, str);
    }
  }
}

/* End of code generation (inv.c) */
