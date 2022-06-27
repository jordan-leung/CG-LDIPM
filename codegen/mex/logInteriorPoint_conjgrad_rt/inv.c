/*
 * inv.c
 *
 * Code generation for function 'inv'
 *
 */

/* Include files */
#include "inv.h"
#include "eml_int_forloop_overflow_check.h"
#include "logInteriorPoint_conjgrad_rt_mexutil.h"
#include "rt_nonfinite.h"
#include "warning.h"
#include "blas.h"
#include "mwmathutil.h"
#include <stddef.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo s_emlrtRSI = { 21,  /* lineNo */
  "inv",                               /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo t_emlrtRSI = { 22,  /* lineNo */
  "inv",                               /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 173, /* lineNo */
  "invNxN",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 190, /* lineNo */
  "invNxN",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo w_emlrtRSI = { 30,  /* lineNo */
  "xgetrf",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgetrf.m"/* pathName */
};

static emlrtRSInfo x_emlrtRSI = { 50,  /* lineNo */
  "xzgetrf",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzgetrf.m"/* pathName */
};

static emlrtRSInfo y_emlrtRSI = { 58,  /* lineNo */
  "xzgetrf",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzgetrf.m"/* pathName */
};

static emlrtRSInfo ab_emlrtRSI = { 21, /* lineNo */
  "eml_int_forloop_overflow_check",    /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"/* pathName */
};

static emlrtRSInfo bb_emlrtRSI = { 45, /* lineNo */
  "xgeru",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+blas/xgeru.m"/* pathName */
};

static emlrtRSInfo cb_emlrtRSI = { 45, /* lineNo */
  "xger",                              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+blas/xger.m"/* pathName */
};

static emlrtRSInfo db_emlrtRSI = { 15, /* lineNo */
  "xger",                              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+refblas/xger.m"/* pathName */
};

static emlrtRSInfo eb_emlrtRSI = { 54, /* lineNo */
  "xgerx",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+refblas/xgerx.m"/* pathName */
};

static emlrtRSInfo fb_emlrtRSI = { 67, /* lineNo */
  "xtrsm",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+blas/xtrsm.m"/* pathName */
};

static emlrtRSInfo gb_emlrtRSI = { 81, /* lineNo */
  "xtrsm_blas",                        /* fcnName */
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

static emlrtRSInfo jc_emlrtRSI = { 53, /* lineNo */
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
static void m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[14]);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[14])
{
  m_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
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

static void m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[14])
{
  static const int32_T dims[2] = { 1, 14 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", false, 2U, dims);
  emlrtImportCharArrayR2015b(sp, src, &ret[0], 14);
  emlrtDestroyArray(&src);
}

void inv(const emlrtStack *sp, const real_T x[100], real_T y[100])
{
  static const int32_T iv[2] = { 1, 6 };

  static const char_T rfmt[6] = { '%', '1', '4', '.', '6', 'e' };

  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  emlrtStack i_st;
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *m;
  real_T b_x[100];
  real_T n1x;
  real_T n1xinv;
  real_T s;
  real_T smax;
  int32_T b;
  int32_T b_tmp;
  int32_T ix;
  int32_T iy;
  int32_T j;
  int32_T jA;
  int32_T jp1j;
  int32_T k;
  int32_T n;
  char_T str[14];
  char_T DIAGA1;
  char_T SIDE1;
  char_T TRANSA1;
  char_T UPLO1;
  int8_T ipiv[10];
  int8_T p[10];
  int8_T i;
  boolean_T exitg1;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &s_emlrtRSI;
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
  h_st.prev = &g_st;
  h_st.tls = g_st.tls;
  i_st.prev = &h_st;
  i_st.tls = h_st.tls;
  memset(&y[0], 0, 100U * sizeof(real_T));
  b_st.site = &u_emlrtRSI;
  memcpy(&b_x[0], &x[0], 100U * sizeof(real_T));
  c_st.site = &w_emlrtRSI;
  for (b = 0; b < 10; b++) {
    ipiv[b] = (int8_T)(b + 1);
  }

  for (j = 0; j < 9; j++) {
    b_tmp = j * 11;
    jp1j = b_tmp + 2;
    n = 10 - j;
    iy = 0;
    ix = b_tmp;
    smax = muDoubleScalarAbs(b_x[b_tmp]);
    for (k = 2; k <= n; k++) {
      ix++;
      s = muDoubleScalarAbs(b_x[ix]);
      if (s > smax) {
        iy = k - 1;
        smax = s;
      }
    }

    if (b_x[b_tmp + iy] != 0.0) {
      if (iy != 0) {
        iy += j;
        ipiv[j] = (int8_T)(iy + 1);
        ix = j;
        for (k = 0; k < 10; k++) {
          smax = b_x[ix];
          b_x[ix] = b_x[iy];
          b_x[iy] = smax;
          ix += 10;
          iy += 10;
        }
      }

      b = (b_tmp - j) + 10;
      d_st.site = &x_emlrtRSI;
      for (n = jp1j; n <= b; n++) {
        b_x[n - 1] /= b_x[b_tmp];
      }
    }

    n = 8 - j;
    iy = b_tmp + 10;
    d_st.site = &y_emlrtRSI;
    e_st.site = &bb_emlrtRSI;
    f_st.site = &cb_emlrtRSI;
    g_st.site = &db_emlrtRSI;
    jA = b_tmp + 12;
    for (jp1j = 0; jp1j <= n; jp1j++) {
      smax = b_x[iy];
      if (b_x[iy] != 0.0) {
        ix = b_tmp + 1;
        b = (jA - j) + 8;
        h_st.site = &eb_emlrtRSI;
        if ((jA <= b) && (b > 2147483646)) {
          i_st.site = &ab_emlrtRSI;
          check_forloop_overflow_error(&i_st);
        }

        for (k = jA; k <= b; k++) {
          b_x[k - 1] += b_x[ix] * -smax;
          ix++;
        }
      }

      iy += 10;
      jA += 10;
    }
  }

  for (b = 0; b < 10; b++) {
    p[b] = (int8_T)(b + 1);
  }

  for (k = 0; k < 9; k++) {
    i = ipiv[k];
    if (i > k + 1) {
      iy = p[i - 1];
      p[i - 1] = p[k];
      p[k] = (int8_T)iy;
    }
  }

  for (k = 0; k < 10; k++) {
    iy = 10 * (p[k] - 1);
    y[k + iy] = 1.0;
    for (j = k + 1; j < 11; j++) {
      b = (j + iy) - 1;
      if (y[b] != 0.0) {
        jA = j + 1;
        for (n = jA; n < 11; n++) {
          jp1j = (n + iy) - 1;
          y[jp1j] -= y[b] * b_x[(n + 10 * (j - 1)) - 1];
        }
      }
    }
  }

  b_st.site = &v_emlrtRSI;
  c_st.site = &fb_emlrtRSI;
  d_st.site = &gb_emlrtRSI;
  smax = 1.0;
  DIAGA1 = 'N';
  TRANSA1 = 'N';
  UPLO1 = 'U';
  SIDE1 = 'L';
  m_t = (ptrdiff_t)10;
  n_t = (ptrdiff_t)10;
  lda_t = (ptrdiff_t)10;
  ldb_t = (ptrdiff_t)10;
  dtrsm(&SIDE1, &UPLO1, &TRANSA1, &DIAGA1, &m_t, &n_t, &smax, &b_x[0], &lda_t,
        &y[0], &ldb_t);
  st.site = &t_emlrtRSI;
  n1x = 0.0;
  j = 0;
  exitg1 = false;
  while ((!exitg1) && (j < 10)) {
    s = 0.0;
    for (n = 0; n < 10; n++) {
      s += muDoubleScalarAbs(x[n + 10 * j]);
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
  while ((!exitg1) && (j < 10)) {
    s = 0.0;
    for (n = 0; n < 10; n++) {
      s += muDoubleScalarAbs(y[n + 10 * j]);
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

  smax = 1.0 / (n1x * n1xinv);
  if ((n1x == 0.0) || (n1xinv == 0.0) || (smax == 0.0)) {
    b_st.site = &hb_emlrtRSI;
    warning(&b_st);
  } else {
    if (muDoubleScalarIsNaN(smax) || (smax < 2.2204460492503131E-16)) {
      b_st.site = &ib_emlrtRSI;
      b_y = NULL;
      m = emlrtCreateCharArray(2, &iv[0]);
      emlrtInitCharArrayR2013a(&b_st, 6, m, &rfmt[0]);
      emlrtAssign(&b_y, m);
      c_st.site = &jc_emlrtRSI;
      emlrt_marshallIn(&c_st, b_sprintf(&c_st, b_y, emlrt_marshallOut(smax),
        &c_emlrtMCI), "<output of sprintf>", str);
      b_st.site = &ib_emlrtRSI;
      b_warning(&b_st, str);
    }
  }
}

/* End of code generation (inv.c) */
