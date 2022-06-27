/*
 * xzlarf.c
 *
 * Code generation for function 'xzlarf'
 *
 */

/* Include files */
#include "xzlarf.h"
#include "eml_int_forloop_overflow_check.h"
#include "projGradSolver_rt_data.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo rb_emlrtRSI = { 75, /* lineNo */
  "xzlarf",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzlarf.m"/* pathName */
};

static emlrtRSInfo sb_emlrtRSI = { 68, /* lineNo */
  "xzlarf",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzlarf.m"/* pathName */
};

static emlrtRSInfo tb_emlrtRSI = { 50, /* lineNo */
  "xzlarf",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzlarf.m"/* pathName */
};

static emlrtRSInfo ub_emlrtRSI = { 103,/* lineNo */
  "ilazlc",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzlarf.m"/* pathName */
};

static emlrtRSInfo vb_emlrtRSI = { 74, /* lineNo */
  "xgemv",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+refblas/xgemv.m"/* pathName */
};

/* Function Definitions */
void xzlarf(const emlrtStack *sp, int32_T m, int32_T n, int32_T iv0, real_T tau,
            real_T C[100], int32_T ic0, real_T work[10])
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  real_T c;
  int32_T b;
  int32_T colbottom;
  int32_T exitg1;
  int32_T i;
  int32_T ia;
  int32_T iac;
  int32_T ix;
  int32_T lastc;
  int32_T lastv;
  boolean_T exitg2;
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
  if (tau != 0.0) {
    lastv = m - 1;
    i = iv0 + m;
    while ((lastv + 1 > 0) && (C[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    st.site = &tb_emlrtRSI;
    lastc = n;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      i = ic0 + (lastc - 1) * 10;
      colbottom = i + lastv;
      b_st.site = &ub_emlrtRSI;
      if ((i <= colbottom) && (colbottom > 2147483646)) {
        c_st.site = &m_emlrtRSI;
        check_forloop_overflow_error(&c_st);
      }

      do {
        exitg1 = 0;
        if (i <= colbottom) {
          if (C[i - 1] != 0.0) {
            exitg1 = 1;
          } else {
            i++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = -1;
    lastc = 0;
  }

  if (lastv + 1 > 0) {
    st.site = &sb_emlrtRSI;
    b_st.site = &jb_emlrtRSI;
    if (lastc != 0) {
      c_st.site = &lb_emlrtRSI;
      if ((1 <= lastc) && (lastc > 2147483646)) {
        d_st.site = &m_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }

      if (0 <= lastc - 1) {
        memset(&work[0], 0, lastc * sizeof(real_T));
      }

      i = 0;
      colbottom = ic0 + 10 * (lastc - 1);
      for (iac = ic0; iac <= colbottom; iac += 10) {
        ix = iv0;
        c = 0.0;
        b = iac + lastv;
        c_st.site = &vb_emlrtRSI;
        if ((iac <= b) && (b > 2147483646)) {
          d_st.site = &m_emlrtRSI;
          check_forloop_overflow_error(&d_st);
        }

        for (ia = iac; ia <= b; ia++) {
          c += C[ia - 1] * C[ix - 1];
          ix++;
        }

        work[i] += c;
        i++;
      }
    }

    st.site = &rb_emlrtRSI;
    b_st.site = &mb_emlrtRSI;
    c_st.site = &nb_emlrtRSI;
    d_st.site = &ob_emlrtRSI;
    if (!(-tau == 0.0)) {
      i = ic0;
      colbottom = 0;
      e_st.site = &pb_emlrtRSI;
      for (iac = 0; iac < lastc; iac++) {
        if (work[colbottom] != 0.0) {
          c = work[colbottom] * -tau;
          ix = iv0;
          b = lastv + i;
          e_st.site = &qb_emlrtRSI;
          if ((i <= b) && (b > 2147483646)) {
            f_st.site = &m_emlrtRSI;
            check_forloop_overflow_error(&f_st);
          }

          for (ia = i; ia <= b; ia++) {
            C[ia - 1] += C[ix - 1] * c;
            ix++;
          }
        }

        colbottom++;
        i += 10;
      }
    }
  }
}

/* End of code generation (xzlarf.c) */
