/*
 * eig.c
 *
 * Code generation for function 'eig'
 *
 */

/* Include files */
#include "eig.h"
#include "anyNonFinite.h"
#include "eml_int_forloop_overflow_check.h"
#include "projGradSolver_rt_data.h"
#include "rt_nonfinite.h"
#include "warning.h"
#include "xdhseqr.h"
#include "xnrm2.h"
#include "xscal.h"
#include "xzlarf.h"
#include "lapacke.h"
#include "mwmathutil.h"
#include <stddef.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo g_emlrtRSI = { 93,  /* lineNo */
  "eig",                               /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/eig.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 139, /* lineNo */
  "eig",                               /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/eig.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 147, /* lineNo */
  "eig",                               /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/eig.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 21,  /* lineNo */
  "eigHermitianStandard",              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/private/eigHermitianStandard.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 35,  /* lineNo */
  "schur",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/schur.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 52,  /* lineNo */
  "schur",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/schur.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 54,  /* lineNo */
  "schur",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/schur.m"/* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 83,  /* lineNo */
  "schur",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/schur.m"/* pathName */
};

static emlrtRSInfo s_emlrtRSI = { 18,  /* lineNo */
  "xgehrd",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgehrd.m"/* pathName */
};

static emlrtRSInfo t_emlrtRSI = { 31,  /* lineNo */
  "xzgehrd",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzgehrd.m"/* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 35,  /* lineNo */
  "xzgehrd",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzgehrd.m"/* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 43,  /* lineNo */
  "xzgehrd",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzgehrd.m"/* pathName */
};

static emlrtRSInfo w_emlrtRSI = { 20,  /* lineNo */
  "xzlarfg",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzlarfg.m"/* pathName */
};

static emlrtRSInfo x_emlrtRSI = { 41,  /* lineNo */
  "xzlarfg",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzlarfg.m"/* pathName */
};

static emlrtRSInfo y_emlrtRSI = { 53,  /* lineNo */
  "xzlarfg",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzlarfg.m"/* pathName */
};

static emlrtRSInfo ab_emlrtRSI = { 68, /* lineNo */
  "xzlarfg",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzlarfg.m"/* pathName */
};

static emlrtRSInfo cb_emlrtRSI = { 81, /* lineNo */
  "xzlarfg",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzlarfg.m"/* pathName */
};

static emlrtRSInfo hb_emlrtRSI = { 84, /* lineNo */
  "xzlarf",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzlarf.m"/* pathName */
};

static emlrtRSInfo ib_emlrtRSI = { 91, /* lineNo */
  "xzlarf",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xzlarf.m"/* pathName */
};

static emlrtRSInfo kb_emlrtRSI = { 58, /* lineNo */
  "xgemv",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+refblas/xgemv.m"/* pathName */
};

static emlrtRSInfo wb_emlrtRSI = { 32, /* lineNo */
  "xhseqr",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/xhseqr.m"/* pathName */
};

static emlrtRSInfo xb_emlrtRSI = { 13, /* lineNo */
  "xdhseqr",                           /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xdhseqr.m"/* pathName */
};

static emlrtRSInfo qc_emlrtRSI = { 59, /* lineNo */
  "eigStandard",                       /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/private/eigStandard.m"/* pathName */
};

static emlrtRSInfo rc_emlrtRSI = { 44, /* lineNo */
  "eigStandard",                       /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/private/eigStandard.m"/* pathName */
};

static emlrtRSInfo sc_emlrtRSI = { 38, /* lineNo */
  "xgeev",                             /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgeev.m"/* pathName */
};

static emlrtRSInfo tc_emlrtRSI = { 143,/* lineNo */
  "ceval_xgeev",                       /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/xgeev.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 47,  /* lineNo */
  13,                                  /* colNo */
  "infocheck",                         /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/infocheck.m"/* pName */
};

static emlrtRTEInfo b_emlrtRTEI = { 44,/* lineNo */
  13,                                  /* colNo */
  "infocheck",                         /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+lapack/infocheck.m"/* pName */
};

/* Function Definitions */
void eig(const emlrtStack *sp, const real_T A[100], creal_T V[10])
{
  static const char_T fname[14] = { 'L', 'A', 'P', 'A', 'C', 'K', 'E', '_', 'd',
    'g', 'e', 'e', 'v', 'x' };

  ptrdiff_t ihi_t;
  ptrdiff_t ilo_t;
  ptrdiff_t info_t;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  emlrtStack i_st;
  emlrtStack j_st;
  emlrtStack k_st;
  emlrtStack st;
  real_T T[100];
  real_T scale[10];
  real_T wimag[10];
  real_T wreal[10];
  real_T tau[9];
  real_T abnrm;
  real_T rconde;
  real_T rcondv;
  real_T vleft;
  real_T vright;
  int32_T alpha1_tmp;
  int32_T b;
  int32_T b_i;
  int32_T exitg1;
  int32_T i;
  int32_T ia;
  int32_T ic0;
  int32_T im1n_tmp;
  int32_T in;
  int32_T ix;
  int32_T j;
  int32_T jy;
  int32_T knt;
  int32_T lastc;
  int32_T lastv;
  boolean_T exitg2;
  boolean_T p;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &g_emlrtRSI;
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
  j_st.prev = &i_st;
  j_st.tls = i_st.tls;
  k_st.prev = &j_st;
  k_st.tls = j_st.tls;
  if (anyNonFinite(A)) {
    for (i = 0; i < 10; i++) {
      V[i].re = rtNaN;
      V[i].im = 0.0;
    }
  } else {
    p = true;
    j = 0;
    exitg2 = false;
    while ((!exitg2) && (j < 10)) {
      i = 0;
      do {
        exitg1 = 0;
        if (i <= j) {
          if (!(A[i + 10 * j] == A[j + 10 * i])) {
            p = false;
            exitg1 = 1;
          } else {
            i++;
          }
        } else {
          j++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }

    if (p) {
      st.site = &h_emlrtRSI;
      b_st.site = &n_emlrtRSI;
      memcpy(&T[0], &A[0], 100U * sizeof(real_T));
      c_st.site = &o_emlrtRSI;
      if (anyNonFinite(A)) {
        for (b_i = 0; b_i < 100; b_i++) {
          T[b_i] = rtNaN;
        }

        b_i = 2;
        for (j = 0; j < 9; j++) {
          if (b_i <= 10) {
            memset(&T[(j * 10 + b_i) + -1], 0, (11 - b_i) * sizeof(real_T));
          }

          b_i++;
        }
      } else {
        c_st.site = &p_emlrtRSI;
        d_st.site = &s_emlrtRSI;
        memset(&scale[0], 0, 10U * sizeof(real_T));
        for (i = 0; i < 9; i++) {
          im1n_tmp = i * 10;
          in = (i + 1) * 10;
          alpha1_tmp = (i + 10 * i) + 1;
          vright = T[alpha1_tmp];
          b_i = i + 3;
          b_i = muIntScalarMin_sint32(b_i, 10) + im1n_tmp;
          e_st.site = &t_emlrtRSI;
          tau[i] = 0.0;
          f_st.site = &w_emlrtRSI;
          vleft = xnrm2(&f_st, 8 - i, T, b_i);
          if (vleft != 0.0) {
            vleft = muDoubleScalarHypot(vright, vleft);
            if (vright >= 0.0) {
              vleft = -vleft;
            }

            if (muDoubleScalarAbs(vleft) < 1.0020841800044864E-292) {
              knt = 0;
              do {
                knt++;
                f_st.site = &x_emlrtRSI;
                xscal(&f_st, 8 - i, 9.9792015476736E+291, T, b_i);
                vleft *= 9.9792015476736E+291;
                vright *= 9.9792015476736E+291;
              } while (!(muDoubleScalarAbs(vleft) >= 1.0020841800044864E-292));

              f_st.site = &y_emlrtRSI;
              vleft = xnrm2(&f_st, 8 - i, T, b_i);
              vleft = muDoubleScalarHypot(vright, vleft);
              if (vright >= 0.0) {
                vleft = -vleft;
              }

              tau[i] = (vleft - vright) / vleft;
              f_st.site = &ab_emlrtRSI;
              xscal(&f_st, 8 - i, 1.0 / (vright - vleft), T, b_i);
              f_st.site = &bb_emlrtRSI;
              if ((1 <= knt) && (knt > 2147483646)) {
                g_st.site = &m_emlrtRSI;
                check_forloop_overflow_error(&g_st);
              }

              for (b_i = 0; b_i < knt; b_i++) {
                vleft *= 1.0020841800044864E-292;
              }

              vright = vleft;
            } else {
              tau[i] = (vleft - vright) / vleft;
              f_st.site = &cb_emlrtRSI;
              xscal(&f_st, 8 - i, 1.0 / (T[(i + 10 * i) + 1] - vleft), T, b_i);
              vright = vleft;
            }
          }

          T[alpha1_tmp] = 1.0;
          jy = (i + im1n_tmp) + 1;
          ic0 = in + 1;
          e_st.site = &u_emlrtRSI;
          if (tau[i] != 0.0) {
            lastv = 8 - i;
            b_i = (jy - i) + 8;
            while ((lastv + 1 > 0) && (T[b_i] == 0.0)) {
              lastv--;
              b_i--;
            }

            lastc = 10;
            exitg2 = false;
            while ((!exitg2) && (lastc > 0)) {
              knt = in + lastc;
              ia = knt;
              do {
                exitg1 = 0;
                if (ia <= knt + lastv * 10) {
                  if (T[ia - 1] != 0.0) {
                    exitg1 = 1;
                  } else {
                    ia += 10;
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
            f_st.site = &hb_emlrtRSI;
            g_st.site = &jb_emlrtRSI;
            if (lastc != 0) {
              h_st.site = &lb_emlrtRSI;
              if (0 <= lastc - 1) {
                memset(&scale[0], 0, lastc * sizeof(real_T));
              }

              ix = jy;
              b_i = (in + 10 * lastv) + 1;
              for (j = ic0; j <= b_i; j += 10) {
                knt = 0;
                b = (j + lastc) - 1;
                h_st.site = &kb_emlrtRSI;
                for (ia = j; ia <= b; ia++) {
                  scale[knt] += T[ia - 1] * T[ix];
                  knt++;
                }

                ix++;
              }
            }

            f_st.site = &ib_emlrtRSI;
            g_st.site = &mb_emlrtRSI;
            h_st.site = &nb_emlrtRSI;
            i_st.site = &ob_emlrtRSI;
            if (!(-tau[i] == 0.0)) {
              knt = in;
              j_st.site = &pb_emlrtRSI;
              for (j = 0; j <= lastv; j++) {
                if (T[jy] != 0.0) {
                  vleft = T[jy] * -tau[i];
                  ix = 0;
                  b_i = knt + 1;
                  b = lastc + knt;
                  j_st.site = &qb_emlrtRSI;
                  if ((knt + 1 <= b) && (b > 2147483646)) {
                    k_st.site = &m_emlrtRSI;
                    check_forloop_overflow_error(&k_st);
                  }

                  for (ic0 = b_i; ic0 <= b; ic0++) {
                    T[ic0 - 1] += scale[ix] * vleft;
                    ix++;
                  }
                }

                jy++;
                knt += 10;
              }
            }
          }

          e_st.site = &v_emlrtRSI;
          xzlarf(&e_st, 9 - i, 9 - i, (i + im1n_tmp) + 2, tau[i], T, (i + in) +
                 2, scale);
          T[alpha1_tmp] = vright;
        }

        c_st.site = &q_emlrtRSI;
        d_st.site = &wb_emlrtRSI;
        e_st.site = &xb_emlrtRSI;
        knt = eml_dlahqr(&e_st, T);
        b_i = 4;
        for (j = 0; j < 7; j++) {
          if (b_i <= 10) {
            memset(&T[(j * 10 + b_i) + -1], 0, (11 - b_i) * sizeof(real_T));
          }

          b_i++;
        }

        if (knt != 0) {
          c_st.site = &r_emlrtRSI;
          warning(&c_st);
        }
      }

      for (b_i = 0; b_i < 10; b_i++) {
        scale[b_i] = T[b_i + 10 * b_i];
      }

      for (i = 0; i < 10; i++) {
        V[i].re = scale[i];
        V[i].im = 0.0;
      }
    } else {
      st.site = &i_emlrtRSI;
      b_st.site = &rc_emlrtRSI;
      c_st.site = &sc_emlrtRSI;
      memcpy(&T[0], &A[0], 100U * sizeof(real_T));
      info_t = LAPACKE_dgeevx(102, 'B', 'N', 'N', 'N', (ptrdiff_t)10, &T[0],
        (ptrdiff_t)10, &wreal[0], &wimag[0], &vleft, (ptrdiff_t)1, &vright,
        (ptrdiff_t)1, &ilo_t, &ihi_t, &scale[0], &abnrm, &rconde, &rcondv);
      knt = (int32_T)info_t;
      d_st.site = &tc_emlrtRSI;
      if (knt < 0) {
        if (knt == -1010) {
          emlrtErrorWithMessageIdR2018a(&d_st, &b_emlrtRTEI, "MATLAB:nomem",
            "MATLAB:nomem", 0);
        } else {
          emlrtErrorWithMessageIdR2018a(&d_st, &emlrtRTEI,
            "Coder:toolbox:LAPACKCallErrorInfo",
            "Coder:toolbox:LAPACKCallErrorInfo", 5, 4, 14, fname, 12, knt);
        }
      }

      for (i = 0; i < 10; i++) {
        V[i].re = wreal[i];
        V[i].im = wimag[i];
      }

      if (knt != 0) {
        b_st.site = &qc_emlrtRSI;
        b_warning(&b_st);
      }
    }
  }
}

/* End of code generation (eig.c) */
