/*
 * logInteriorPoint_rt.c
 *
 * Code generation for function 'logInteriorPoint_rt'
 *
 */

/* Include files */
#include "logInteriorPoint_rt.h"
#include "chol.h"
#include "inv.h"
#include "logInteriorPoint_rt_data.h"
#include "logInteriorPoint_rt_types.h"
#include "mldivide.h"
#include "mtimes.h"
#include "rt_nonfinite.h"
#include "tic.h"
#include "toc.h"
#include "mwmathutil.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 13,    /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 24,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 33,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 36,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 43,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 44,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 48,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 49,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 58,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 68,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 88,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 91,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 107, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 109, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 113, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 114, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 118, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 119, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo s_emlrtRSI = { 128, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo t_emlrtRSI = { 138, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 142, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 144, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo rb_emlrtRSI = { 34, /* lineNo */
  "chol",                              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/chol.m"/* pathName */
};

static emlrtRSInfo xb_emlrtRSI = { 224,/* lineNo */
  "solveNewtonStep_decomp",            /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo kc_emlrtRSI = { 188,/* lineNo */
  "muStarSolve",                       /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo nc_emlrtRSI = { 204,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 13,  /* lineNo */
  9,                                   /* colNo */
  "sqrt",                              /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/elfun/sqrt.m"/* pName */
};

/* Function Declarations */
static void muStarSolve(const emlrtStack *sp, const real_T d0[120], const real_T
  d1[120], real_T dinfmax, real_T mu_f, real_T *muStar, real_T d[120]);
static void solveNewtonStep(logInteriorPoint_rtStackData *SD, const emlrtStack
  *sp, real_T mu, const real_T v[120], const real_T const_W[14400], const real_T
  const_c[120], const real_T const_A[14400], const real_T const_b[120], real_T
  x[120]);
static void solveNewtonStep_decomp(logInteriorPoint_rtStackData *SD, const
  emlrtStack *sp, real_T mu, const real_T v[120], const real_T const_c[120],
  const real_T const_A[14400], const real_T const_b[120], const real_T dG[14400],
  real_T x[120], real_T d[120]);

/* Function Definitions */
static void muStarSolve(const emlrtStack *sp, const real_T d0[120], const real_T
  d1[120], real_T dinfmax, real_T mu_f, real_T *muStar, real_T d[120])
{
  emlrtStack st;
  real_T lower_bound;
  real_T lower_bound_i;
  real_T temp;
  real_T upper_bound;
  real_T upper_bound_i;
  int32_T i;
  boolean_T exitg1;
  st.prev = sp;
  st.tls = sp->tls;

  /*  // Returns the largest k satisfying | d0 + k * d1 |_{infty} \le dinfmax */
  upper_bound = 1.0E+14;
  lower_bound = 0.0;
  i = 0;
  exitg1 = false;
  while ((!exitg1) && (i < 120)) {
    upper_bound_i = (dinfmax - d0[i]) / d1[i];
    lower_bound_i = (-dinfmax - d0[i]) / d1[i];

    /*  Neither is either lower or upper necessarily, so just reorder */
    /*  according to whichever is larger */
    if (lower_bound_i > upper_bound_i) {
      temp = upper_bound_i;
      upper_bound_i = lower_bound_i;
      lower_bound_i = temp;
    }

    /*  Check whether or not this is larger/small than the previous upper and */
    /*  lower bounds */
    if (upper_bound_i < upper_bound) {
      upper_bound = upper_bound_i;
    }

    if (lower_bound_i > lower_bound) {
      lower_bound = lower_bound_i;
    }

    /*  Finally, check if lower_bound > upper_bound, indicating that  the */
    /*  solution is infeasible */
    if (lower_bound > upper_bound) {
      upper_bound_i = 0.0;
      exitg1 = true;
    } else {
      upper_bound_i = upper_bound;

      /*  note k = 1/sqrt(mu), so we want the largest k */
      i++;
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }
  }

  upper_bound_i = 1.0 / upper_bound_i;
  *muStar = upper_bound_i * upper_bound_i;

  /*  returns muStar = Inf if infeasible */
  if (*muStar < mu_f) {
    /*  make sure muStar isn't too small or too large for numerical reasons */
    *muStar = mu_f;
  } else {
    if (*muStar > 1.0E+10) {
      *muStar = rtInf;
    }
  }

  st.site = &kc_emlrtRSI;
  if (*muStar < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  upper_bound_i = 1.0 / muDoubleScalarSqrt(*muStar);
  for (i = 0; i < 120; i++) {
    d[i] = d0[i] + upper_bound_i * d1[i];
  }
}

static void solveNewtonStep(logInteriorPoint_rtStackData *SD, const emlrtStack
  *sp, real_T mu, const real_T v[120], const real_T const_W[14400], const real_T
  const_c[120], const real_T const_A[14400], const real_T const_b[120], real_T
  x[120])
{
  emlrtStack st;
  real_T Qvec[120];
  real_T y;
  int32_T i;
  int32_T k;
  st.prev = sp;
  st.tls = sp->tls;

  /*  Construct Q matrix */
  for (k = 0; k < 120; k++) {
    Qvec[k] = muDoubleScalarExp(2.0 * v[k]);
  }

  memset(&SD->u1.f1.Q[0], 0, 14400U * sizeof(real_T));
  for (k = 0; k < 120; k++) {
    SD->u1.f1.Q[k + 120 * k] = Qvec[k];
  }

  /*  Solve for x */
  mtimes(const_A, SD->u1.f1.Q, SD->u1.f1.y);
  b_mtimes(SD->u1.f1.y, const_A, SD->u1.f1.b_y);
  st.site = &nc_emlrtRSI;
  if (mu < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  y = 2.0 * muDoubleScalarSqrt(mu);
  for (k = 0; k < 120; k++) {
    Qvec[k] = muDoubleScalarExp(v[k]);
    for (i = 0; i < 120; i++) {
      SD->u1.f1.y[i + 120 * k] = y * const_A[k + 120 * i];
    }
  }

  c_mtimes(SD->u1.f1.y, Qvec, x);
  mtimes(const_A, SD->u1.f1.Q, SD->u1.f1.y);
  c_mtimes(SD->u1.f1.y, const_b, Qvec);
  for (i = 0; i < 120; i++) {
    x[i] -= const_c[i] + Qvec[i];
  }

  for (i = 0; i < 14400; i++) {
    SD->u1.f1.b_y[i] += const_W[i];
  }

  st.site = &nc_emlrtRSI;
  mldivide(&st, SD->u1.f1.b_y, x);

  /*  Solve for d */
  c_mtimes(const_A, x, Qvec);
}

static void solveNewtonStep_decomp(logInteriorPoint_rtStackData *SD, const
  emlrtStack *sp, real_T mu, const real_T v[120], const real_T const_c[120],
  const real_T const_A[14400], const real_T const_b[120], const real_T dG[14400],
  real_T x[120], real_T d[120])
{
  emlrtStack st;
  real_T Qvec[120];
  real_T y;
  real_T y_tmp;
  int32_T i;
  int32_T k;
  st.prev = sp;
  st.tls = sp->tls;

  /*  Construct Q matrix */
  for (k = 0; k < 120; k++) {
    Qvec[k] = muDoubleScalarExp(2.0 * v[k]);
  }

  memset(&SD->u1.f0.Q[0], 0, 14400U * sizeof(real_T));
  for (k = 0; k < 120; k++) {
    SD->u1.f0.Q[k + 120 * k] = Qvec[k];
  }

  /*  Solve for x */
  /*  x =  dG\(2*sqrt(mu)*A'*exp(v) - (c + A'*Q*b)); */
  st.site = &xb_emlrtRSI;
  if (mu < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  y_tmp = muDoubleScalarSqrt(mu);
  y = 2.0 * y_tmp;
  for (k = 0; k < 120; k++) {
    Qvec[k] = muDoubleScalarExp(v[k]);
    for (i = 0; i < 120; i++) {
      SD->u1.f0.y[i + 120 * k] = y * const_A[k + 120 * i];
    }
  }

  c_mtimes(SD->u1.f0.y, Qvec, x);
  mtimes(const_A, SD->u1.f0.Q, SD->u1.f0.y);
  c_mtimes(SD->u1.f0.y, const_b, Qvec);
  for (i = 0; i < 120; i++) {
    x[i] -= const_c[i] + Qvec[i];
    for (k = 0; k < 120; k++) {
      SD->u1.f0.Q[k + 120 * i] = dG[i + 120 * k];
    }
  }

  st.site = &xb_emlrtRSI;
  mldivide(&st, SD->u1.f0.Q, x);
  st.site = &xb_emlrtRSI;
  mldivide(&st, dG, x);
  y_tmp = 1.0 / y_tmp;
  for (k = 0; k < 120; k++) {
    d[k] = muDoubleScalarExp(v[k]);
  }

  c_mtimes(const_A, x, Qvec);
  for (i = 0; i < 120; i++) {
    d[i] = 1.0 - y_tmp * d[i] * (Qvec[i] + const_b[i]);
  }
}

void logInteriorPoint_rt(logInteriorPoint_rtStackData *SD, const emlrtStack *sp,
  const real_T W[14400], const real_T c[120], real_T A[14400], const real_T
  bineq[120], real_T mu_f, real_T mu_0, const real_T v0[120], real_T maxIter,
  real_T x[120], real_T *mu, real_T *execTime, real_T *numIter)
{
  emlrtStack b_st;
  emlrtStack st;
  real_T Qvec[120];
  real_T b_c0[120];
  real_T d[120];
  real_T d1Hat[120];
  real_T d2Hat[120];
  real_T v[120];
  real_T a;
  real_T c0;
  real_T k2;
  real_T mu1;
  real_T mu2;
  int32_T exitg1;
  int32_T k;
  int32_T startFlag;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;

  /*  min 0.5*x'*W*x + c'*x   subject to:  A*x <= b */
  /*  Get size variables */
  /*  First, change variables to Ax + b >= 0... This is just for uniformity */
  /*  with quadprog's inputs. */
  for (startFlag = 0; startFlag < 14400; startFlag++) {
    A[startFlag] = -A[startFlag];
  }

  /*  Pack */
  st.site = &emlrtRSI;
  inv(&st, W, SD->f2.Q);

  /*  Max amount of bisection iterations */
  *numIter = 0.0;

  /*  number of newton iterations performed */
  /*  --------------------- FIRST NEWTON ITERATION --------------------- */
  /*  If we provide a warm-start */
  st.site = &b_emlrtRSI;
  tic(&st);

  /*  On initialization pass we first check if we can find a muStar value */
  /*  within N iterations of biscection... If not, then we just initialize with */
  /*  mu0 */
  /*  Construct Q matrix */
  for (k = 0; k < 120; k++) {
    Qvec[k] = muDoubleScalarExp(2.0 * v0[k]);
  }

  memset(&SD->f2.Q[0], 0, 14400U * sizeof(real_T));
  for (startFlag = 0; startFlag < 120; startFlag++) {
    SD->f2.Q[startFlag + 120 * startFlag] = Qvec[startFlag];
  }

  st.site = &c_emlrtRSI;
  mtimes(A, SD->f2.Q, SD->f2.y);
  st.site = &c_emlrtRSI;
  b_mtimes(SD->f2.y, A, SD->f2.Q);

  /*  dG = decomposition(G,'lu'); */
  /*  dG = decomposition(1/2*(G+G'),'chol'); */
  st.site = &d_emlrtRSI;
  for (startFlag = 0; startFlag < 14400; startFlag++) {
    SD->f2.Q[startFlag] += W[startFlag];
  }

  b_st.site = &rb_emlrtRSI;
  cholesky(&b_st, SD->f2.Q);

  /*  First, we sample two points of (mu,d) and solve the linear system */
  /*  generated by d = d0 + k*d1, where k = 1/sqrt(mu) */
  /*  Sample point 1 */
  st.site = &e_emlrtRSI;
  if (mu_0 < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  st.site = &f_emlrtRSI;
  solveNewtonStep_decomp(SD, &st, mu_0, v0, c, A, bineq, SD->f2.Q, Qvec, d1Hat);

  /*  Sample point 2 */
  mu2 = 0.1 * mu_0;
  st.site = &g_emlrtRSI;
  k2 = 1.0 / muDoubleScalarSqrt(mu2);
  st.site = &h_emlrtRSI;
  solveNewtonStep_decomp(SD, &st, mu2, v0, c, A, bineq, SD->f2.Q, Qvec, d2Hat);

  /*  Obtain affine representation of d = d0 + k*d1 */
  mu2 = 1.0 / muDoubleScalarSqrt(mu_0) - k2;
  c0 = -k2 / mu2;
  a = 1.0 / mu2;

  /*  Solve for muStar */
  for (startFlag = 0; startFlag < 120; startFlag++) {
    mu2 = d1Hat[startFlag];
    k2 = d2Hat[startFlag];
    b_c0[startFlag] = c0 * mu2 + (1.0 - c0) * k2;
    mu2 = a * (mu2 - k2);
    d1Hat[startFlag] = mu2;
  }

  st.site = &i_emlrtRSI;
  muStarSolve(&st, b_c0, d1Hat, 1.0, mu_f, mu, d);

  /*  If we found a muStar (or a feasible point), then use these values */
  if (!muDoubleScalarIsInf(*mu)) {
    startFlag = 1;

    /*  Update v - this updates regardless of whether or not we find a */
    /*  feasible mu... Either just iterates off mu0, or to muStar */
    st.site = &j_emlrtRSI;
    a = 0.0;
    for (k = 0; k < 120; k++) {
      c0 = muDoubleScalarAbs(d[k]);
      if (muDoubleScalarIsNaN(c0) || (c0 > a)) {
        a = c0;
      }
    }

    a = muDoubleScalarMin(1.0, 1.0 / (a * a));
    for (k = 0; k < 120; k++) {
      v[k] = v0[k] + a * d[k];
    }

    *numIter = 1.0;
  } else {
    /*  Otherwise, truly give  up and cold start */
    *mu = mu_0;

    /*  under the update at the end */
    for (k = 0; k < 120; k++) {
      d[k] *= rtInf;
      v[k] = 0.0;
    }

    startFlag = 0;
  }

  /*  If we have a pair (mu,d) with d > 1 and we iterate the Newton */
  /*  algorithm until we obtain convergence. If we found a muStar from our */
  /*  initialization, this loop will be skipped. */
  if (startFlag < 1) {
    do {
      exitg1 = 0;
      k2 = 0.0;
      for (k = 0; k < 120; k++) {
        c0 = muDoubleScalarAbs(d[k]);
        if (muDoubleScalarIsNaN(c0) || (c0 > k2)) {
          k2 = c0;
        }
      }

      if ((k2 > 1.0) && (!(*numIter >= maxIter))) {
        /*  Solve for d */
        st.site = &k_emlrtRSI;

        /*  Construct Q matrix */
        for (k = 0; k < 120; k++) {
          Qvec[k] = muDoubleScalarExp(2.0 * v[k]);
        }

        memset(&SD->f2.Q[0], 0, 14400U * sizeof(real_T));
        for (startFlag = 0; startFlag < 120; startFlag++) {
          SD->f2.Q[startFlag + 120 * startFlag] = Qvec[startFlag];
        }

        /*  Solve for x */
        mtimes(A, SD->f2.Q, SD->f2.y);
        b_mtimes(SD->f2.y, A, SD->f2.b_y);
        b_st.site = &nc_emlrtRSI;
        if (*mu < 0.0) {
          emlrtErrorWithMessageIdR2018a(&b_st, &emlrtRTEI,
            "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError",
            3, 4, 4, "sqrt");
        }

        mu2 = muDoubleScalarSqrt(*mu);
        k2 = 2.0 * mu2;
        for (k = 0; k < 120; k++) {
          Qvec[k] = muDoubleScalarExp(v[k]);
          for (startFlag = 0; startFlag < 120; startFlag++) {
            SD->f2.y[startFlag + 120 * k] = k2 * A[k + 120 * startFlag];
          }
        }

        c_mtimes(SD->f2.y, Qvec, d2Hat);
        mtimes(A, SD->f2.Q, SD->f2.y);
        c_mtimes(SD->f2.y, bineq, Qvec);
        for (startFlag = 0; startFlag < 120; startFlag++) {
          d2Hat[startFlag] -= c[startFlag] + Qvec[startFlag];
        }

        for (startFlag = 0; startFlag < 14400; startFlag++) {
          SD->f2.b_y[startFlag] += W[startFlag];
        }

        b_st.site = &nc_emlrtRSI;
        mldivide(&b_st, SD->f2.b_y, d2Hat);

        /*  Solve for d */
        a = 1.0 / mu2;
        for (k = 0; k < 120; k++) {
          d[k] = muDoubleScalarExp(v[k]);
        }

        c_mtimes(A, d2Hat, Qvec);
        for (startFlag = 0; startFlag < 120; startFlag++) {
          d[startFlag] = 1.0 - a * d[startFlag] * (Qvec[startFlag] +
            bineq[startFlag]);
        }

        /*  Update v */
        st.site = &l_emlrtRSI;
        a = 0.0;
        for (k = 0; k < 120; k++) {
          c0 = muDoubleScalarAbs(d[k]);
          if (muDoubleScalarIsNaN(c0) || (c0 > a)) {
            a = c0;
          }
        }

        a = muDoubleScalarMin(1.0, 1.0 / (a * a));
        for (startFlag = 0; startFlag < 120; startFlag++) {
          v[startFlag] += a * d[startFlag];
        }

        (*numIter)++;
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(sp);
        }
      } else {
        exitg1 = 1;
      }
    } while (exitg1 == 0);

    /*  We want at least aa feasible soluton */
  }

  /*  --------------------- MAIN NEWTON ITERATION LOOP --------------------- */
  /*  Then, we finally run the main loop which selects uStar */
  /*  Construct Q matrix */
  mu2 = 1.0;
  while (((*mu > mu_f) || (mu2 > 1.0)) && (*numIter < maxIter)) {
    /*  First, we sample two points of (mu,d) and solve the linear system */
    /*  generated by d = d0 + k*d1, where k = 1/sqrt(mu) */
    for (k = 0; k < 120; k++) {
      d2Hat[k] = muDoubleScalarExp(2.0 * v[k]);
    }

    memset(&SD->f2.Q[0], 0, 14400U * sizeof(real_T));
    for (startFlag = 0; startFlag < 120; startFlag++) {
      SD->f2.Q[startFlag + 120 * startFlag] = d2Hat[startFlag];
    }

    st.site = &m_emlrtRSI;
    mtimes(A, SD->f2.Q, SD->f2.y);
    st.site = &m_emlrtRSI;
    b_mtimes(SD->f2.y, A, SD->f2.Q);

    /*      dG = decomposition(1/2*(G+G'),'chol'); */
    st.site = &n_emlrtRSI;
    for (startFlag = 0; startFlag < 14400; startFlag++) {
      SD->f2.Q[startFlag] += W[startFlag];
    }

    b_st.site = &rb_emlrtRSI;
    cholesky(&b_st, SD->f2.Q);

    /*  Sample point 1 */
    mu1 = *mu;
    st.site = &o_emlrtRSI;
    if (*mu < 0.0) {
      emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    st.site = &p_emlrtRSI;
    solveNewtonStep_decomp(SD, &st, *mu, v, c, A, bineq, SD->f2.Q, Qvec, d1Hat);

    /*  Sample point 2 */
    mu2 = 0.1 * *mu;
    st.site = &q_emlrtRSI;
    k2 = 1.0 / muDoubleScalarSqrt(mu2);
    st.site = &r_emlrtRSI;
    solveNewtonStep_decomp(SD, &st, mu2, v, c, A, bineq, SD->f2.Q, Qvec, d2Hat);

    /*  Obtain affine representation of d = d0 + k*d1 */
    mu2 = 1.0 / muDoubleScalarSqrt(*mu) - k2;
    c0 = -k2 / mu2;
    a = 1.0 / mu2;

    /*  Solve for muStar using bisection */
    for (startFlag = 0; startFlag < 120; startFlag++) {
      mu2 = d1Hat[startFlag];
      k2 = d2Hat[startFlag];
      b_c0[startFlag] = c0 * mu2 + (1.0 - c0) * k2;
      Qvec[startFlag] = a * (mu2 - k2);
    }

    st.site = &s_emlrtRSI;
    muStarSolve(&st, b_c0, Qvec, 1.0, mu_f, mu, d);

    /*  Make sure that muStar is finite and catch it if not */
    if (muDoubleScalarIsInf(*mu)) {
      *mu = mu1;
      memcpy(&d[0], &d1Hat[0], 120U * sizeof(real_T));
    }

    /*  Update x, v, d */
    mu2 = 0.0;
    for (k = 0; k < 120; k++) {
      c0 = muDoubleScalarAbs(d[k]);
      if (muDoubleScalarIsNaN(c0) || (c0 > mu2)) {
        mu2 = c0;
      }
    }

    st.site = &t_emlrtRSI;
    a = muDoubleScalarMin(1.0, 1.0 / (mu2 * mu2));
    for (startFlag = 0; startFlag < 120; startFlag++) {
      v[startFlag] += a * d[startFlag];
    }

    (*numIter)++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  st.site = &u_emlrtRSI;
  *execTime = toc(&st);

  /*  Get x  */
  st.site = &v_emlrtRSI;
  solveNewtonStep(SD, &st, *mu, v, W, c, A, bineq, x);
}

/* End of code generation (logInteriorPoint_rt.c) */
