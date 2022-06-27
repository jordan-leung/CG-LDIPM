/*
 * logInteriorPoint_rt.c
 *
 * Code generation for function 'logInteriorPoint_rt'
 *
 */

/* Include files */
#include "logInteriorPoint_rt.h"
#include "chol.h"
#include "indexShapeCheck.h"
#include "inv.h"
#include "logInteriorPoint_rt_data.h"
#include "logInteriorPoint_rt_emxutil.h"
#include "logInteriorPoint_rt_types.h"
#include "mldivide.h"
#include "mtimes.h"
#include "norm.h"
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

static emlrtRSInfo m_emlrtRSI = { 101, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 111, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 113, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 117, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 118, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 122, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo s_emlrtRSI = { 123, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo t_emlrtRSI = { 132, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 135, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 147, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo w_emlrtRSI = { 151, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo x_emlrtRSI = { 152, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo xb_emlrtRSI = { 34, /* lineNo */
  "chol",                              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/chol.m"/* pathName */
};

static emlrtRSInfo dc_emlrtRSI = { 232,/* lineNo */
  "solveNewtonStep_decomp",            /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo ec_emlrtRSI = { 233,/* lineNo */
  "solveNewtonStep_decomp",            /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo oc_emlrtRSI = { 196,/* lineNo */
  "muStarSolve",                       /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo rc_emlrtRSI = { 212,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo sc_emlrtRSI = { 215,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  104,                                 /* lineNo */
  12,                                  /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  152,                                 /* lineNo */
  25,                                  /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  152,                                 /* lineNo */
  27,                                  /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo emlrtDCI = { 102,   /* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = { 102, /* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m",/* pName */
  4                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  137,                                 /* lineNo */
  5,                                   /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo c_emlrtRTEI = { 13,/* lineNo */
  9,                                   /* colNo */
  "sqrt",                              /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/elfun/sqrt.m"/* pName */
};

static emlrtRTEInfo i_emlrtRTEI = { 102,/* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pName */
};

static emlrtRTEInfo j_emlrtRTEI = { 152,/* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pName */
};

/* Function Declarations */
static void muStarSolve(const emlrtStack *sp, const real_T d0[20], const real_T
  d1[20], real_T mu_f, real_T *muStar, real_T d[20]);

/* Function Definitions */
static void muStarSolve(const emlrtStack *sp, const real_T d0[20], const real_T
  d1[20], real_T mu_f, real_T *muStar, real_T d[20])
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
  while ((!exitg1) && (i < 20)) {
    upper_bound_i = (1.0 - d0[i]) / d1[i];
    lower_bound_i = (-1.0 - d0[i]) / d1[i];

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

  st.site = &oc_emlrtRSI;
  if (*muStar < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  upper_bound_i = 1.0 / muDoubleScalarSqrt(*muStar);
  for (i = 0; i < 20; i++) {
    d[i] = d0[i] + upper_bound_i * d1[i];
  }
}

void logInteriorPoint_rt(const emlrtStack *sp, const real_T W[100], const real_T
  c[10], real_T A[200], const real_T bineq[20], real_T mu_f, real_T mu_0, const
  real_T v0[20], real_T maxIter, const real_T xStar[10], real_T xTol, real_T x
  [10], emxArray_real_T *xError_vec, real_T *execTime, real_T *numIter)
{
  emlrtStack b_st;
  emlrtStack st;
  real_T Q[400];
  real_T y[200];
  real_T y_tmp[200];
  real_T G[100];
  real_T dG[100];
  real_T invW[100];
  real_T Qvec[20];
  real_T b_A[20];
  real_T b_Qvec[20];
  real_T b_a[20];
  real_T c0[20];
  real_T v[20];
  real_T b_c[10];
  real_T b_y[10];
  real_T a;
  real_T d;
  real_T d1;
  real_T k2;
  real_T mu;
  real_T mu1;
  real_T mu2;
  real_T xError;
  int32_T iv[2];
  int32_T exitg1;
  int32_T i;
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
  for (i = 0; i < 200; i++) {
    A[i] = -A[i];
  }

  /*  Pack */
  st.site = &emlrtRSI;
  inv(&st, W, invW);

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
  for (k = 0; k < 20; k++) {
    d = 2.0 * v0[k];
    Qvec[k] = d;
    b_Qvec[k] = muDoubleScalarExp(d);
  }

  memset(&Q[0], 0, 400U * sizeof(real_T));
  for (startFlag = 0; startFlag < 20; startFlag++) {
    Q[startFlag + 20 * startFlag] = b_Qvec[startFlag];
  }

  st.site = &c_emlrtRSI;
  mtimes(A, Q, y);
  st.site = &c_emlrtRSI;
  b_mtimes(y, A, G);
  for (i = 0; i < 100; i++) {
    G[i] += W[i];
  }

  /*  dG = decomposition(G,'lu'); */
  /*  dG = decomposition(1/2*(G+G'),'chol'); */
  st.site = &d_emlrtRSI;
  for (i = 0; i < 10; i++) {
    for (k = 0; k < 10; k++) {
      startFlag = k + 10 * i;
      dG[startFlag] = 0.5 * (G[startFlag] + G[i + 10 * k]);
    }
  }

  b_st.site = &xb_emlrtRSI;
  cholesky(&b_st, dG);

  /*  First, we sample two points of (mu,d) and solve the linear system */
  /*  generated by d = d0 + k*d1, where k = 1/sqrt(mu) */
  /*  Sample point 1 */
  st.site = &e_emlrtRSI;
  if (mu_0 < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  st.site = &f_emlrtRSI;

  /*  Construct Q matrix */
  for (k = 0; k < 20; k++) {
    b_Qvec[k] = muDoubleScalarExp(Qvec[k]);
  }

  memset(&Q[0], 0, 400U * sizeof(real_T));
  for (startFlag = 0; startFlag < 20; startFlag++) {
    Q[startFlag + 20 * startFlag] = b_Qvec[startFlag];
  }

  /*  Solve for x */
  /*  x =  dG\(2*sqrt(mu)*A'*exp(v) - (c + A'*Q*b)); */
  b_st.site = &dc_emlrtRSI;
  if (mu_0 < 0.0) {
    emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  xError = 2.0 * muDoubleScalarSqrt(mu_0);
  for (k = 0; k < 20; k++) {
    b_Qvec[k] = muDoubleScalarExp(v0[k]);
  }

  mtimes(A, Q, y);
  for (i = 0; i < 10; i++) {
    d = 0.0;
    d1 = 0.0;
    for (k = 0; k < 20; k++) {
      d1 += xError * A[k + 20 * i] * b_Qvec[k];
      d += y[i + 10 * k] * bineq[k];
    }

    b_y[i] = d1 - (c[i] + d);
    for (k = 0; k < 10; k++) {
      G[k + 10 * i] = dG[i + 10 * k];
    }
  }

  b_st.site = &dc_emlrtRSI;
  mldivide(&b_st, G, b_y);
  b_st.site = &dc_emlrtRSI;
  mldivide(&b_st, dG, b_y);
  b_st.site = &ec_emlrtRSI;
  if (mu_0 < 0.0) {
    emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  a = 1.0 / muDoubleScalarSqrt(mu_0);
  for (k = 0; k < 20; k++) {
    d = 0.0;
    for (i = 0; i < 10; i++) {
      d += A[k + 20 * i] * b_y[i];
    }

    b_A[k] = 1.0 - a * muDoubleScalarExp(v0[k]) * (d + bineq[k]);
  }

  /*  Sample point 2 */
  mu2 = 0.1 * mu_0;
  st.site = &g_emlrtRSI;
  if (mu2 < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  k2 = 1.0 / muDoubleScalarSqrt(mu2);
  st.site = &h_emlrtRSI;

  /*  Construct Q matrix */
  for (k = 0; k < 20; k++) {
    Qvec[k] = muDoubleScalarExp(Qvec[k]);
  }

  memset(&Q[0], 0, 400U * sizeof(real_T));
  for (startFlag = 0; startFlag < 20; startFlag++) {
    Q[startFlag + 20 * startFlag] = Qvec[startFlag];
  }

  /*  Solve for x */
  /*  x =  dG\(2*sqrt(mu)*A'*exp(v) - (c + A'*Q*b)); */
  b_st.site = &dc_emlrtRSI;
  if (mu2 < 0.0) {
    emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  xError = 2.0 * muDoubleScalarSqrt(mu2);
  for (k = 0; k < 20; k++) {
    b_Qvec[k] = muDoubleScalarExp(v0[k]);
  }

  mtimes(A, Q, y);
  for (i = 0; i < 10; i++) {
    d = 0.0;
    d1 = 0.0;
    for (k = 0; k < 20; k++) {
      d1 += xError * A[k + 20 * i] * b_Qvec[k];
      d += y[i + 10 * k] * bineq[k];
    }

    b_y[i] = d1 - (c[i] + d);
    for (k = 0; k < 10; k++) {
      G[k + 10 * i] = dG[i + 10 * k];
    }
  }

  b_st.site = &dc_emlrtRSI;
  mldivide(&b_st, G, b_y);
  b_st.site = &dc_emlrtRSI;
  mldivide(&b_st, dG, b_y);
  b_st.site = &ec_emlrtRSI;
  if (mu2 < 0.0) {
    emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  a = 1.0 / muDoubleScalarSqrt(mu2);
  for (k = 0; k < 20; k++) {
    d = 0.0;
    for (i = 0; i < 10; i++) {
      d += A[k + 20 * i] * b_y[i];
    }

    b_Qvec[k] = 1.0 - a * muDoubleScalarExp(v0[k]) * (d + bineq[k]);
  }

  /*  Obtain affine representation of d = d0 + k*d1 */
  mu2 = 1.0 / muDoubleScalarSqrt(mu_0) - k2;
  k2 = -k2 / mu2;
  a = 1.0 / mu2;

  /*  Solve for muStar */
  for (i = 0; i < 20; i++) {
    d = b_A[i];
    d1 = b_Qvec[i];
    c0[i] = k2 * d + (1.0 - k2) * d1;
    d = a * (d - d1);
    b_A[i] = d;
  }

  st.site = &i_emlrtRSI;
  muStarSolve(&st, c0, b_A, mu_f, &mu, Qvec);

  /*  If we found a muStar (or a feasible point), then use these values */
  if (!muDoubleScalarIsInf(mu)) {
    startFlag = 1;

    /*  Update v - this updates regardless of whether or not we find a */
    /*  feasible mu... Either just iterates off mu0, or to muStar */
    st.site = &j_emlrtRSI;
    a = 0.0;
    for (k = 0; k < 20; k++) {
      k2 = muDoubleScalarAbs(Qvec[k]);
      if (muDoubleScalarIsNaN(k2) || (k2 > a)) {
        a = k2;
      }
    }

    a = muDoubleScalarMin(1.0, 1.0 / (a * a));
    for (k = 0; k < 20; k++) {
      v[k] = v0[k] + a * Qvec[k];
    }

    *numIter = 1.0;
  } else {
    /*  Otherwise, truly give  up and cold start */
    mu = mu_0;

    /*  under the update at the end */
    for (k = 0; k < 20; k++) {
      Qvec[k] *= rtInf;
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
      xError = 0.0;
      for (k = 0; k < 20; k++) {
        k2 = muDoubleScalarAbs(Qvec[k]);
        if (muDoubleScalarIsNaN(k2) || (k2 > xError)) {
          xError = k2;
        }
      }

      if ((xError > 1.0) && (!(*numIter >= maxIter))) {
        /*  Solve for d */
        st.site = &k_emlrtRSI;

        /*  Construct Q matrix */
        for (k = 0; k < 20; k++) {
          b_Qvec[k] = muDoubleScalarExp(2.0 * v[k]);
        }

        memset(&Q[0], 0, 400U * sizeof(real_T));
        for (startFlag = 0; startFlag < 20; startFlag++) {
          Q[startFlag + 20 * startFlag] = b_Qvec[startFlag];
        }

        /*  Solve for x */
        mtimes(A, Q, y);
        b_mtimes(y, A, G);
        b_st.site = &rc_emlrtRSI;
        if (mu < 0.0) {
          emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
            "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError",
            3, 4, 4, "sqrt");
        }

        xError = 2.0 * muDoubleScalarSqrt(mu);
        for (k = 0; k < 20; k++) {
          b_Qvec[k] = muDoubleScalarExp(v[k]);
        }

        mtimes(A, Q, y);
        for (i = 0; i < 10; i++) {
          d = 0.0;
          d1 = 0.0;
          for (k = 0; k < 20; k++) {
            d1 += xError * A[k + 20 * i] * b_Qvec[k];
            d += y[i + 10 * k] * bineq[k];
          }

          x[i] = d1 - (c[i] + d);
        }

        for (i = 0; i < 100; i++) {
          G[i] += W[i];
        }

        b_st.site = &rc_emlrtRSI;
        mldivide(&b_st, G, x);

        /*  Solve for d */
        b_st.site = &sc_emlrtRSI;
        if (mu < 0.0) {
          emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
            "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError",
            3, 4, 4, "sqrt");
        }

        a = 1.0 / muDoubleScalarSqrt(mu);
        for (k = 0; k < 20; k++) {
          d = 0.0;
          for (i = 0; i < 10; i++) {
            d += A[k + 20 * i] * x[i];
          }

          Qvec[k] = 1.0 - a * muDoubleScalarExp(v[k]) * (d + bineq[k]);
        }

        /*  Update v */
        st.site = &l_emlrtRSI;
        a = 0.0;
        for (k = 0; k < 20; k++) {
          k2 = muDoubleScalarAbs(Qvec[k]);
          if (muDoubleScalarIsNaN(k2) || (k2 > a)) {
            a = k2;
          }
        }

        a = muDoubleScalarMin(1.0, 1.0 / (a * a));
        for (i = 0; i < 20; i++) {
          v[i] += a * Qvec[i];
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
  for (k = 0; k < 20; k++) {
    b_Qvec[k] = muDoubleScalarExp(v[k]);
  }

  st.site = &m_emlrtRSI;
  if (mu < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  xError = muDoubleScalarSqrt(mu);
  for (i = 0; i < 20; i++) {
    for (k = 0; k < 10; k++) {
      y_tmp[k + 10 * i] = A[i + 20 * k];
    }

    d = b_Qvec[i];
    c0[i] = d + d * Qvec[i];
  }

  for (i = 0; i < 10; i++) {
    d = 0.0;
    for (k = 0; k < 20; k++) {
      d += xError * y_tmp[i + 10 * k] * c0[k];
    }

    b_y[i] = d - c[i];
  }

  for (i = 0; i < 10; i++) {
    d = 0.0;
    for (k = 0; k < 10; k++) {
      d += invW[i + 10 * k] * b_y[k];
    }

    x[i] = d;
  }

  if (!(maxIter >= 0.0)) {
    emlrtNonNegativeCheckR2012b(maxIter, &b_emlrtDCI, sp);
  }

  d = (int32_T)muDoubleScalarFloor(maxIter);
  if (maxIter != d) {
    emlrtIntegerCheckR2012b(maxIter, &emlrtDCI, sp);
  }

  i = xError_vec->size[0];
  xError_vec->size[0] = (int32_T)maxIter;
  emxEnsureCapacity_real_T(sp, xError_vec, i, &i_emlrtRTEI);
  if (maxIter != d) {
    emlrtIntegerCheckR2012b(maxIter, &emlrtDCI, sp);
  }

  startFlag = (int32_T)maxIter;
  for (i = 0; i < startFlag; i++) {
    xError_vec->data[i] = 0.0;
  }

  for (i = 0; i < 10; i++) {
    b_c[i] = x[i] - xStar[i];
  }

  xError = b_norm(b_c);
  if (((int32_T)*numIter < 1) || ((int32_T)*numIter > (int32_T)maxIter)) {
    emlrtDynamicBoundsCheckR2012b((int32_T)*numIter, 1, (int32_T)maxIter,
      &emlrtBCI, sp);
  }

  xError_vec->data[(int32_T)*numIter - 1] = xError;
  mu2 = 1.0;
  while (((xError > xTol) || (mu2 > 1.0)) && (*numIter < maxIter)) {
    /*  First, we sample two points of (mu,d) and solve the linear system */
    /*  generated by d = d0 + k*d1, where k = 1/sqrt(mu) */
    for (k = 0; k < 20; k++) {
      d = 2.0 * v[k];
      Qvec[k] = d;
      b_Qvec[k] = muDoubleScalarExp(d);
    }

    memset(&Q[0], 0, 400U * sizeof(real_T));
    for (startFlag = 0; startFlag < 20; startFlag++) {
      Q[startFlag + 20 * startFlag] = b_Qvec[startFlag];
    }

    st.site = &n_emlrtRSI;
    mtimes(A, Q, y);
    st.site = &n_emlrtRSI;
    b_mtimes(y, A, G);
    for (i = 0; i < 100; i++) {
      G[i] += W[i];
    }

    /*      dG = decomposition(1/2*(G+G'),'chol'); */
    st.site = &o_emlrtRSI;
    for (i = 0; i < 10; i++) {
      for (k = 0; k < 10; k++) {
        startFlag = k + 10 * i;
        dG[startFlag] = 0.5 * (G[startFlag] + G[i + 10 * k]);
      }
    }

    b_st.site = &xb_emlrtRSI;
    cholesky(&b_st, dG);

    /*  Sample point 1 */
    mu1 = mu;
    st.site = &p_emlrtRSI;
    if (mu < 0.0) {
      emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    st.site = &q_emlrtRSI;

    /*  Construct Q matrix */
    for (k = 0; k < 20; k++) {
      b_Qvec[k] = muDoubleScalarExp(Qvec[k]);
    }

    memset(&Q[0], 0, 400U * sizeof(real_T));
    for (startFlag = 0; startFlag < 20; startFlag++) {
      Q[startFlag + 20 * startFlag] = b_Qvec[startFlag];
    }

    /*  Solve for x */
    /*  x =  dG\(2*sqrt(mu)*A'*exp(v) - (c + A'*Q*b)); */
    b_st.site = &dc_emlrtRSI;
    if (mu < 0.0) {
      emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    xError = 2.0 * muDoubleScalarSqrt(mu);
    for (k = 0; k < 20; k++) {
      b_Qvec[k] = muDoubleScalarExp(v[k]);
    }

    mtimes(A, Q, y);
    for (i = 0; i < 10; i++) {
      d = 0.0;
      d1 = 0.0;
      for (k = 0; k < 20; k++) {
        d1 += xError * A[k + 20 * i] * b_Qvec[k];
        d += y[i + 10 * k] * bineq[k];
      }

      b_y[i] = d1 - (c[i] + d);
      for (k = 0; k < 10; k++) {
        G[k + 10 * i] = dG[i + 10 * k];
      }
    }

    b_st.site = &dc_emlrtRSI;
    mldivide(&b_st, G, b_y);
    b_st.site = &dc_emlrtRSI;
    mldivide(&b_st, dG, b_y);
    b_st.site = &ec_emlrtRSI;
    if (mu < 0.0) {
      emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    a = 1.0 / muDoubleScalarSqrt(mu);
    for (k = 0; k < 20; k++) {
      d = 0.0;
      for (i = 0; i < 10; i++) {
        d += A[k + 20 * i] * b_y[i];
      }

      b_A[k] = 1.0 - a * muDoubleScalarExp(v[k]) * (d + bineq[k]);
    }

    /*  Sample point 2 */
    mu2 = 0.1 * mu;
    st.site = &r_emlrtRSI;
    if (mu2 < 0.0) {
      emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    k2 = 1.0 / muDoubleScalarSqrt(mu2);
    st.site = &s_emlrtRSI;

    /*  Construct Q matrix */
    for (k = 0; k < 20; k++) {
      Qvec[k] = muDoubleScalarExp(Qvec[k]);
    }

    memset(&Q[0], 0, 400U * sizeof(real_T));
    for (startFlag = 0; startFlag < 20; startFlag++) {
      Q[startFlag + 20 * startFlag] = Qvec[startFlag];
    }

    /*  Solve for x */
    /*  x =  dG\(2*sqrt(mu)*A'*exp(v) - (c + A'*Q*b)); */
    b_st.site = &dc_emlrtRSI;
    if (mu2 < 0.0) {
      emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    xError = 2.0 * muDoubleScalarSqrt(mu2);
    for (k = 0; k < 20; k++) {
      b_Qvec[k] = muDoubleScalarExp(v[k]);
    }

    mtimes(A, Q, y);
    for (i = 0; i < 10; i++) {
      d = 0.0;
      d1 = 0.0;
      for (k = 0; k < 20; k++) {
        d1 += xError * A[k + 20 * i] * b_Qvec[k];
        d += y[i + 10 * k] * bineq[k];
      }

      b_y[i] = d1 - (c[i] + d);
      for (k = 0; k < 10; k++) {
        G[k + 10 * i] = dG[i + 10 * k];
      }
    }

    b_st.site = &dc_emlrtRSI;
    mldivide(&b_st, G, b_y);
    b_st.site = &dc_emlrtRSI;
    mldivide(&b_st, dG, b_y);
    b_st.site = &ec_emlrtRSI;
    if (mu2 < 0.0) {
      emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    a = 1.0 / muDoubleScalarSqrt(mu2);
    for (k = 0; k < 20; k++) {
      d = 0.0;
      for (i = 0; i < 10; i++) {
        d += A[k + 20 * i] * b_y[i];
      }

      b_Qvec[k] = 1.0 - a * muDoubleScalarExp(v[k]) * (d + bineq[k]);
    }

    /*  Obtain affine representation of d = d0 + k*d1 */
    mu2 = 1.0 / muDoubleScalarSqrt(mu) - k2;
    k2 = -k2 / mu2;
    a = 1.0 / mu2;

    /*  Solve for muStar using bisection */
    for (i = 0; i < 20; i++) {
      d = b_A[i];
      d1 = b_Qvec[i];
      c0[i] = k2 * d + (1.0 - k2) * d1;
      b_a[i] = a * (d - d1);
    }

    st.site = &t_emlrtRSI;
    muStarSolve(&st, c0, b_a, mu_f, &mu, Qvec);

    /*  Get xError */
    for (k = 0; k < 20; k++) {
      b_Qvec[k] = muDoubleScalarExp(v[k]);
    }

    st.site = &u_emlrtRSI;
    if (mu < 0.0) {
      emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    xError = muDoubleScalarSqrt(mu);
    for (i = 0; i < 20; i++) {
      d = b_Qvec[i];
      c0[i] = d + d * Qvec[i];
    }

    for (i = 0; i < 10; i++) {
      d = 0.0;
      for (k = 0; k < 20; k++) {
        d += xError * y_tmp[i + 10 * k] * c0[k];
      }

      b_y[i] = d - c[i];
    }

    for (i = 0; i < 10; i++) {
      d = 0.0;
      for (k = 0; k < 10; k++) {
        d += invW[i + 10 * k] * b_y[k];
      }

      x[i] = d;
      b_c[i] = d - xStar[i];
    }

    xError = b_norm(b_c);
    if (((int32_T)*numIter < 1) || ((int32_T)*numIter > xError_vec->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)*numIter, 1, xError_vec->size[0],
        &d_emlrtBCI, sp);
    }

    xError_vec->data[(int32_T)*numIter - 1] = xError;

    /*  Make sure that muStar is finite and catch it if not */
    if (muDoubleScalarIsInf(mu)) {
      mu = mu1;
      memcpy(&Qvec[0], &b_A[0], 20U * sizeof(real_T));
    }

    /*  Update x, v, d */
    mu2 = 0.0;
    for (k = 0; k < 20; k++) {
      k2 = muDoubleScalarAbs(Qvec[k]);
      if (muDoubleScalarIsNaN(k2) || (k2 > mu2)) {
        mu2 = k2;
      }
    }

    st.site = &v_emlrtRSI;
    a = muDoubleScalarMin(1.0, 1.0 / (mu2 * mu2));
    for (i = 0; i < 20; i++) {
      v[i] += a * Qvec[i];
    }

    (*numIter)++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  st.site = &w_emlrtRSI;
  *execTime = toc(&st);
  if (1.0 > *numIter - 1.0) {
    i = 0;
  } else {
    if (1 > xError_vec->size[0]) {
      emlrtDynamicBoundsCheckR2012b(1, 1, xError_vec->size[0], &b_emlrtBCI, sp);
    }

    if (((int32_T)((uint32_T)*numIter - 1U) < 1) || ((int32_T)((uint32_T)
          *numIter - 1U) > xError_vec->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)*numIter - 1U), 1,
        xError_vec->size[0], &c_emlrtBCI, sp);
    }

    i = (int32_T)((uint32_T)*numIter - 1U);
  }

  iv[0] = 1;
  iv[1] = i;
  st.site = &x_emlrtRSI;
  indexShapeCheck(&st, xError_vec->size[0], iv);
  k = xError_vec->size[0];
  xError_vec->size[0] = i;
  emxEnsureCapacity_real_T(sp, xError_vec, k, &j_emlrtRTEI);
}

/* End of code generation (logInteriorPoint_rt.c) */
