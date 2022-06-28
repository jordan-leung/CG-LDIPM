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

static emlrtRSInfo k_emlrtRSI = { 77,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 92,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 97,  /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 113, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 115, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 119, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 120, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 124, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo s_emlrtRSI = { 125, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo t_emlrtRSI = { 134, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 137, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 149, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo w_emlrtRSI = { 153, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo x_emlrtRSI = { 154, /* lineNo */
  "logInteriorPoint_rt",               /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo rb_emlrtRSI = { 34, /* lineNo */
  "chol",                              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/matfun/chol.m"/* pathName */
};

static emlrtRSInfo vb_emlrtRSI = { 234,/* lineNo */
  "solveNewtonStep_decomp",            /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo wb_emlrtRSI = { 235,/* lineNo */
  "solveNewtonStep_decomp",            /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo ic_emlrtRSI = { 198,/* lineNo */
  "muStarSolve",                       /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo lc_emlrtRSI = { 214,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtRSInfo mc_emlrtRSI = { 217,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pathName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  80,                                  /* lineNo */
  12,                                  /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  154,                                 /* lineNo */
  25,                                  /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  154,                                 /* lineNo */
  27,                                  /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo emlrtDCI = { 78,    /* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = { 78,  /* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m",/* pName */
  4                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  94,                                  /* lineNo */
  9,                                   /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  139,                                 /* lineNo */
  5,                                   /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo d_emlrtRTEI = { 13,/* lineNo */
  9,                                   /* colNo */
  "sqrt",                              /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/elfun/sqrt.m"/* pName */
};

static emlrtRTEInfo h_emlrtRTEI = { 78,/* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pName */
};

static emlrtRTEInfo i_emlrtRTEI = { 154,/* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_rt",               /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_rt.m"/* pName */
};

/* Function Declarations */
static void muStarSolve(const emlrtStack *sp, const real_T d0[50], const real_T
  d1[50], real_T mu_f, real_T *muStar, real_T d[50]);

/* Function Definitions */
static void muStarSolve(const emlrtStack *sp, const real_T d0[50], const real_T
  d1[50], real_T mu_f, real_T *muStar, real_T d[50])
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
  while ((!exitg1) && (i < 50)) {
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

  st.site = &ic_emlrtRSI;
  if (*muStar < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  upper_bound_i = 1.0 / muDoubleScalarSqrt(*muStar);
  for (i = 0; i < 50; i++) {
    d[i] = d0[i] + upper_bound_i * d1[i];
  }
}

void logInteriorPoint_rt(logInteriorPoint_rtStackData *SD, const emlrtStack *sp,
  const real_T W[10000], const real_T c[100], real_T A[5000], const real_T
  bineq[50], real_T mu_f, real_T mu_0, const real_T v0[50], real_T maxIter,
  const real_T xStar[100], real_T xTol, real_T x[100], emxArray_real_T
  *xError_vec, real_T *execTime, real_T *numIter)
{
  emlrtStack b_st;
  emlrtStack st;
  real_T y[5000];
  real_T y_tmp[5000];
  real_T Q[2500];
  real_T b_y[100];
  real_T c_y[100];
  real_T Qvec[50];
  real_T b[50];
  real_T b_Qvec[50];
  real_T c0[50];
  real_T d[50];
  real_T d_y[50];
  real_T v[50];
  real_T a;
  real_T k2;
  real_T mu;
  real_T mu1;
  real_T mu2;
  real_T xError;
  int32_T iv[2];
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
  for (i = 0; i < 5000; i++) {
    A[i] = -A[i];
  }

  /*  Pack */
  st.site = &emlrtRSI;
  inv(&st, W, SD->f0.invW);

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
  for (k = 0; k < 50; k++) {
    xError = 2.0 * v0[k];
    Qvec[k] = xError;
    b_Qvec[k] = muDoubleScalarExp(xError);
  }

  memset(&Q[0], 0, 2500U * sizeof(real_T));
  for (k = 0; k < 50; k++) {
    Q[k + 50 * k] = b_Qvec[k];
  }

  st.site = &c_emlrtRSI;
  mtimes(A, Q, y);
  st.site = &c_emlrtRSI;
  b_mtimes(y, A, SD->f0.G);

  /*  dG = decomposition(G,'lu'); */
  /*  dG = decomposition(1/2*(G+G'),'chol'); */
  st.site = &d_emlrtRSI;
  for (i = 0; i < 10000; i++) {
    SD->f0.G[i] += W[i];
  }

  b_st.site = &rb_emlrtRSI;
  cholesky(&b_st, SD->f0.G);

  /*  First, we sample two points of (mu,d) and solve the linear system */
  /*  generated by d = d0 + k*d1, where k = 1/sqrt(mu) */
  /*  Sample point 1 */
  st.site = &e_emlrtRSI;
  if (mu_0 < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  st.site = &f_emlrtRSI;

  /*  Construct Q matrix */
  for (k = 0; k < 50; k++) {
    b_Qvec[k] = muDoubleScalarExp(Qvec[k]);
  }

  memset(&Q[0], 0, 2500U * sizeof(real_T));
  for (k = 0; k < 50; k++) {
    Q[k + 50 * k] = b_Qvec[k];
  }

  /*  Solve for x */
  /*  x =  dG\(2*sqrt(mu)*A'*exp(v) - (c + A'*Q*b)); */
  b_st.site = &vb_emlrtRSI;
  if (mu_0 < 0.0) {
    emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  xError = 2.0 * muDoubleScalarSqrt(mu_0);
  for (k = 0; k < 50; k++) {
    b[k] = muDoubleScalarExp(v0[k]);
    for (i = 0; i < 100; i++) {
      y[i + 100 * k] = xError * A[k + 50 * i];
    }
  }

  c_mtimes(y, b, b_y);
  mtimes(A, Q, y);
  c_mtimes(y, bineq, c_y);
  for (i = 0; i < 100; i++) {
    b_y[i] -= c[i] + c_y[i];
    for (k = 0; k < 100; k++) {
      SD->f0.b_G[k + 100 * i] = SD->f0.G[i + 100 * k];
    }
  }

  b_st.site = &vb_emlrtRSI;
  mldivide(&b_st, SD->f0.b_G, b_y);
  b_st.site = &vb_emlrtRSI;
  mldivide(&b_st, SD->f0.G, b_y);
  b_st.site = &wb_emlrtRSI;
  if (mu_0 < 0.0) {
    emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  a = 1.0 / muDoubleScalarSqrt(mu_0);
  for (k = 0; k < 50; k++) {
    b[k] = muDoubleScalarExp(v0[k]);
  }

  d_mtimes(A, b_y, d_y);
  for (i = 0; i < 50; i++) {
    b[i] = 1.0 - a * b[i] * (d_y[i] + bineq[i]);
  }

  /*  Sample point 2 */
  mu2 = 0.1 * mu_0;
  st.site = &g_emlrtRSI;
  if (mu2 < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  k2 = 1.0 / muDoubleScalarSqrt(mu2);
  st.site = &h_emlrtRSI;

  /*  Construct Q matrix */
  for (k = 0; k < 50; k++) {
    Qvec[k] = muDoubleScalarExp(Qvec[k]);
  }

  memset(&Q[0], 0, 2500U * sizeof(real_T));
  for (k = 0; k < 50; k++) {
    Q[k + 50 * k] = Qvec[k];
  }

  /*  Solve for x */
  /*  x =  dG\(2*sqrt(mu)*A'*exp(v) - (c + A'*Q*b)); */
  b_st.site = &vb_emlrtRSI;
  if (mu2 < 0.0) {
    emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  xError = 2.0 * muDoubleScalarSqrt(mu2);
  for (k = 0; k < 50; k++) {
    b_Qvec[k] = muDoubleScalarExp(v0[k]);
    for (i = 0; i < 100; i++) {
      y[i + 100 * k] = xError * A[k + 50 * i];
    }
  }

  c_mtimes(y, b_Qvec, b_y);
  mtimes(A, Q, y);
  c_mtimes(y, bineq, c_y);
  for (i = 0; i < 100; i++) {
    b_y[i] -= c[i] + c_y[i];
    for (k = 0; k < 100; k++) {
      SD->f0.b_G[k + 100 * i] = SD->f0.G[i + 100 * k];
    }
  }

  b_st.site = &vb_emlrtRSI;
  mldivide(&b_st, SD->f0.b_G, b_y);
  b_st.site = &vb_emlrtRSI;
  mldivide(&b_st, SD->f0.G, b_y);
  b_st.site = &wb_emlrtRSI;
  if (mu2 < 0.0) {
    emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  a = 1.0 / muDoubleScalarSqrt(mu2);
  for (k = 0; k < 50; k++) {
    b_Qvec[k] = muDoubleScalarExp(v0[k]);
  }

  d_mtimes(A, b_y, d_y);
  for (i = 0; i < 50; i++) {
    b_Qvec[i] = 1.0 - a * b_Qvec[i] * (d_y[i] + bineq[i]);
  }

  /*  Obtain affine representation of d = d0 + k*d1 */
  mu2 = 1.0 / muDoubleScalarSqrt(mu_0) - k2;
  k2 = -k2 / mu2;
  a = 1.0 / mu2;

  /*  Solve for muStar */
  for (i = 0; i < 50; i++) {
    xError = b[i];
    mu2 = b_Qvec[i];
    c0[i] = k2 * xError + (1.0 - k2) * mu2;
    d_y[i] = a * (xError - mu2);
  }

  st.site = &i_emlrtRSI;
  muStarSolve(&st, c0, d_y, mu_f, &mu, d);

  /*  If we found a muStar (or a feasible point), then use these values */
  if (!muDoubleScalarIsInf(mu)) {
    startFlag = 1;

    /*  Update v - this updates regardless of whether or not we find a */
    /*  feasible mu... Either just iterates off mu0, or to muStar */
    st.site = &j_emlrtRSI;
    a = b_norm(d);
    a = muDoubleScalarMin(1.0, 1.0 / (a * a));
    for (k = 0; k < 50; k++) {
      v[k] = v0[k] + a * d[k];
    }

    *numIter = 1.0;
  } else {
    /*  Otherwise, truly give  up and cold start */
    mu = mu_0;

    /*  under the update at the end */
    for (k = 0; k < 50; k++) {
      d[k] *= rtInf;
      v[k] = 0.0;
    }

    startFlag = 0;
  }

  st.site = &k_emlrtRSI;
  if (mu < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  xError = muDoubleScalarSqrt(mu);
  for (i = 0; i < 50; i++) {
    for (k = 0; k < 100; k++) {
      y_tmp[k + 100 * i] = A[i + 50 * k];
    }
  }

  st.site = &k_emlrtRSI;
  for (k = 0; k < 50; k++) {
    mu2 = muDoubleScalarExp(v[k]);
    b_Qvec[k] = mu2;
    Qvec[k] = mu2;
  }

  for (i = 0; i < 5000; i++) {
    y[i] = xError * y_tmp[i];
  }

  for (i = 0; i < 50; i++) {
    c0[i] = b_Qvec[i] + Qvec[i] * d[i];
  }

  c_mtimes(y, c0, b_y);
  st.site = &k_emlrtRSI;
  for (i = 0; i < 100; i++) {
    b_y[i] -= c[i];
  }

  e_mtimes(SD->f0.invW, b_y, x);
  if (!(maxIter >= 0.0)) {
    emlrtNonNegativeCheckR2012b(maxIter, &b_emlrtDCI, sp);
  }

  xError = (int32_T)muDoubleScalarFloor(maxIter);
  if (maxIter != xError) {
    emlrtIntegerCheckR2012b(maxIter, &emlrtDCI, sp);
  }

  i = xError_vec->size[0];
  xError_vec->size[0] = (int32_T)maxIter;
  emxEnsureCapacity_real_T(sp, xError_vec, i, &h_emlrtRTEI);
  if (maxIter != xError) {
    emlrtIntegerCheckR2012b(maxIter, &emlrtDCI, sp);
  }

  k = (int32_T)maxIter;
  for (i = 0; i < k; i++) {
    xError_vec->data[i] = 0.0;
  }

  for (i = 0; i < 100; i++) {
    c_y[i] = x[i] - xStar[i];
  }

  xError = c_norm(c_y);
  if (1 > (int32_T)maxIter) {
    emlrtDynamicBoundsCheckR2012b(1, 1, (int32_T)maxIter, &emlrtBCI, sp);
  }

  xError_vec->data[0] = xError;

  /*  If we have a pair (mu,d) with d > 1 and we iterate the Newton */
  /*  algorithm until we obtain convergence. If we found a muStar from our */
  /*  initialization, this loop will be skipped. */
  if (startFlag < 1) {
    while ((b_norm(d) > 1.0) && (!(*numIter >= maxIter))) {
      /*  Solve for d */
      st.site = &l_emlrtRSI;

      /*  Construct Q matrix */
      for (k = 0; k < 50; k++) {
        b_Qvec[k] = muDoubleScalarExp(2.0 * v[k]);
      }

      memset(&Q[0], 0, 2500U * sizeof(real_T));
      for (k = 0; k < 50; k++) {
        Q[k + 50 * k] = b_Qvec[k];
      }

      /*  Solve for x */
      mtimes(A, Q, y);
      b_mtimes(y, A, SD->f0.G);
      b_st.site = &lc_emlrtRSI;
      if (mu < 0.0) {
        emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
          "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3,
          4, 4, "sqrt");
      }

      xError = 2.0 * muDoubleScalarSqrt(mu);
      for (k = 0; k < 50; k++) {
        b[k] = muDoubleScalarExp(v[k]);
        for (i = 0; i < 100; i++) {
          y[i + 100 * k] = xError * A[k + 50 * i];
        }
      }

      c_mtimes(y, b, x);
      mtimes(A, Q, y);
      c_mtimes(y, bineq, b_y);
      for (i = 0; i < 100; i++) {
        x[i] -= c[i] + b_y[i];
      }

      for (i = 0; i < 10000; i++) {
        SD->f0.G[i] += W[i];
      }

      b_st.site = &lc_emlrtRSI;
      mldivide(&b_st, SD->f0.G, x);

      /*  Solve for d */
      b_st.site = &mc_emlrtRSI;
      if (mu < 0.0) {
        emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
          "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3,
          4, 4, "sqrt");
      }

      a = 1.0 / muDoubleScalarSqrt(mu);
      for (k = 0; k < 50; k++) {
        d[k] = muDoubleScalarExp(v[k]);
      }

      d_mtimes(A, x, d_y);
      for (i = 0; i < 50; i++) {
        d[i] = 1.0 - a * d[i] * (d_y[i] + bineq[i]);
      }

      for (i = 0; i < 100; i++) {
        c_y[i] = x[i] - xStar[i];
      }

      xError = c_norm(c_y);
      if (((int32_T)(*numIter + 1.0) < 1) || ((int32_T)(*numIter + 1.0) >
           xError_vec->size[0])) {
        emlrtDynamicBoundsCheckR2012b((int32_T)(*numIter + 1.0), 1,
          xError_vec->size[0], &d_emlrtBCI, sp);
      }

      xError_vec->data[(int32_T)(*numIter + 1.0) - 1] = xError;

      /*  Update v */
      st.site = &m_emlrtRSI;
      a = b_norm(d);
      a = muDoubleScalarMin(1.0, 1.0 / (a * a));
      for (i = 0; i < 50; i++) {
        v[i] += a * d[i];
      }

      (*numIter)++;
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }

    /*  We want at least aa feasible soluton */
  }

  /*  --------------------- MAIN NEWTON ITERATION LOOP --------------------- */
  /*  Then, we finally run the main loop which selects uStar */
  /*  Construct Q matrix */
  mu2 = 1.0;
  while (((xError > xTol) || (mu2 > 1.0)) && (*numIter < maxIter)) {
    /*  First, we sample two points of (mu,d) and solve the linear system */
    /*  generated by d = d0 + k*d1, where k = 1/sqrt(mu) */
    for (k = 0; k < 50; k++) {
      xError = 2.0 * v[k];
      Qvec[k] = xError;
      d_y[k] = muDoubleScalarExp(xError);
    }

    memset(&Q[0], 0, 2500U * sizeof(real_T));
    for (k = 0; k < 50; k++) {
      Q[k + 50 * k] = d_y[k];
    }

    st.site = &n_emlrtRSI;
    mtimes(A, Q, y);
    st.site = &n_emlrtRSI;
    b_mtimes(y, A, SD->f0.G);

    /*      dG = decomposition(1/2*(G+G'),'chol'); */
    st.site = &o_emlrtRSI;
    for (i = 0; i < 10000; i++) {
      SD->f0.G[i] += W[i];
    }

    b_st.site = &rb_emlrtRSI;
    cholesky(&b_st, SD->f0.G);

    /*  Sample point 1 */
    mu1 = mu;
    st.site = &p_emlrtRSI;
    if (mu < 0.0) {
      emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    st.site = &q_emlrtRSI;

    /*  Construct Q matrix */
    for (k = 0; k < 50; k++) {
      b_Qvec[k] = muDoubleScalarExp(Qvec[k]);
    }

    memset(&Q[0], 0, 2500U * sizeof(real_T));
    for (k = 0; k < 50; k++) {
      Q[k + 50 * k] = b_Qvec[k];
    }

    /*  Solve for x */
    /*  x =  dG\(2*sqrt(mu)*A'*exp(v) - (c + A'*Q*b)); */
    b_st.site = &vb_emlrtRSI;
    if (mu < 0.0) {
      emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    xError = 2.0 * muDoubleScalarSqrt(mu);
    for (k = 0; k < 50; k++) {
      b[k] = muDoubleScalarExp(v[k]);
      for (i = 0; i < 100; i++) {
        y[i + 100 * k] = xError * A[k + 50 * i];
      }
    }

    c_mtimes(y, b, b_y);
    mtimes(A, Q, y);
    c_mtimes(y, bineq, c_y);
    for (i = 0; i < 100; i++) {
      b_y[i] -= c[i] + c_y[i];
      for (k = 0; k < 100; k++) {
        SD->f0.b_G[k + 100 * i] = SD->f0.G[i + 100 * k];
      }
    }

    b_st.site = &vb_emlrtRSI;
    mldivide(&b_st, SD->f0.b_G, b_y);
    b_st.site = &vb_emlrtRSI;
    mldivide(&b_st, SD->f0.G, b_y);
    b_st.site = &wb_emlrtRSI;
    if (mu < 0.0) {
      emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    a = 1.0 / muDoubleScalarSqrt(mu);
    for (k = 0; k < 50; k++) {
      b[k] = muDoubleScalarExp(v[k]);
    }

    d_mtimes(A, b_y, d_y);
    for (i = 0; i < 50; i++) {
      b[i] = 1.0 - a * b[i] * (d_y[i] + bineq[i]);
    }

    /*  Sample point 2 */
    mu2 = 0.1 * mu;
    st.site = &r_emlrtRSI;
    if (mu2 < 0.0) {
      emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    k2 = 1.0 / muDoubleScalarSqrt(mu2);
    st.site = &s_emlrtRSI;

    /*  Construct Q matrix */
    for (k = 0; k < 50; k++) {
      Qvec[k] = muDoubleScalarExp(Qvec[k]);
    }

    memset(&Q[0], 0, 2500U * sizeof(real_T));
    for (k = 0; k < 50; k++) {
      Q[k + 50 * k] = Qvec[k];
    }

    /*  Solve for x */
    /*  x =  dG\(2*sqrt(mu)*A'*exp(v) - (c + A'*Q*b)); */
    b_st.site = &vb_emlrtRSI;
    if (mu2 < 0.0) {
      emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    xError = 2.0 * muDoubleScalarSqrt(mu2);
    for (k = 0; k < 50; k++) {
      b_Qvec[k] = muDoubleScalarExp(v[k]);
      for (i = 0; i < 100; i++) {
        y[i + 100 * k] = xError * A[k + 50 * i];
      }
    }

    c_mtimes(y, b_Qvec, b_y);
    mtimes(A, Q, y);
    c_mtimes(y, bineq, c_y);
    for (i = 0; i < 100; i++) {
      b_y[i] -= c[i] + c_y[i];
      for (k = 0; k < 100; k++) {
        SD->f0.b_G[k + 100 * i] = SD->f0.G[i + 100 * k];
      }
    }

    b_st.site = &vb_emlrtRSI;
    mldivide(&b_st, SD->f0.b_G, b_y);
    b_st.site = &vb_emlrtRSI;
    mldivide(&b_st, SD->f0.G, b_y);
    b_st.site = &wb_emlrtRSI;
    if (mu2 < 0.0) {
      emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    a = 1.0 / muDoubleScalarSqrt(mu2);
    for (k = 0; k < 50; k++) {
      b_Qvec[k] = muDoubleScalarExp(v[k]);
    }

    d_mtimes(A, b_y, d_y);
    for (i = 0; i < 50; i++) {
      b_Qvec[i] = 1.0 - a * b_Qvec[i] * (d_y[i] + bineq[i]);
    }

    /*  Obtain affine representation of d = d0 + k*d1 */
    mu2 = 1.0 / muDoubleScalarSqrt(mu) - k2;
    k2 = -k2 / mu2;
    a = 1.0 / mu2;

    /*  Solve for muStar using bisection */
    for (i = 0; i < 50; i++) {
      xError = b[i];
      mu2 = b_Qvec[i];
      c0[i] = k2 * xError + (1.0 - k2) * mu2;
      d_y[i] = a * (xError - mu2);
    }

    st.site = &t_emlrtRSI;
    muStarSolve(&st, c0, d_y, mu_f, &mu, d);

    /*  Get xError */
    st.site = &u_emlrtRSI;
    if (mu < 0.0) {
      emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    xError = muDoubleScalarSqrt(mu);
    st.site = &u_emlrtRSI;
    for (k = 0; k < 50; k++) {
      mu2 = muDoubleScalarExp(v[k]);
      b_Qvec[k] = mu2;
      Qvec[k] = mu2;
    }

    for (i = 0; i < 5000; i++) {
      y[i] = xError * y_tmp[i];
    }

    for (i = 0; i < 50; i++) {
      c0[i] = b_Qvec[i] + Qvec[i] * d[i];
    }

    c_mtimes(y, c0, b_y);
    st.site = &u_emlrtRSI;
    for (i = 0; i < 100; i++) {
      b_y[i] -= c[i];
    }

    e_mtimes(SD->f0.invW, b_y, x);
    for (i = 0; i < 100; i++) {
      c_y[i] = x[i] - xStar[i];
    }

    xError = c_norm(c_y);
    if (((int32_T)(*numIter + 1.0) < 1) || ((int32_T)(*numIter + 1.0) >
         xError_vec->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(*numIter + 1.0), 1,
        xError_vec->size[0], &e_emlrtBCI, sp);
    }

    xError_vec->data[(int32_T)(*numIter + 1.0) - 1] = xError;

    /*  Make sure that muStar is finite and catch it if not */
    if (muDoubleScalarIsInf(mu)) {
      mu = mu1;
      memcpy(&d[0], &b[0], 50U * sizeof(real_T));
    }

    /*  Update x, v, d */
    mu2 = b_norm(d);
    st.site = &v_emlrtRSI;
    a = muDoubleScalarMin(1.0, 1.0 / (mu2 * mu2));
    for (i = 0; i < 50; i++) {
      v[i] += a * d[i];
    }

    (*numIter)++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  st.site = &w_emlrtRSI;
  *execTime = toc(&st);
  if (1.0 > *numIter) {
    i = 0;
  } else {
    if (1 > xError_vec->size[0]) {
      emlrtDynamicBoundsCheckR2012b(1, 1, xError_vec->size[0], &b_emlrtBCI, sp);
    }

    if (((int32_T)*numIter < 1) || ((int32_T)*numIter > xError_vec->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)*numIter, 1, xError_vec->size[0],
        &c_emlrtBCI, sp);
    }

    i = (int32_T)*numIter;
  }

  iv[0] = 1;
  iv[1] = i;
  st.site = &x_emlrtRSI;
  indexShapeCheck(&st, xError_vec->size[0], iv);
  k = xError_vec->size[0];
  xError_vec->size[0] = i;
  emxEnsureCapacity_real_T(sp, xError_vec, k, &i_emlrtRTEI);
}

/* End of code generation (logInteriorPoint_rt.c) */
