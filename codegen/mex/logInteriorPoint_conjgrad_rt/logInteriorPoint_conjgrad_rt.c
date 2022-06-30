/*
 * logInteriorPoint_conjgrad_rt.c
 *
 * Code generation for function 'logInteriorPoint_conjgrad_rt'
 *
 */

/* Include files */
#include "logInteriorPoint_conjgrad_rt.h"
#include "find.h"
#include "inv.h"
#include "logInteriorPoint_conjgrad_rt_data.h"
#include "logInteriorPoint_conjgrad_rt_types.h"
#include "mtimes.h"
#include "norm.h"
#include "rt_nonfinite.h"
#include "tic.h"
#include "toc.h"
#include "mwmathutil.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 12,    /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 28,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 37,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 55,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 59,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 69,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 79,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 105, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 109, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 126, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 130, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 140, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 144, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 152, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 157, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 159, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo kb_emlrtRSI = { 192,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo ec_emlrtRSI = { 324,/* lineNo */
  "solveNewtonStep_warmStart",         /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo gc_emlrtRSI = { 347,/* lineNo */
  "solveNewtonStep_warmStart",         /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo jc_emlrtRSI = { 412,/* lineNo */
  "muStarSolve",                       /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRTEInfo d_emlrtRTEI = { 13,/* lineNo */
  9,                                   /* colNo */
  "sqrt",                              /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/elfun/sqrt.m"/* pName */
};

static emlrtRTEInfo f_emlrtRTEI = { 13,/* lineNo */
  13,                                  /* colNo */
  "toLogicalCheck",                    /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/toLogicalCheck.m"/* pName */
};

/* Function Declarations */
static void b_solveNewtonStep(c_logInteriorPoint_conjgrad_rtS *SD, const
  emlrtStack *sp, const real_T v[120], const real_T const_invW[14400], const
  real_T const_c[120], const real_T const_b[120], real_T const_maxCGIter, real_T
  const_CGTol, const real_T const_GDiag[120], const real_T const_A[14400], const
  real_T d0[120], real_T preCondFlag, real_T d[120], real_T dNext[120], real_T
  fNext[120], real_T MTilde[120], real_T *numIter, real_T *resStar, real_T
  *applyPreCond);
static void muStarSolve(const emlrtStack *sp, const real_T d0[120], const real_T
  d1[120], real_T mu_f, real_T *muStar, real_T d[120]);
static void solveNewtonStep(c_logInteriorPoint_conjgrad_rtS *SD, const
  emlrtStack *sp, const real_T v[120], const real_T const_invW[14400], const
  real_T const_c[120], const real_T const_b[120], real_T const_maxCGIter, real_T
  const_CGTol, const real_T const_GDiag[120], const real_T const_A[14400],
  real_T preCondFlag, real_T d[120], real_T dNext[120], real_T fNext[120],
  real_T MTilde[120], real_T *numIter, real_T *resStar, real_T *applyPreCond);
static real_T solveNewtonStep_warmStart(const emlrtStack *sp, const real_T b[120],
  const real_T v[120], const real_T const_invW[14400], real_T const_maxCGIter,
  real_T const_CGTol, const real_T const_A[14400], real_T d0[120], const real_T
  MTilde[120], real_T applyPreCond);

/* Function Definitions */
static void b_solveNewtonStep(c_logInteriorPoint_conjgrad_rtS *SD, const
  emlrtStack *sp, const real_T v[120], const real_T const_invW[14400], const
  real_T const_c[120], const real_T const_b[120], real_T const_maxCGIter, real_T
  const_CGTol, const real_T const_GDiag[120], const real_T const_A[14400], const
  real_T d0[120], real_T preCondFlag, real_T d[120], real_T dNext[120], real_T
  fNext[120], real_T MTilde[120], real_T *numIter, real_T *resStar, real_T
  *applyPreCond)
{
  real_T y_tmp[14400];
  real_T D[120];
  real_T b[120];
  real_T b_b[120];
  real_T p[120];
  real_T r[120];
  real_T rPrev[120];
  real_T x[120];
  real_T z[120];
  real_T zPrev[120];
  real_T a;
  real_T absxk;
  real_T alpha;
  real_T bHat_i;
  real_T bSum;
  real_T b_d;
  real_T beta;
  real_T res;
  real_T scale;
  real_T t;
  int32_T tmp_data[120];
  int32_T tmp_size[1];
  int32_T i;
  int32_T k;
  boolean_T b_v[120];

  /*  W = const.W; */
  /*  invW = const.invW; */
  /*  Define the preconditioner MTilde */
  /*  Define the RHS vector b */
  for (k = 0; k < 120; k++) {
    b_d = muDoubleScalarExp(v[k]);
    b[k] = b_d;
    MTilde[k] = b_d * const_GDiag[k] * b_d + 1.0;
    for (i = 0; i < 120; i++) {
      y_tmp[i + 120 * k] = const_A[k + 120 * i];
    }
  }

  for (k = 0; k < 120; k++) {
    b_b[k] = muDoubleScalarExp(v[k]);
  }

  for (i = 0; i < 14400; i++) {
    SD->u1.f0.dv[i] = 1.4142135623730951 * y_tmp[i];
  }

  mtimes(SD->u1.f0.dv, b_b, rPrev);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (i = 0; i < 120; i++) {
    rPrev[i] -= const_c[i];
  }

  mtimes(const_invW, rPrev, b_b);
  mtimes(const_A, b_b, rPrev);
  for (k = 0; k < 120; k++) {
    b_b[k] = muDoubleScalarExp(v[k]);
  }

  mtimes(y_tmp, b_b, zPrev);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (i = 0; i < 120; i++) {
    zPrev[i] -= const_c[i];
  }

  mtimes(const_invW, zPrev, b_b);
  mtimes(const_A, b_b, fNext);
  for (i = 0; i < 120; i++) {
    fNext[i] = 1.0 - b[i] * (fNext[i] + const_b[i]);
  }

  /*  First, determine whether or not to apply the diagonal preconditioner. Use */
  /*  a criterion than at least 1/4 of the variables have dropped below */
  /*  vTresh... This means that many elements of exp(v) will be near zero */
  for (i = 0; i < 120; i++) {
    b_v[i] = (v[i] < -4.0);
  }

  eml_find(b_v, tmp_data, tmp_size);
  if (preCondFlag == 1.0) {
    *applyPreCond = ((int8_T)tmp_size[0] > 30);
  } else {
    *applyPreCond = 0.0;
  }

  /*  --------------- CONJUGATE GRADIENT --------------- */
  /*  Initialize and redefine the problem such that x0 = 0. */
  /*  correct at the end by d = x + d0 */
  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  for (k = 0; k < 120; k++) {
    D[k] = muDoubleScalarExp(v[k]);
  }

  for (i = 0; i < 120; i++) {
    b_b[i] = D[i] * d0[i];
  }

  b_mtimes(const_A, b_b, zPrev);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  mtimes(const_invW, zPrev, b_b);
  mtimes(const_A, b_b, zPrev);
  for (i = 0; i < 120; i++) {
    D[i] = d0[i] + D[i] * zPrev[i];
  }

  /*  Run the first iteration of CG and iniialize iteration variables */
  /*  initialize */
  *resStar = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 120; k++) {
    b_d = (1.0 - 0.70710678118654746 * b[k] * (rPrev[k] + const_b[k])) - D[k];
    b[k] = b_d;
    d[k] = 0.0;
    absxk = muDoubleScalarAbs(b_d);
    if (absxk > scale) {
      t = scale / absxk;
      *resStar = *resStar * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      *resStar += t * t;
    }
  }

  *resStar = scale * muDoubleScalarSqrt(*resStar);

  /*  initialize */
  /*  Calculate iteration constants */
  if (*applyPreCond != 0.0) {
    for (k = 0; k < 120; k++) {
      z[k] = b[k] / MTilde[k];
    }

    /*  preconditioner step */
  } else {
    memcpy(&z[0], &b[0], 120U * sizeof(real_T));
  }

  memcpy(&p[0], &z[0], 120U * sizeof(real_T));

  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  for (k = 0; k < 120; k++) {
    zPrev[k] = muDoubleScalarExp(v[k]);
  }

  for (i = 0; i < 120; i++) {
    b_b[i] = zPrev[i] * z[i];
  }

  b_mtimes(const_A, b_b, rPrev);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  mtimes(const_invW, rPrev, b_b);
  mtimes(const_A, b_b, rPrev);
  for (i = 0; i < 120; i++) {
    zPrev[i] = z[i] + zPrev[i] * rPrev[i];
  }

  t = 0.0;
  absxk = 0.0;

  /*  Define/initialize xNext and the associated b vectors */
  for (i = 0; i < 120; i++) {
    b_d = z[i];
    t += b[i] * b_d;
    absxk += b_d * zPrev[i];
    D[i] = fNext[i] - D[i];
  }

  alpha = t / absxk;
  t = 0.0;
  for (i = 0; i < 120; i++) {
    t += b[i] * D[i];
  }

  bHat_i = t / (*resStar * *resStar);

  /*  update */
  bSum = bHat_i;

  /*  initialize */
  a = alpha * bHat_i;

  /*  Update x and r */
  /*  r1 */
  res = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 120; k++) {
    b_d = b[k];
    D[k] -= bHat_i * b_d;
    dNext[k] = a * b_d;
    x[k] = alpha * z[k];
    rPrev[k] = b_d;
    b_d -= alpha * zPrev[k];
    r[k] = b_d;
    absxk = muDoubleScalarAbs(b_d);
    if (absxk > scale) {
      t = scale / absxk;
      res = res * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      res += t * t;
    }
  }

  res = scale * muDoubleScalarSqrt(res);
  if (res < *resStar) {
    /*  Store d as minimal residual solution */
    memcpy(&d[0], &x[0], 120U * sizeof(real_T));
    *resStar = res;
  }

  *numIter = 1.0;

  /*  Iterate...  */
  /*  bNorm = norm(b,2); */
  while ((*numIter <= const_maxCGIter) && (res > const_CGTol)) {
    memcpy(&zPrev[0], &z[0], 120U * sizeof(real_T));
    if (*applyPreCond != 0.0) {
      for (k = 0; k < 120; k++) {
        z[k] = r[k] / MTilde[k];
      }
    } else {
      memcpy(&z[0], &r[0], 120U * sizeof(real_T));
    }

    t = 0.0;
    absxk = 0.0;
    for (i = 0; i < 120; i++) {
      t += r[i] * z[i];
      absxk += rPrev[i] * zPrev[i];
    }

    beta = t / absxk;
    for (i = 0; i < 120; i++) {
      p[i] = z[i] + beta * p[i];
    }

    /*  Function to evaluate M(v)*x */
    /*  Returns zOut = M(v)*zIn */
    for (k = 0; k < 120; k++) {
      zPrev[k] = muDoubleScalarExp(v[k]);
    }

    for (i = 0; i < 120; i++) {
      b_b[i] = zPrev[i] * p[i];
    }

    b_mtimes(const_A, b_b, rPrev);

    /*  This is a placeholder function for when we eventually use Riccatti */
    /*  Returns zOut = invW*zIn */
    mtimes(const_invW, rPrev, b_b);
    mtimes(const_A, b_b, rPrev);
    for (i = 0; i < 120; i++) {
      zPrev[i] = p[i] + zPrev[i] * rPrev[i];
    }

    absxk = 0.0;
    for (i = 0; i < 120; i++) {
      absxk += p[i] * zPrev[i];
    }

    alpha = t / absxk;

    /*  Update xNext and bNext */
    t = 0.0;
    for (i = 0; i < 120; i++) {
      t += r[i] * D[i];
    }

    bHat_i = t / (res * res);
    bSum += bHat_i;
    a = alpha * bSum;

    /*  Update x */
    res = 0.0;
    scale = 3.3121686421112381E-170;
    for (k = 0; k < 120; k++) {
      b_d = r[k];
      D[k] -= bHat_i * b_d;
      absxk = beta * b[k] + b_d;
      b[k] = absxk;
      dNext[k] += a * absxk;
      x[k] += alpha * p[k];
      rPrev[k] = b_d;
      b_d -= alpha * zPrev[k];
      r[k] = b_d;
      absxk = muDoubleScalarAbs(b_d);
      if (absxk > scale) {
        t = scale / absxk;
        res = res * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        res += t * t;
      }
    }

    res = scale * muDoubleScalarSqrt(res);
    if (res < *resStar) {
      memcpy(&d[0], &x[0], 120U * sizeof(real_T));
      *resStar = res;
    }

    /*  i++ */
    (*numIter)++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  /*  Undo the change of variables */
  for (i = 0; i < 120; i++) {
    d[i] += d0[i];
  }

  if (*applyPreCond != 0.0) {
    for (i = 0; i < 120; i++) {
      dNext[i] = dNext[i] / MTilde[i] + d0[i];
    }

    /*  note the rescaling for the precond case */
  } else {
    for (i = 0; i < 120; i++) {
      dNext[i] += d0[i];
    }
  }
}

static void muStarSolve(const emlrtStack *sp, const real_T d0[120], const real_T
  d1[120], real_T mu_f, real_T *muStar, real_T d[120])
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

  st.site = &jc_emlrtRSI;
  if (*muStar < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  upper_bound_i = 1.0 / muDoubleScalarSqrt(*muStar);
  for (i = 0; i < 120; i++) {
    d[i] = d0[i] + upper_bound_i * d1[i];
  }
}

static void solveNewtonStep(c_logInteriorPoint_conjgrad_rtS *SD, const
  emlrtStack *sp, const real_T v[120], const real_T const_invW[14400], const
  real_T const_c[120], const real_T const_b[120], real_T const_maxCGIter, real_T
  const_CGTol, const real_T const_GDiag[120], const real_T const_A[14400],
  real_T preCondFlag, real_T d[120], real_T dNext[120], real_T fNext[120],
  real_T MTilde[120], real_T *numIter, real_T *resStar, real_T *applyPreCond)
{
  real_T y_tmp[14400];
  real_T D[120];
  real_T b[120];
  real_T b_b[120];
  real_T p[120];
  real_T r[120];
  real_T rPrev[120];
  real_T x[120];
  real_T z[120];
  real_T zPrev[120];
  real_T a;
  real_T absxk;
  real_T alpha;
  real_T bHat_i;
  real_T bSum;
  real_T b_d;
  real_T beta;
  real_T res;
  real_T scale;
  real_T t;
  int32_T tmp_data[120];
  int32_T tmp_size[1];
  int32_T i;
  int32_T k;
  boolean_T b_v[120];

  /*  W = const.W; */
  /*  invW = const.invW; */
  /*  Define the preconditioner MTilde */
  /*  Define the RHS vector b */
  for (k = 0; k < 120; k++) {
    b_d = muDoubleScalarExp(v[k]);
    b[k] = b_d;
    MTilde[k] = b_d * const_GDiag[k] * b_d + 1.0;
    for (i = 0; i < 120; i++) {
      y_tmp[i + 120 * k] = const_A[k + 120 * i];
    }
  }

  for (k = 0; k < 120; k++) {
    b_b[k] = muDoubleScalarExp(v[k]);
  }

  for (i = 0; i < 14400; i++) {
    SD->u1.f1.dv[i] = 1.4142135623730951 * y_tmp[i];
  }

  mtimes(SD->u1.f1.dv, b_b, rPrev);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (i = 0; i < 120; i++) {
    rPrev[i] -= const_c[i];
  }

  mtimes(const_invW, rPrev, b_b);
  mtimes(const_A, b_b, rPrev);
  for (k = 0; k < 120; k++) {
    b_b[k] = muDoubleScalarExp(v[k]);
  }

  mtimes(y_tmp, b_b, zPrev);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (i = 0; i < 120; i++) {
    zPrev[i] -= const_c[i];
  }

  mtimes(const_invW, zPrev, b_b);
  mtimes(const_A, b_b, fNext);
  for (i = 0; i < 120; i++) {
    fNext[i] = 1.0 - b[i] * (fNext[i] + const_b[i]);
  }

  /*  First, determine whether or not to apply the diagonal preconditioner. Use */
  /*  a criterion than at least 1/4 of the variables have dropped below */
  /*  vTresh... This means that many elements of exp(v) will be near zero */
  for (i = 0; i < 120; i++) {
    b_v[i] = (v[i] < -4.0);
  }

  eml_find(b_v, tmp_data, tmp_size);
  if (preCondFlag == 1.0) {
    *applyPreCond = ((int8_T)tmp_size[0] > 30);
  } else {
    *applyPreCond = 0.0;
  }

  /*  --------------- CONJUGATE GRADIENT --------------- */
  /*  Initialize and redefine the problem such that x0 = 0. */
  /*  correct at the end by d = x + d0 */
  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  for (k = 0; k < 120; k++) {
    D[k] = muDoubleScalarExp(v[k]);
  }

  for (i = 0; i < 120; i++) {
    b_b[i] = D[i] * 0.0;
  }

  b_mtimes(const_A, b_b, zPrev);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  mtimes(const_invW, zPrev, b_b);
  mtimes(const_A, b_b, zPrev);
  for (i = 0; i < 120; i++) {
    D[i] *= zPrev[i];
  }

  /*  Run the first iteration of CG and iniialize iteration variables */
  /*  initialize */
  *resStar = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 120; k++) {
    b_d = (1.0 - 0.70710678118654746 * b[k] * (rPrev[k] + const_b[k])) - D[k];
    b[k] = b_d;
    d[k] = 0.0;
    absxk = muDoubleScalarAbs(b_d);
    if (absxk > scale) {
      t = scale / absxk;
      *resStar = *resStar * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      *resStar += t * t;
    }
  }

  *resStar = scale * muDoubleScalarSqrt(*resStar);

  /*  initialize */
  /*  Calculate iteration constants */
  if (*applyPreCond != 0.0) {
    for (k = 0; k < 120; k++) {
      z[k] = b[k] / MTilde[k];
    }

    /*  preconditioner step */
  } else {
    memcpy(&z[0], &b[0], 120U * sizeof(real_T));
  }

  memcpy(&p[0], &z[0], 120U * sizeof(real_T));

  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  for (k = 0; k < 120; k++) {
    zPrev[k] = muDoubleScalarExp(v[k]);
  }

  for (i = 0; i < 120; i++) {
    b_b[i] = zPrev[i] * z[i];
  }

  b_mtimes(const_A, b_b, rPrev);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  mtimes(const_invW, rPrev, b_b);
  mtimes(const_A, b_b, rPrev);
  for (i = 0; i < 120; i++) {
    zPrev[i] = z[i] + zPrev[i] * rPrev[i];
  }

  t = 0.0;
  absxk = 0.0;

  /*  Define/initialize xNext and the associated b vectors */
  for (i = 0; i < 120; i++) {
    b_d = z[i];
    t += b[i] * b_d;
    absxk += b_d * zPrev[i];
    D[i] = fNext[i] - D[i];
  }

  alpha = t / absxk;
  t = 0.0;
  for (i = 0; i < 120; i++) {
    t += b[i] * D[i];
  }

  bHat_i = t / (*resStar * *resStar);

  /*  update */
  bSum = bHat_i;

  /*  initialize */
  a = alpha * bHat_i;

  /*  Update x and r */
  /*  r1 */
  res = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 120; k++) {
    b_d = b[k];
    D[k] -= bHat_i * b_d;
    dNext[k] = a * b_d;
    x[k] = alpha * z[k];
    rPrev[k] = b_d;
    b_d -= alpha * zPrev[k];
    r[k] = b_d;
    absxk = muDoubleScalarAbs(b_d);
    if (absxk > scale) {
      t = scale / absxk;
      res = res * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      res += t * t;
    }
  }

  res = scale * muDoubleScalarSqrt(res);
  if (res < *resStar) {
    /*  Store d as minimal residual solution */
    memcpy(&d[0], &x[0], 120U * sizeof(real_T));
    *resStar = res;
  }

  *numIter = 1.0;

  /*  Iterate...  */
  /*  bNorm = norm(b,2); */
  while ((*numIter <= const_maxCGIter) && (res > const_CGTol)) {
    memcpy(&zPrev[0], &z[0], 120U * sizeof(real_T));
    if (*applyPreCond != 0.0) {
      for (k = 0; k < 120; k++) {
        z[k] = r[k] / MTilde[k];
      }
    } else {
      memcpy(&z[0], &r[0], 120U * sizeof(real_T));
    }

    t = 0.0;
    absxk = 0.0;
    for (i = 0; i < 120; i++) {
      t += r[i] * z[i];
      absxk += rPrev[i] * zPrev[i];
    }

    beta = t / absxk;
    for (i = 0; i < 120; i++) {
      p[i] = z[i] + beta * p[i];
    }

    /*  Function to evaluate M(v)*x */
    /*  Returns zOut = M(v)*zIn */
    for (k = 0; k < 120; k++) {
      zPrev[k] = muDoubleScalarExp(v[k]);
    }

    for (i = 0; i < 120; i++) {
      b_b[i] = zPrev[i] * p[i];
    }

    b_mtimes(const_A, b_b, rPrev);

    /*  This is a placeholder function for when we eventually use Riccatti */
    /*  Returns zOut = invW*zIn */
    mtimes(const_invW, rPrev, b_b);
    mtimes(const_A, b_b, rPrev);
    for (i = 0; i < 120; i++) {
      zPrev[i] = p[i] + zPrev[i] * rPrev[i];
    }

    absxk = 0.0;
    for (i = 0; i < 120; i++) {
      absxk += p[i] * zPrev[i];
    }

    alpha = t / absxk;

    /*  Update xNext and bNext */
    t = 0.0;
    for (i = 0; i < 120; i++) {
      t += r[i] * D[i];
    }

    bHat_i = t / (res * res);
    bSum += bHat_i;
    a = alpha * bSum;

    /*  Update x */
    res = 0.0;
    scale = 3.3121686421112381E-170;
    for (k = 0; k < 120; k++) {
      b_d = r[k];
      D[k] -= bHat_i * b_d;
      absxk = beta * b[k] + b_d;
      b[k] = absxk;
      dNext[k] += a * absxk;
      x[k] += alpha * p[k];
      rPrev[k] = b_d;
      b_d -= alpha * zPrev[k];
      r[k] = b_d;
      absxk = muDoubleScalarAbs(b_d);
      if (absxk > scale) {
        t = scale / absxk;
        res = res * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        res += t * t;
      }
    }

    res = scale * muDoubleScalarSqrt(res);
    if (res < *resStar) {
      memcpy(&d[0], &x[0], 120U * sizeof(real_T));
      *resStar = res;
    }

    /*  i++ */
    (*numIter)++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  /*  Undo the change of variables */
  if (*applyPreCond != 0.0) {
    for (i = 0; i < 120; i++) {
      dNext[i] /= MTilde[i];
    }

    /*  note the rescaling for the precond case */
  }
}

static real_T solveNewtonStep_warmStart(const emlrtStack *sp, const real_T b[120],
  const real_T v[120], const real_T const_invW[14400], real_T const_maxCGIter,
  real_T const_CGTol, const real_T const_A[14400], real_T d0[120], const real_T
  MTilde[120], real_T applyPreCond)
{
  emlrtStack st;
  real_T D[120];
  real_T b_D[120];
  real_T p[120];
  real_T rPrev[120];
  real_T x[120];
  real_T z[120];
  real_T zPrev[120];
  real_T absxk;
  real_T alpha;
  real_T numIter;
  real_T res;
  real_T resStar;
  real_T scale;
  real_T t;
  int32_T k;
  st.prev = sp;
  st.tls = sp->tls;

  /*  --------------- CONJUGATE GRADIENT --------------- */
  /*  Initialize */
  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  for (k = 0; k < 120; k++) {
    D[k] = muDoubleScalarExp(v[k]);
  }

  for (k = 0; k < 120; k++) {
    rPrev[k] = D[k] * d0[k];
  }

  b_mtimes(const_A, rPrev, zPrev);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  mtimes(const_invW, zPrev, rPrev);
  mtimes(const_A, rPrev, zPrev);
  resStar = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 120; k++) {
    t = b[k] - (d0[k] + D[k] * zPrev[k]);
    D[k] = t;
    absxk = muDoubleScalarAbs(t);
    if (absxk > scale) {
      t = scale / absxk;
      resStar = resStar * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      resStar += t * t;
    }
  }

  resStar = scale * muDoubleScalarSqrt(resStar);

  /*  bNorm = norm(b,2); */
  st.site = &ec_emlrtRSI;
  if (muDoubleScalarIsNaN(applyPreCond)) {
    emlrtErrorWithMessageIdR2018a(&st, &f_emlrtRTEI, "MATLAB:nologicalnan",
      "MATLAB:nologicalnan", 0);
  }

  if (applyPreCond != 0.0) {
    for (k = 0; k < 120; k++) {
      z[k] = D[k] / MTilde[k];
    }

    /*  preconditioner step */
  } else {
    memcpy(&z[0], &D[0], 120U * sizeof(real_T));
  }

  memcpy(&p[0], &z[0], 120U * sizeof(real_T));

  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  for (k = 0; k < 120; k++) {
    b_D[k] = muDoubleScalarExp(v[k]);
  }

  for (k = 0; k < 120; k++) {
    rPrev[k] = b_D[k] * z[k];
  }

  b_mtimes(const_A, rPrev, zPrev);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  mtimes(const_invW, zPrev, rPrev);
  mtimes(const_A, rPrev, zPrev);
  for (k = 0; k < 120; k++) {
    b_D[k] = z[k] + b_D[k] * zPrev[k];
  }

  scale = 0.0;
  absxk = 0.0;
  for (k = 0; k < 120; k++) {
    t = z[k];
    scale += D[k] * t;
    absxk += t * b_D[k];
  }

  alpha = scale / absxk;
  res = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 120; k++) {
    x[k] = d0[k] + alpha * z[k];
    t = D[k];
    rPrev[k] = t;
    t -= alpha * b_D[k];
    D[k] = t;
    absxk = muDoubleScalarAbs(t);
    if (absxk > scale) {
      t = scale / absxk;
      res = res * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      res += t * t;
    }
  }

  res = scale * muDoubleScalarSqrt(res);
  if (res < resStar) {
    /*  Store d as minimal residual solution */
    memcpy(&d0[0], &x[0], 120U * sizeof(real_T));
    resStar = res;
  }

  numIter = 1.0;

  /*  Iterate...  */
  while ((numIter < const_maxCGIter) && (res > const_CGTol)) {
    memcpy(&zPrev[0], &z[0], 120U * sizeof(real_T));
    st.site = &gc_emlrtRSI;
    if (muDoubleScalarIsNaN(applyPreCond)) {
      emlrtErrorWithMessageIdR2018a(&st, &f_emlrtRTEI, "MATLAB:nologicalnan",
        "MATLAB:nologicalnan", 0);
    }

    if (applyPreCond != 0.0) {
      for (k = 0; k < 120; k++) {
        z[k] = D[k] / MTilde[k];
      }
    } else {
      memcpy(&z[0], &D[0], 120U * sizeof(real_T));
    }

    scale = 0.0;
    absxk = 0.0;
    for (k = 0; k < 120; k++) {
      scale += D[k] * z[k];
      absxk += rPrev[k] * zPrev[k];
    }

    absxk = scale / absxk;
    for (k = 0; k < 120; k++) {
      p[k] = z[k] + absxk * p[k];
    }

    /*  Function to evaluate M(v)*x */
    /*  Returns zOut = M(v)*zIn */
    for (k = 0; k < 120; k++) {
      b_D[k] = muDoubleScalarExp(v[k]);
    }

    for (k = 0; k < 120; k++) {
      rPrev[k] = b_D[k] * p[k];
    }

    b_mtimes(const_A, rPrev, zPrev);

    /*  This is a placeholder function for when we eventually use Riccatti */
    /*  Returns zOut = invW*zIn */
    mtimes(const_invW, zPrev, rPrev);
    mtimes(const_A, rPrev, zPrev);
    for (k = 0; k < 120; k++) {
      b_D[k] = p[k] + b_D[k] * zPrev[k];
    }

    absxk = 0.0;
    for (k = 0; k < 120; k++) {
      absxk += p[k] * b_D[k];
    }

    alpha = scale / absxk;
    res = 0.0;
    scale = 3.3121686421112381E-170;
    for (k = 0; k < 120; k++) {
      x[k] += alpha * p[k];
      t = D[k];
      rPrev[k] = t;
      t -= alpha * b_D[k];
      D[k] = t;
      absxk = muDoubleScalarAbs(t);
      if (absxk > scale) {
        t = scale / absxk;
        res = res * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        res += t * t;
      }
    }

    res = scale * muDoubleScalarSqrt(res);
    if (res < resStar) {
      memcpy(&d0[0], &x[0], 120U * sizeof(real_T));
      resStar = res;
    }

    /*  i++ */
    numIter++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  return numIter;
}

void logInteriorPoint_conjgrad_rt(c_logInteriorPoint_conjgrad_rtS *SD, const
  emlrtStack *sp, const real_T W[14400], const real_T c[120], real_T Aineq[14400],
  const real_T bineq[120], real_T mu_f, real_T mu_0, const real_T v0[120],
  real_T maxIter, real_T maxCGIter, real_T CGTol, real_T preCondFlag, real_T x
  [120], real_T *mu, real_T *execTime, real_T *numIter, real_T *totalCGIters)
{
  emlrtStack b_st;
  emlrtStack st;
  real_T D[120];
  real_T GDiag[120];
  real_T b[120];
  real_T b_Aineq[120];
  real_T b_b[120];
  real_T d[120];
  real_T d1[120];
  real_T p[120];
  real_T r[120];
  real_T unusedU3[120];
  real_T v[120];
  real_T xStar[120];
  real_T z[120];
  real_T zPrev[120];
  real_T CGIters;
  real_T alpha;
  real_T appliedPreCond;
  real_T b_d;
  real_T resStar;
  int32_T tmp_data[120];
  int32_T tmp_size[1];
  int32_T applyPreCond;
  int32_T i;
  int32_T k;
  int32_T startFlag;
  boolean_T b_v[120];
  boolean_T exitg1;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;

  /*  min 0.5*x'*W*x + c'*x   subject to:  A*x <= b */
  /*  Get size variables */
  /*  First, change variables to Ax + b >= 0... This is just for uniformity */
  /*  with quadprog's inputs. */
  for (i = 0; i < 14400; i++) {
    Aineq[i] = -Aineq[i];
  }

  /*  Pack */
  st.site = &emlrtRSI;
  inv(&st, W, SD->f2.invW);

  /*  Compute the diagonal values of G which we use for preconditioning in CG */
  for (k = 0; k < 120; k++) {
    st.site = &b_emlrtRSI;

    /*  This is a placeholder function for when we eventually use Riccatti */
    /*  Returns zOut = invW*zIn */
    for (i = 0; i < 120; i++) {
      b_Aineq[i] = Aineq[k + 120 * i];
    }

    mtimes(SD->f2.invW, b_Aineq, b);
    alpha = 0.0;
    for (i = 0; i < 120; i++) {
      alpha += Aineq[k + 120 * i] * b[i];
    }

    GDiag[k] = alpha;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  /*  Max amount of bisection iterations */
  *numIter = 0.0;

  /*  number of newton iterations performed */
  /*  --------------------- FIRST NEWTON ITERATION --------------------- */
  /*  If we provide a warm-start */
  st.site = &c_emlrtRSI;
  tic(&st);

  /*  On initialization pass we first check if we can find a muStar value */
  /*  within N iterations of biscection... If not, then we just initialize with */
  /*  mu0 */
  /*  First, we sample two points of (mu,d) and solve the linear system */
  /*  generated by d = d0 + k*d1, where k = 1/sqrt(mu) */
  /*  Sample point 1 */
  /*  Sample point 2 */
  /*  Run first Newton system, generate a warm-start d2Hat for the 2nd */
  st.site = &d_emlrtRSI;
  solveNewtonStep(SD, &st, v0, SD->f2.invW, c, bineq, maxCGIter, CGTol, GDiag,
                  Aineq, preCondFlag, xStar, D, zPrev, b, &CGIters, &alpha,
                  &appliedPreCond);

  /*  Run seconds Newton system starting at warm-start d2Hat */
  st.site = &e_emlrtRSI;
  alpha = solveNewtonStep_warmStart(&st, zPrev, v0, SD->f2.invW, maxCGIter,
    CGTol, Aineq, D, b, appliedPreCond);
  *totalCGIters = CGIters + alpha;

  /*  Obtain affine representation of d = d0 + k*d1 */
  /*  Solve for muStar */
  for (i = 0; i < 120; i++) {
    b_d = xStar[i];
    alpha = D[i];
    b_Aineq[i] = 3.4142135623730945 * b_d + -2.4142135623730945 * alpha;
    d1[i] = -3.4142135623730945 * (b_d - alpha);
  }

  st.site = &f_emlrtRSI;
  muStarSolve(&st, b_Aineq, d1, mu_f, mu, d);

  /*  If we found a muStar (or a feasible point), then use these values */
  if (!muDoubleScalarIsInf(*mu)) {
    startFlag = 1;

    /*  Update v - this updates regardless of whether or not we find a */
    /*  feasible mu... Either just iterates off mu0, or to muStar */
    st.site = &g_emlrtRSI;
    appliedPreCond = c_norm(d);
    appliedPreCond = muDoubleScalarMin(1.0, 1.0 / (appliedPreCond *
      appliedPreCond));
    for (k = 0; k < 120; k++) {
      v[k] = v0[k] + appliedPreCond * d[k];
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
    startFlag = 1;
    exitg1 = false;
    while ((!exitg1) && ((c_norm(d) > 1.0) || (startFlag == 1))) {
      if (startFlag == 1) {
        startFlag = 0;
        memset(&d[0], 0, 120U * sizeof(real_T));
      }

      /*  Solve for d */
      if (*numIter >= maxIter) {
        /*  We want at least aa feasible soluton */
        exitg1 = true;
      } else {
        /*          [d,cg,res,CGpcflag] = solveNewtonStep(mu,v,const,zeros(m,1)); */
        st.site = &h_emlrtRSI;

        /*  W = const.W; */
        /*  invW = const.invW; */
        /*  Define the preconditioner MTilde */
        for (k = 0; k < 120; k++) {
          b_d = muDoubleScalarExp(v[k]);
          b_b[k] = b_d;
          unusedU3[k] = b_d * GDiag[k] * b_d + 1.0;
        }

        /*  Define the RHS vector b */
        b_st.site = &kb_emlrtRSI;
        if (*mu < 0.0) {
          emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
            "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError",
            3, 4, 4, "sqrt");
        }

        appliedPreCond = 1.0 / muDoubleScalarSqrt(*mu);
        b_st.site = &kb_emlrtRSI;
        if (*mu < 0.0) {
          emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
            "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError",
            3, 4, 4, "sqrt");
        }

        alpha = muDoubleScalarSqrt(*mu);
        for (k = 0; k < 120; k++) {
          b[k] = muDoubleScalarExp(v[k]);
          for (i = 0; i < 120; i++) {
            SD->f2.y[i + 120 * k] = alpha * Aineq[k + 120 * i];
          }
        }

        mtimes(SD->f2.y, b, d1);

        /*  This is a placeholder function for when we eventually use Riccatti */
        /*  Returns zOut = invW*zIn */
        for (i = 0; i < 120; i++) {
          d1[i] -= c[i];
        }

        mtimes(SD->f2.invW, d1, b);
        mtimes(Aineq, b, d1);

        /*  First, determine whether or not to apply the diagonal preconditioner. Use */
        /*  a criterion than at least 1/4 of the variables have dropped below */
        /*  vTresh... This means that many elements of exp(v) will be near zero */
        for (i = 0; i < 120; i++) {
          b_v[i] = (v[i] < -4.0);
        }

        eml_find(b_v, tmp_data, tmp_size);
        if (preCondFlag == 1.0) {
          applyPreCond = ((int8_T)tmp_size[0] > 30);
        } else {
          applyPreCond = 0;
        }

        /*  --------------- CONJUGATE GRADIENT --------------- */
        /*  Initialize and redefine the problem such that x0 = 0. */
        /*  correct at the end by d = x + d0 */
        /*  Function to evaluate M(v)*x */
        /*  Returns zOut = M(v)*zIn */
        for (k = 0; k < 120; k++) {
          D[k] = muDoubleScalarExp(v[k]);
        }

        for (i = 0; i < 120; i++) {
          b_Aineq[i] = D[i] * d[i];
        }

        b_mtimes(Aineq, b_Aineq, zPrev);

        /*  This is a placeholder function for when we eventually use Riccatti */
        /*  Returns zOut = invW*zIn */
        mtimes(SD->f2.invW, zPrev, b);
        mtimes(Aineq, b, zPrev);

        /*  Run the first iteration of CG and iniialize iteration variables */
        for (k = 0; k < 120; k++) {
          b_b[k] = (1.0 - appliedPreCond * b_b[k] * (d1[k] + bineq[k])) - (d[k]
            + D[k] * zPrev[k]);
          xStar[k] = 0.0;
        }

        /*  initialize */
        resStar = b_norm(b_b);

        /*  initialize */
        /*  Calculate iteration constants */
        if (applyPreCond != 0) {
          for (k = 0; k < 120; k++) {
            z[k] = b_b[k] / unusedU3[k];
          }

          /*  preconditioner step */
        } else {
          memcpy(&z[0], &b_b[0], 120U * sizeof(real_T));
        }

        memcpy(&p[0], &z[0], 120U * sizeof(real_T));

        /*  Function to evaluate M(v)*x */
        /*  Returns zOut = M(v)*zIn */
        for (k = 0; k < 120; k++) {
          D[k] = muDoubleScalarExp(v[k]);
        }

        for (i = 0; i < 120; i++) {
          b_Aineq[i] = D[i] * z[i];
        }

        b_mtimes(Aineq, b_Aineq, d1);

        /*  This is a placeholder function for when we eventually use Riccatti */
        /*  Returns zOut = invW*zIn */
        mtimes(SD->f2.invW, d1, b);
        mtimes(Aineq, b, d1);
        for (i = 0; i < 120; i++) {
          D[i] = z[i] + D[i] * d1[i];
        }

        appliedPreCond = 0.0;
        alpha = 0.0;
        for (i = 0; i < 120; i++) {
          b_d = z[i];
          appliedPreCond += b_b[i] * b_d;
          alpha += b_d * D[i];
        }

        alpha = appliedPreCond / alpha;

        /*  Define/initialize xNext and the associated b vectors */
        /*  Update x and r */
        for (k = 0; k < 120; k++) {
          x[k] = alpha * z[k];
          r[k] = b_b[k] - alpha * D[k];
        }

        /*  r1 */
        alpha = b_norm(r);
        if (alpha < resStar) {
          /*  Store d as minimal residual solution */
          memcpy(&xStar[0], &x[0], 120U * sizeof(real_T));
          resStar = alpha;
        }

        CGIters = 1.0;

        /*  Iterate...  */
        /*  bNorm = norm(b,2); */
        while ((CGIters <= maxCGIter) && (alpha > CGTol)) {
          memcpy(&zPrev[0], &z[0], 120U * sizeof(real_T));
          if (applyPreCond != 0) {
            for (k = 0; k < 120; k++) {
              z[k] = r[k] / unusedU3[k];
            }
          } else {
            memcpy(&z[0], &r[0], 120U * sizeof(real_T));
          }

          alpha = 0.0;
          appliedPreCond = 0.0;
          for (i = 0; i < 120; i++) {
            alpha += r[i] * z[i];
            appliedPreCond += b_b[i] * zPrev[i];
          }

          appliedPreCond = alpha / appliedPreCond;
          for (i = 0; i < 120; i++) {
            p[i] = z[i] + appliedPreCond * p[i];
          }

          /*  Function to evaluate M(v)*x */
          /*  Returns zOut = M(v)*zIn */
          for (k = 0; k < 120; k++) {
            D[k] = muDoubleScalarExp(v[k]);
          }

          for (i = 0; i < 120; i++) {
            b_Aineq[i] = D[i] * p[i];
          }

          b_mtimes(Aineq, b_Aineq, d1);

          /*  This is a placeholder function for when we eventually use Riccatti */
          /*  Returns zOut = invW*zIn */
          mtimes(SD->f2.invW, d1, b);
          mtimes(Aineq, b, d1);
          for (i = 0; i < 120; i++) {
            D[i] = p[i] + D[i] * d1[i];
          }

          appliedPreCond = 0.0;
          for (i = 0; i < 120; i++) {
            appliedPreCond += p[i] * D[i];
          }

          alpha /= appliedPreCond;

          /*  Update xNext and bNext */
          /*  Update x */
          for (k = 0; k < 120; k++) {
            x[k] += alpha * p[k];
            b_d = r[k];
            b_b[k] = b_d;
            b_d -= alpha * D[k];
            r[k] = b_d;
          }

          alpha = b_norm(r);
          if (alpha < resStar) {
            memcpy(&xStar[0], &x[0], 120U * sizeof(real_T));
            resStar = alpha;
          }

          /*  i++ */
          CGIters++;
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b(&st);
          }
        }

        /*  Undo the change of variables */
        for (i = 0; i < 120; i++) {
          d[i] += xStar[i];
        }

        *totalCGIters += CGIters;

        /*  Update v */
        st.site = &i_emlrtRSI;
        appliedPreCond = c_norm(d);
        appliedPreCond = muDoubleScalarMin(1.0, 1.0 / (appliedPreCond *
          appliedPreCond));
        for (i = 0; i < 120; i++) {
          v[i] += appliedPreCond * d[i];
        }

        (*numIter)++;
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(sp);
        }
      }
    }
  }

  /*  --------------------- MAIN NEWTON ITERATION LOOP --------------------- */
  /*  Then, we finally run the main loop which selects muStar */
  memcpy(&zPrev[0], &v[0], 120U * sizeof(real_T));
  memcpy(&xStar[0], &d[0], 120U * sizeof(real_T));

  /*  initialize for use in WS the first step */
  alpha = 1.0;
  while (((*mu > mu_f) || (alpha > 1.0)) && (*numIter < maxIter)) {
    /*  First, we sample two points of (mu,d) and solve the linear system */
    /*  generated by d = d0 + k*d1, where k = 1/sqrt(mu) */
    /*  Run the first Newton system */
    st.site = &j_emlrtRSI;
    b_solveNewtonStep(SD, &st, v, SD->f2.invW, c, bineq, maxCGIter, CGTol, GDiag,
                      Aineq, xStar, preCondFlag, d1, D, zPrev, b, &CGIters,
                      &alpha, &appliedPreCond);
    memcpy(&xStar[0], &d1[0], 120U * sizeof(real_T));
    *totalCGIters += CGIters;

    /*  Run the second */
    st.site = &k_emlrtRSI;
    CGIters = solveNewtonStep_warmStart(&st, zPrev, v, SD->f2.invW, maxCGIter,
      CGTol, Aineq, D, b, appliedPreCond);
    *totalCGIters += CGIters;

    /*  Obtain affine representation of d = d0 + k*d1 */
    for (k = 0; k < 120; k++) {
      b_d = d1[k];
      alpha = D[k];
      zPrev[k] = 3.4142135623730945 * b_d + -2.4142135623730945 * alpha;
      b_d = -3.4142135623730945 * (b_d - alpha);
      d1[k] = b_d;
    }

    /*  Solve for muStar using bisection */
    st.site = &l_emlrtRSI;
    muStarSolve(&st, zPrev, d1, mu_f, &alpha, d);

    /*  Make sure that muStar is finite and catch it if not */
    if (muDoubleScalarIsInf(alpha)) {
      st.site = &m_emlrtRSI;
      if (*mu < 0.0) {
        emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
          "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3,
          4, 4, "sqrt");
      }

      appliedPreCond = 1.0 / muDoubleScalarSqrt(*mu);
      for (k = 0; k < 120; k++) {
        d[k] = zPrev[k] + appliedPreCond * d1[k];
      }

      /*  mu is unchanged, change d */
    } else {
      *mu = alpha;

      /*  mu is changed, we use the calculated d */
    }

    memcpy(&zPrev[0], &v[0], 120U * sizeof(real_T));

    /*  Update x, v, d */
    alpha = c_norm(d);
    st.site = &n_emlrtRSI;
    appliedPreCond = muDoubleScalarMin(1.0, 1.0 / (alpha * alpha));
    for (i = 0; i < 120; i++) {
      v[i] += appliedPreCond * d[i];
    }

    (*numIter)++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  st.site = &o_emlrtRSI;
  *execTime = toc(&st);
  st.site = &p_emlrtRSI;
  if (*mu < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  alpha = muDoubleScalarSqrt(*mu);
  st.site = &p_emlrtRSI;
  for (k = 0; k < 120; k++) {
    b_d = muDoubleScalarExp(zPrev[k]);
    zPrev[k] = b_d;
    for (i = 0; i < 120; i++) {
      SD->f2.y[i + 120 * k] = alpha * Aineq[k + 120 * i];
    }

    x[k] = b_d + b_d * d[k];
  }

  mtimes(SD->f2.y, x, d1);
  st.site = &p_emlrtRSI;
  for (i = 0; i < 120; i++) {
    d1[i] -= c[i];
  }

  mtimes(SD->f2.invW, d1, x);
}

/* End of code generation (logInteriorPoint_conjgrad_rt.c) */
