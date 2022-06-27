/*
 * logInteriorPoint_conjgrad_rt.c
 *
 * Code generation for function 'logInteriorPoint_conjgrad_rt'
 *
 */

/* Include files */
#include "logInteriorPoint_conjgrad_rt.h"
#include "find.h"
#include "indexShapeCheck.h"
#include "inv.h"
#include "logInteriorPoint_conjgrad_rt_data.h"
#include "logInteriorPoint_conjgrad_rt_emxutil.h"
#include "logInteriorPoint_conjgrad_rt_types.h"
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

static emlrtRSInfo b_emlrtRSI = { 36,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 54,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 57,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 66,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 76,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 90,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 106, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 107, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 112, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 128, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 131, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 140, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 144, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 150, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 156, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 161, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 162, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo ob_emlrtRSI = { 205,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo pb_emlrtRSI = { 228,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo qb_emlrtRSI = { 240,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo rb_emlrtRSI = { 264,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo sb_emlrtRSI = { 276,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo tb_emlrtRSI = { 301,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo ub_emlrtRSI = { 39, /* lineNo */
  "find",                              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/elmat/find.m"/* pathName */
};

static emlrtRSInfo yb_emlrtRSI = { 45, /* lineNo */
  "mpower",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/ops/mpower.m"/* pathName */
};

static emlrtRSInfo bc_emlrtRSI = { 325,/* lineNo */
  "solveNewtonStep_warmStart",         /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo cc_emlrtRSI = { 348,/* lineNo */
  "solveNewtonStep_warmStart",         /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo ec_emlrtRSI = { 413,/* lineNo */
  "muStarSolve",                       /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo fc_emlrtRSI = { 195,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  162,                                 /* lineNo */
  27,                                  /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_conjgrad_rt",      /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  162,                                 /* lineNo */
  25,                                  /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_conjgrad_rt",      /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  92,                                  /* lineNo */
  12,                                  /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_conjgrad_rt",      /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo emlrtDCI = { 89,    /* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_conjgrad_rt",      /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = { 89,  /* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_conjgrad_rt",      /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m",/* pName */
  4                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  152,                                 /* lineNo */
  5,                                   /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_conjgrad_rt",      /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  109,                                 /* lineNo */
  9,                                   /* colNo */
  "xError_vec",                        /* aName */
  "logInteriorPoint_conjgrad_rt",      /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m",/* pName */
  0                                    /* checkKind */
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

static emlrtRTEInfo g_emlrtRTEI = { 89,/* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_conjgrad_rt",      /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pName */
};

static emlrtRTEInfo h_emlrtRTEI = { 162,/* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_conjgrad_rt",      /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pName */
};

/* Function Declarations */
static void b_solveNewtonStep(const emlrtStack *sp, const real_T v[20], const
  real_T const_invW[100], const real_T const_c[10], const real_T const_A[200],
  const real_T const_b[20], real_T const_maxCGIter, const real_T const_GDiag[20],
  const real_T d0[20], real_T preCondFlag, real_T d[20], real_T dNext[20],
  real_T fNext[20], real_T MTilde[20], real_T *numIter, real_T *resStar, real_T *
  applyPreCond);
static void muStarSolve(const emlrtStack *sp, const real_T d0[20], const real_T
  d1[20], real_T mu_f, real_T *muStar, real_T d[20]);
static void solveNewtonStep(const emlrtStack *sp, const real_T v[20], const
  real_T const_invW[100], const real_T const_c[10], const real_T const_A[200],
  const real_T const_b[20], real_T const_maxCGIter, const real_T const_GDiag[20],
  real_T preCondFlag, real_T d[20], real_T dNext[20], real_T fNext[20], real_T
  MTilde[20], real_T *numIter, real_T *resStar, real_T *applyPreCond);
static void solveNewtonStep_warmStart(const emlrtStack *sp, const real_T b[20],
  const real_T v[20], const real_T const_invW[100], const real_T const_A[200],
  real_T const_maxCGIter, real_T d0[20], const real_T MTilde[20], real_T
  applyPreCond);

/* Function Definitions */
static void b_solveNewtonStep(const emlrtStack *sp, const real_T v[20], const
  real_T const_invW[100], const real_T const_c[10], const real_T const_A[200],
  const real_T const_b[20], real_T const_maxCGIter, const real_T const_GDiag[20],
  const real_T d0[20], real_T preCondFlag, real_T d[20], real_T dNext[20],
  real_T fNext[20], real_T MTilde[20], real_T *numIter, real_T *resStar, real_T *
  applyPreCond)
{
  emlrtStack b_st;
  emlrtStack st;
  real_T y_tmp[200];
  real_T D[20];
  real_T Md0[20];
  real_T b[20];
  real_T b_D[20];
  real_T p[20];
  real_T r[20];
  real_T w[20];
  real_T x[20];
  real_T z[20];
  real_T b_const_invW[10];
  real_T b_y_tmp[10];
  real_T a;
  real_T absxk;
  real_T alpha;
  real_T bHat_i;
  real_T bNorm;
  real_T bSum;
  real_T b_d;
  real_T beta;
  real_T d1;
  real_T maxIter;
  real_T res;
  real_T scale;
  real_T t;
  int32_T tmp_data[20];
  int32_T tmp_size[1];
  int32_T i;
  int32_T k;
  boolean_T b_v[20];
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;

  /*  W = const.W; */
  /*  invW = const.invW; */
  if ((const_maxCGIter > 20.0) || muDoubleScalarIsNaN(const_maxCGIter)) {
    maxIter = 20.0;
  } else {
    maxIter = const_maxCGIter;
  }

  /*  Define the preconditioner MTilde */
  /*  Define the RHS vector b */
  for (k = 0; k < 20; k++) {
    b_d = muDoubleScalarExp(v[k]);
    w[k] = b_d;
    MTilde[k] = b_d * const_GDiag[k] * b_d + 1.0;
    for (i = 0; i < 10; i++) {
      y_tmp[i + 10 * k] = const_A[k + 20 * i];
    }
  }

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 20; k++) {
      b_d += y_tmp[i + 10 * k] * w[k];
    }

    b_y_tmp[i] = b_d - const_c[i];
  }

  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_invW[i + 10 * k] * b_y_tmp[k];
    }

    b_const_invW[i] = b_d;
  }

  for (i = 0; i < 20; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_A[i + 20 * k] * b_const_invW[k];
    }

    fNext[i] = 1.0 - w[i] * (b_d + const_b[i]);
  }

  /*  First, determine whether or not to apply the diagonal preconditioner. Use */
  /*  a criterion than at least 1/4 of the variables have dropped below */
  /*  vTresh... This means that many elements of exp(v) will be near zero */
  st.site = &ob_emlrtRSI;
  for (i = 0; i < 20; i++) {
    b_v[i] = (v[i] < -4.0);
  }

  b_st.site = &ub_emlrtRSI;
  eml_find(&b_st, b_v, tmp_data, tmp_size);
  if (preCondFlag == 1.0) {
    *applyPreCond = ((int8_T)tmp_size[0] > 5);
  } else {
    *applyPreCond = 0.0;
  }

  /*  --------------- CONJUGATE GRADIENT --------------- */
  /*  Initialize and redefine the problem such that x0 = 0. */
  /*  correct at the end by d = x + d0 */
  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (k = 0; k < 20; k++) {
    b_d = muDoubleScalarExp(v[k]);
    D[k] = b_d;
    b_D[k] = b_d * d0[k];
  }

  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 20; k++) {
      b_d += const_A[k + 20 * i] * b_D[k];
    }

    b_y_tmp[i] = b_d;
  }

  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_invW[i + 10 * k] * b_y_tmp[k];
    }

    b_const_invW[i] = b_d;
  }

  for (i = 0; i < 20; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_A[i + 20 * k] * b_const_invW[k];
    }

    Md0[i] = d0[i] + D[i] * b_d;
  }

  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 20; k++) {
      b_d += 1.4142135623730951 * y_tmp[i + 10 * k] * w[k];
    }

    b_y_tmp[i] = b_d - const_c[i];
  }

  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_invW[i + 10 * k] * b_y_tmp[k];
    }

    b_const_invW[i] = b_d;
  }

  /*  Run the first iteration of CG and iniialize iteration variables */
  /*  initialize */
  *resStar = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 20; k++) {
    b_d = 0.0;
    for (i = 0; i < 10; i++) {
      b_d += const_A[k + 20 * i] * b_const_invW[i];
    }

    b_d = (1.0 - 0.70710678118654746 * w[k] * (b_d + const_b[k])) - Md0[k];
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
  st.site = &pb_emlrtRSI;
  if (*applyPreCond != 0.0) {
    for (k = 0; k < 20; k++) {
      z[k] = b[k] / MTilde[k];
    }

    /*  preconditioner step */
  } else {
    memcpy(&z[0], &b[0], 20U * sizeof(real_T));
  }

  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (k = 0; k < 20; k++) {
    b_d = z[k];
    p[k] = b_d;
    d1 = muDoubleScalarExp(v[k]);
    D[k] = d1;
    b_D[k] = d1 * b_d;
  }

  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 20; k++) {
      b_d += const_A[k + 20 * i] * b_D[k];
    }

    b_y_tmp[i] = b_d;
  }

  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_invW[i + 10 * k] * b_y_tmp[k];
    }

    b_const_invW[i] = b_d;
  }

  t = 0.0;
  absxk = 0.0;

  /*  Define/initialize xNext and the associated b vectors */
  for (i = 0; i < 20; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_A[i + 20 * k] * b_const_invW[k];
    }

    d1 = z[i];
    b_d = d1 + D[i] * b_d;
    w[i] = b_d;
    t += b[i] * d1;
    absxk += d1 * b_d;
    Md0[i] = fNext[i] - Md0[i];
  }

  alpha = t / absxk;
  st.site = &qb_emlrtRSI;
  t = 0.0;
  for (i = 0; i < 20; i++) {
    t += b[i] * Md0[i];
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
  for (k = 0; k < 20; k++) {
    b_d = b[k];
    Md0[k] -= bHat_i * b_d;
    dNext[k] = a * b_d;
    x[k] = alpha * z[k];
    D[k] = b_d;
    b_d -= alpha * w[k];
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
    memcpy(&d[0], &x[0], 20U * sizeof(real_T));
    *resStar = res;
  }

  *numIter = 1.0;

  /*  Iterate...  */
  bNorm = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 20; k++) {
    absxk = muDoubleScalarAbs(b[k]);
    if (absxk > scale) {
      t = scale / absxk;
      bNorm = bNorm * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      bNorm += t * t;
    }
  }

  bNorm = scale * muDoubleScalarSqrt(bNorm);
  while ((*numIter <= maxIter) && (res / bNorm > 1.0E-12)) {
    memcpy(&w[0], &z[0], 20U * sizeof(real_T));
    st.site = &rb_emlrtRSI;
    if (*applyPreCond != 0.0) {
      for (k = 0; k < 20; k++) {
        z[k] = r[k] / MTilde[k];
      }
    } else {
      memcpy(&z[0], &r[0], 20U * sizeof(real_T));
    }

    t = 0.0;
    absxk = 0.0;
    for (i = 0; i < 20; i++) {
      t += r[i] * z[i];
      absxk += D[i] * w[i];
    }

    beta = t / absxk;

    /*  Function to evaluate M(v)*x */
    /*  Returns zOut = M(v)*zIn */
    /*  This is a placeholder function for when we eventually use Riccatti */
    /*  Returns zOut = invW*zIn */
    for (k = 0; k < 20; k++) {
      b_d = z[k] + beta * p[k];
      p[k] = b_d;
      d1 = muDoubleScalarExp(v[k]);
      D[k] = d1;
      b_D[k] = d1 * b_d;
    }

    for (i = 0; i < 10; i++) {
      b_d = 0.0;
      for (k = 0; k < 20; k++) {
        b_d += const_A[k + 20 * i] * b_D[k];
      }

      b_y_tmp[i] = b_d;
    }

    for (i = 0; i < 10; i++) {
      b_d = 0.0;
      for (k = 0; k < 10; k++) {
        b_d += const_invW[i + 10 * k] * b_y_tmp[k];
      }

      b_const_invW[i] = b_d;
    }

    absxk = 0.0;
    for (i = 0; i < 20; i++) {
      b_d = 0.0;
      for (k = 0; k < 10; k++) {
        b_d += const_A[i + 20 * k] * b_const_invW[k];
      }

      d1 = p[i];
      b_d = d1 + D[i] * b_d;
      w[i] = b_d;
      absxk += d1 * b_d;
    }

    alpha = t / absxk;

    /*  Update xNext and bNext */
    st.site = &sb_emlrtRSI;
    t = 0.0;
    for (i = 0; i < 20; i++) {
      t += r[i] * Md0[i];
    }

    bHat_i = t / (res * res);
    bSum += bHat_i;
    a = alpha * bSum;

    /*  Update x */
    res = 0.0;
    scale = 3.3121686421112381E-170;
    for (k = 0; k < 20; k++) {
      b_d = r[k];
      Md0[k] -= bHat_i * b_d;
      d1 = beta * b[k] + b_d;
      b[k] = d1;
      dNext[k] += a * d1;
      x[k] += alpha * p[k];
      D[k] = b_d;
      b_d -= alpha * w[k];
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
      memcpy(&d[0], &x[0], 20U * sizeof(real_T));
      *resStar = res;
    }

    /*  i++ */
    (*numIter)++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  /*  Undo the change of variables */
  for (i = 0; i < 20; i++) {
    d[i] += d0[i];
  }

  st.site = &tb_emlrtRSI;
  if (*applyPreCond != 0.0) {
    for (i = 0; i < 20; i++) {
      dNext[i] = dNext[i] / MTilde[i] + d0[i];
    }

    /*  note the rescaling for the precond case */
  } else {
    for (i = 0; i < 20; i++) {
      dNext[i] += d0[i];
    }
  }
}

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

  st.site = &ec_emlrtRSI;
  if (*muStar < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  upper_bound_i = 1.0 / muDoubleScalarSqrt(*muStar);
  for (i = 0; i < 20; i++) {
    d[i] = d0[i] + upper_bound_i * d1[i];
  }
}

static void solveNewtonStep(const emlrtStack *sp, const real_T v[20], const
  real_T const_invW[100], const real_T const_c[10], const real_T const_A[200],
  const real_T const_b[20], real_T const_maxCGIter, const real_T const_GDiag[20],
  real_T preCondFlag, real_T d[20], real_T dNext[20], real_T fNext[20], real_T
  MTilde[20], real_T *numIter, real_T *resStar, real_T *applyPreCond)
{
  emlrtStack b_st;
  emlrtStack st;
  real_T y_tmp[200];
  real_T D[20];
  real_T Md0[20];
  real_T b[20];
  real_T p[20];
  real_T r[20];
  real_T w[20];
  real_T x[20];
  real_T z[20];
  real_T b_const_invW[10];
  real_T b_y_tmp[10];
  real_T a;
  real_T absxk;
  real_T alpha;
  real_T bHat_i;
  real_T bNorm;
  real_T bSum;
  real_T b_d;
  real_T beta;
  real_T d1;
  real_T maxIter;
  real_T res;
  real_T scale;
  real_T t;
  int32_T tmp_data[20];
  int32_T tmp_size[1];
  int32_T i;
  int32_T k;
  boolean_T b_v[20];
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;

  /*  W = const.W; */
  /*  invW = const.invW; */
  if ((const_maxCGIter > 20.0) || muDoubleScalarIsNaN(const_maxCGIter)) {
    maxIter = 20.0;
  } else {
    maxIter = const_maxCGIter;
  }

  /*  Define the preconditioner MTilde */
  /*  Define the RHS vector b */
  for (k = 0; k < 20; k++) {
    b_d = muDoubleScalarExp(v[k]);
    w[k] = b_d;
    MTilde[k] = b_d * const_GDiag[k] * b_d + 1.0;
    for (i = 0; i < 10; i++) {
      y_tmp[i + 10 * k] = const_A[k + 20 * i];
    }
  }

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 20; k++) {
      b_d += y_tmp[i + 10 * k] * w[k];
    }

    b_y_tmp[i] = b_d - const_c[i];
  }

  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_invW[i + 10 * k] * b_y_tmp[k];
    }

    b_const_invW[i] = b_d;
  }

  for (i = 0; i < 20; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_A[i + 20 * k] * b_const_invW[k];
    }

    fNext[i] = 1.0 - w[i] * (b_d + const_b[i]);
  }

  /*  First, determine whether or not to apply the diagonal preconditioner. Use */
  /*  a criterion than at least 1/4 of the variables have dropped below */
  /*  vTresh... This means that many elements of exp(v) will be near zero */
  st.site = &ob_emlrtRSI;
  for (i = 0; i < 20; i++) {
    b_v[i] = (v[i] < -4.0);
  }

  b_st.site = &ub_emlrtRSI;
  eml_find(&b_st, b_v, tmp_data, tmp_size);
  if (preCondFlag == 1.0) {
    *applyPreCond = ((int8_T)tmp_size[0] > 5);
  } else {
    *applyPreCond = 0.0;
  }

  /*  --------------- CONJUGATE GRADIENT --------------- */
  /*  Initialize and redefine the problem such that x0 = 0. */
  /*  correct at the end by d = x + d0 */
  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  for (k = 0; k < 20; k++) {
    D[k] = muDoubleScalarExp(v[k]);
  }

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 20; k++) {
      b_d += const_A[k + 20 * i] * (D[k] * 0.0);
    }

    b_y_tmp[i] = b_d;
  }

  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_invW[i + 10 * k] * b_y_tmp[k];
    }

    b_const_invW[i] = b_d;
  }

  for (i = 0; i < 20; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_A[i + 20 * k] * b_const_invW[k];
    }

    Md0[i] = D[i] * b_d;
  }

  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 20; k++) {
      b_d += 1.4142135623730951 * y_tmp[i + 10 * k] * w[k];
    }

    b_y_tmp[i] = b_d - const_c[i];
  }

  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_invW[i + 10 * k] * b_y_tmp[k];
    }

    b_const_invW[i] = b_d;
  }

  /*  Run the first iteration of CG and iniialize iteration variables */
  /*  initialize */
  *resStar = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 20; k++) {
    b_d = 0.0;
    for (i = 0; i < 10; i++) {
      b_d += const_A[k + 20 * i] * b_const_invW[i];
    }

    b_d = (1.0 - 0.70710678118654746 * w[k] * (b_d + const_b[k])) - Md0[k];
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
  st.site = &pb_emlrtRSI;
  if (*applyPreCond != 0.0) {
    for (k = 0; k < 20; k++) {
      z[k] = b[k] / MTilde[k];
    }

    /*  preconditioner step */
  } else {
    memcpy(&z[0], &b[0], 20U * sizeof(real_T));
  }

  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (k = 0; k < 20; k++) {
    b_d = z[k];
    p[k] = b_d;
    d1 = muDoubleScalarExp(v[k]);
    D[k] = d1;
    w[k] = d1 * b_d;
  }

  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 20; k++) {
      b_d += const_A[k + 20 * i] * w[k];
    }

    b_y_tmp[i] = b_d;
  }

  for (i = 0; i < 10; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_invW[i + 10 * k] * b_y_tmp[k];
    }

    b_const_invW[i] = b_d;
  }

  t = 0.0;
  absxk = 0.0;

  /*  Define/initialize xNext and the associated b vectors */
  for (i = 0; i < 20; i++) {
    b_d = 0.0;
    for (k = 0; k < 10; k++) {
      b_d += const_A[i + 20 * k] * b_const_invW[k];
    }

    d1 = z[i];
    b_d = d1 + D[i] * b_d;
    w[i] = b_d;
    t += b[i] * d1;
    absxk += d1 * b_d;
    Md0[i] = fNext[i] - Md0[i];
  }

  alpha = t / absxk;
  st.site = &qb_emlrtRSI;
  t = 0.0;
  for (i = 0; i < 20; i++) {
    t += b[i] * Md0[i];
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
  for (k = 0; k < 20; k++) {
    b_d = b[k];
    Md0[k] -= bHat_i * b_d;
    dNext[k] = a * b_d;
    x[k] = alpha * z[k];
    D[k] = b_d;
    b_d -= alpha * w[k];
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
    memcpy(&d[0], &x[0], 20U * sizeof(real_T));
    *resStar = res;
  }

  *numIter = 1.0;

  /*  Iterate...  */
  bNorm = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 20; k++) {
    absxk = muDoubleScalarAbs(b[k]);
    if (absxk > scale) {
      t = scale / absxk;
      bNorm = bNorm * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      bNorm += t * t;
    }
  }

  bNorm = scale * muDoubleScalarSqrt(bNorm);
  while ((*numIter <= maxIter) && (res / bNorm > 1.0E-12)) {
    memcpy(&w[0], &z[0], 20U * sizeof(real_T));
    st.site = &rb_emlrtRSI;
    if (*applyPreCond != 0.0) {
      for (k = 0; k < 20; k++) {
        z[k] = r[k] / MTilde[k];
      }
    } else {
      memcpy(&z[0], &r[0], 20U * sizeof(real_T));
    }

    t = 0.0;
    absxk = 0.0;
    for (i = 0; i < 20; i++) {
      t += r[i] * z[i];
      absxk += D[i] * w[i];
    }

    beta = t / absxk;

    /*  Function to evaluate M(v)*x */
    /*  Returns zOut = M(v)*zIn */
    /*  This is a placeholder function for when we eventually use Riccatti */
    /*  Returns zOut = invW*zIn */
    for (k = 0; k < 20; k++) {
      b_d = z[k] + beta * p[k];
      p[k] = b_d;
      d1 = muDoubleScalarExp(v[k]);
      D[k] = d1;
      w[k] = d1 * b_d;
    }

    for (i = 0; i < 10; i++) {
      b_d = 0.0;
      for (k = 0; k < 20; k++) {
        b_d += const_A[k + 20 * i] * w[k];
      }

      b_y_tmp[i] = b_d;
    }

    for (i = 0; i < 10; i++) {
      b_d = 0.0;
      for (k = 0; k < 10; k++) {
        b_d += const_invW[i + 10 * k] * b_y_tmp[k];
      }

      b_const_invW[i] = b_d;
    }

    absxk = 0.0;
    for (i = 0; i < 20; i++) {
      b_d = 0.0;
      for (k = 0; k < 10; k++) {
        b_d += const_A[i + 20 * k] * b_const_invW[k];
      }

      d1 = p[i];
      b_d = d1 + D[i] * b_d;
      w[i] = b_d;
      absxk += d1 * b_d;
    }

    alpha = t / absxk;

    /*  Update xNext and bNext */
    st.site = &sb_emlrtRSI;
    t = 0.0;
    for (i = 0; i < 20; i++) {
      t += r[i] * Md0[i];
    }

    bHat_i = t / (res * res);
    bSum += bHat_i;
    a = alpha * bSum;

    /*  Update x */
    res = 0.0;
    scale = 3.3121686421112381E-170;
    for (k = 0; k < 20; k++) {
      b_d = r[k];
      Md0[k] -= bHat_i * b_d;
      d1 = beta * b[k] + b_d;
      b[k] = d1;
      dNext[k] += a * d1;
      x[k] += alpha * p[k];
      D[k] = b_d;
      b_d -= alpha * w[k];
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
      memcpy(&d[0], &x[0], 20U * sizeof(real_T));
      *resStar = res;
    }

    /*  i++ */
    (*numIter)++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  /*  Undo the change of variables */
  st.site = &tb_emlrtRSI;
  if (*applyPreCond != 0.0) {
    for (i = 0; i < 20; i++) {
      dNext[i] /= MTilde[i];
    }

    /*  note the rescaling for the precond case */
  }
}

static void solveNewtonStep_warmStart(const emlrtStack *sp, const real_T b[20],
  const real_T v[20], const real_T const_invW[100], const real_T const_A[200],
  real_T const_maxCGIter, real_T d0[20], const real_T MTilde[20], real_T
  applyPreCond)
{
  emlrtStack st;
  real_T D[20];
  real_T p[20];
  real_T r[20];
  real_T rPrev[20];
  real_T w[20];
  real_T x[20];
  real_T z[20];
  real_T b_const_A[10];
  real_T b_const_invW[10];
  real_T absxk;
  real_T bNorm;
  real_T b_scale;
  real_T d;
  real_T maxIter;
  real_T res;
  real_T resStar;
  real_T scale;
  real_T t;
  int32_T i;
  int32_T k;
  int32_T numIter;
  st.prev = sp;
  st.tls = sp->tls;
  if ((const_maxCGIter > 20.0) || muDoubleScalarIsNaN(const_maxCGIter)) {
    maxIter = 20.0;
  } else {
    maxIter = const_maxCGIter;
  }

  /*  --------------- CONJUGATE GRADIENT --------------- */
  /*  Initialize */
  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (k = 0; k < 20; k++) {
    absxk = muDoubleScalarExp(v[k]);
    D[k] = absxk;
    w[k] = absxk * d0[k];
  }

  for (i = 0; i < 10; i++) {
    absxk = 0.0;
    for (k = 0; k < 20; k++) {
      absxk += const_A[k + 20 * i] * w[k];
    }

    b_const_A[i] = absxk;
  }

  for (i = 0; i < 10; i++) {
    absxk = 0.0;
    for (k = 0; k < 10; k++) {
      absxk += const_invW[i + 10 * k] * b_const_A[k];
    }

    b_const_invW[i] = absxk;
  }

  resStar = 0.0;
  scale = 3.3121686421112381E-170;
  bNorm = 0.0;
  b_scale = 3.3121686421112381E-170;
  for (k = 0; k < 20; k++) {
    absxk = 0.0;
    for (i = 0; i < 10; i++) {
      absxk += const_A[k + 20 * i] * b_const_invW[i];
    }

    d = b[k];
    absxk = d - (d0[k] + D[k] * absxk);
    r[k] = absxk;
    absxk = muDoubleScalarAbs(absxk);
    if (absxk > scale) {
      t = scale / absxk;
      resStar = resStar * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      resStar += t * t;
    }

    absxk = muDoubleScalarAbs(d);
    if (absxk > b_scale) {
      t = b_scale / absxk;
      bNorm = bNorm * t * t + 1.0;
      b_scale = absxk;
    } else {
      t = absxk / b_scale;
      bNorm += t * t;
    }
  }

  resStar = scale * muDoubleScalarSqrt(resStar);
  bNorm = b_scale * muDoubleScalarSqrt(bNorm);
  st.site = &bc_emlrtRSI;
  if (muDoubleScalarIsNaN(applyPreCond)) {
    emlrtErrorWithMessageIdR2018a(&st, &f_emlrtRTEI, "MATLAB:nologicalnan",
      "MATLAB:nologicalnan", 0);
  }

  if (applyPreCond != 0.0) {
    for (k = 0; k < 20; k++) {
      z[k] = r[k] / MTilde[k];
    }

    /*  preconditioner step */
  } else {
    memcpy(&z[0], &r[0], 20U * sizeof(real_T));
  }

  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (k = 0; k < 20; k++) {
    absxk = z[k];
    p[k] = absxk;
    d = muDoubleScalarExp(v[k]);
    D[k] = d;
    w[k] = d * absxk;
  }

  for (i = 0; i < 10; i++) {
    absxk = 0.0;
    for (k = 0; k < 20; k++) {
      absxk += const_A[k + 20 * i] * w[k];
    }

    b_const_A[i] = absxk;
  }

  for (i = 0; i < 10; i++) {
    absxk = 0.0;
    for (k = 0; k < 10; k++) {
      absxk += const_invW[i + 10 * k] * b_const_A[k];
    }

    b_const_invW[i] = absxk;
  }

  res = 0.0;
  b_scale = 0.0;
  for (i = 0; i < 20; i++) {
    absxk = 0.0;
    for (k = 0; k < 10; k++) {
      absxk += const_A[i + 20 * k] * b_const_invW[k];
    }

    d = z[i];
    absxk = d + D[i] * absxk;
    w[i] = absxk;
    res += r[i] * d;
    b_scale += d * absxk;
  }

  b_scale = res / b_scale;
  res = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 20; k++) {
    x[k] = d0[k] + b_scale * z[k];
    absxk = r[k];
    rPrev[k] = absxk;
    absxk -= b_scale * w[k];
    r[k] = absxk;
    absxk = muDoubleScalarAbs(absxk);
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
    memcpy(&d0[0], &x[0], 20U * sizeof(real_T));
    resStar = res;
  }

  numIter = 1;

  /*  Iterate...  */
  while ((numIter < maxIter) && (res / bNorm > 1.0E-12)) {
    memcpy(&D[0], &z[0], 20U * sizeof(real_T));
    st.site = &cc_emlrtRSI;
    if (muDoubleScalarIsNaN(applyPreCond)) {
      emlrtErrorWithMessageIdR2018a(&st, &f_emlrtRTEI, "MATLAB:nologicalnan",
        "MATLAB:nologicalnan", 0);
    }

    if (applyPreCond != 0.0) {
      for (k = 0; k < 20; k++) {
        z[k] = r[k] / MTilde[k];
      }
    } else {
      memcpy(&z[0], &r[0], 20U * sizeof(real_T));
    }

    res = 0.0;
    b_scale = 0.0;
    for (i = 0; i < 20; i++) {
      res += r[i] * z[i];
      b_scale += rPrev[i] * D[i];
    }

    b_scale = res / b_scale;

    /*  Function to evaluate M(v)*x */
    /*  Returns zOut = M(v)*zIn */
    /*  This is a placeholder function for when we eventually use Riccatti */
    /*  Returns zOut = invW*zIn */
    for (k = 0; k < 20; k++) {
      absxk = z[k] + b_scale * p[k];
      p[k] = absxk;
      d = muDoubleScalarExp(v[k]);
      D[k] = d;
      w[k] = d * absxk;
    }

    for (i = 0; i < 10; i++) {
      absxk = 0.0;
      for (k = 0; k < 20; k++) {
        absxk += const_A[k + 20 * i] * w[k];
      }

      b_const_A[i] = absxk;
    }

    for (i = 0; i < 10; i++) {
      absxk = 0.0;
      for (k = 0; k < 10; k++) {
        absxk += const_invW[i + 10 * k] * b_const_A[k];
      }

      b_const_invW[i] = absxk;
    }

    b_scale = 0.0;
    for (i = 0; i < 20; i++) {
      absxk = 0.0;
      for (k = 0; k < 10; k++) {
        absxk += const_A[i + 20 * k] * b_const_invW[k];
      }

      d = p[i];
      absxk = d + D[i] * absxk;
      w[i] = absxk;
      b_scale += d * absxk;
    }

    b_scale = res / b_scale;
    res = 0.0;
    scale = 3.3121686421112381E-170;
    for (k = 0; k < 20; k++) {
      x[k] += b_scale * p[k];
      absxk = r[k];
      rPrev[k] = absxk;
      absxk -= b_scale * w[k];
      r[k] = absxk;
      absxk = muDoubleScalarAbs(absxk);
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
      memcpy(&d0[0], &x[0], 20U * sizeof(real_T));
      resStar = res;
    }

    /*  i++ */
    numIter++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }
}

void logInteriorPoint_conjgrad_rt(const emlrtStack *sp, const real_T W[100],
  const real_T c[10], real_T Aineq[200], const real_T bineq[20], real_T mu_f,
  real_T mu_0, const real_T v0[20], real_T maxIter, real_T maxCGIter, real_T
  preCondFlag, const real_T xStar[10], real_T xTol, real_T x[10],
  emxArray_real_T *xError_vec, real_T *execTime, real_T *numIter)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  real_T y_tmp[200];
  real_T invW[100];
  real_T GDiag[20];
  real_T MTilde[20];
  real_T b[20];
  real_T b_d[20];
  real_T b_x[20];
  real_T b_xStar[20];
  real_T d2Hat[20];
  real_T p[20];
  real_T r[20];
  real_T v[20];
  real_T w[20];
  real_T z[20];
  real_T zPrev[20];
  real_T b_invW[10];
  real_T y[10];
  real_T a;
  real_T absx;
  real_T alpha;
  real_T bNorm;
  real_T b_maxIter;
  real_T d;
  real_T d1;
  real_T mu;
  real_T resStar;
  real_T xError;
  int32_T tmp_data[20];
  int32_T iv[2];
  int32_T tmp_size[1];
  int32_T b_i;
  int32_T b_numIter;
  int32_T exitg1;
  int32_T i;
  int32_T i1;
  int32_T init;
  int32_T startFlag;
  boolean_T b_v[20];
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;

  /*  min 0.5*x'*W*x + c'*x   subject to:  A*x <= b */
  /*  Get size variables */
  /*  First, change variables to Ax + b >= 0... This is just for uniformity */
  /*  with quadprog's inputs. */
  for (i = 0; i < 200; i++) {
    Aineq[i] = -Aineq[i];
  }

  /*  Pack */
  st.site = &emlrtRSI;
  inv(&st, W, invW);

  /*  Compute the diagonal values of G which we use for preconditioning in CG */
  for (b_i = 0; b_i < 20; b_i++) {
    /*  This is a placeholder function for when we eventually use Riccatti */
    /*  Returns zOut = invW*zIn */
    absx = 0.0;
    for (i = 0; i < 10; i++) {
      d = 0.0;
      for (i1 = 0; i1 < 10; i1++) {
        d += invW[i + 10 * i1] * Aineq[b_i + 20 * i1];
      }

      absx += Aineq[b_i + 20 * i] * d;
    }

    GDiag[b_i] = absx;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

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
  /*  First, we sample two points of (mu,d) and solve the linear system */
  /*  generated by d = d0 + k*d1, where k = 1/sqrt(mu) */
  /*  Sample point 1 */
  /*  Sample point 2 */
  /*  Run first Newton system, generate a warm-start d2Hat for the 2nd */
  st.site = &c_emlrtRSI;
  solveNewtonStep(&st, v0, invW, c, Aineq, bineq, maxCGIter, GDiag, preCondFlag,
                  b_xStar, d2Hat, zPrev, MTilde, &absx, &alpha, &xError);

  /*  Run seconds Newton system starting at warm-start d2Hat */
  st.site = &d_emlrtRSI;
  solveNewtonStep_warmStart(&st, zPrev, v0, invW, Aineq, maxCGIter, d2Hat,
    MTilde, xError);

  /*  Obtain affine representation of d = d0 + k*d1 */
  /*  Solve for muStar */
  for (i = 0; i < 20; i++) {
    d = b_xStar[i];
    d1 = d2Hat[i];
    w[i] = 3.4142135623730945 * d + -2.4142135623730945 * d1;
    b[i] = -3.4142135623730945 * (d - d1);
  }

  st.site = &e_emlrtRSI;
  muStarSolve(&st, w, b, mu_f, &mu, b_d);

  /*  If we found a muStar (or a feasible point), then use these values */
  if (!muDoubleScalarIsInf(mu)) {
    startFlag = 1;

    /*  Update v - this updates regardless of whether or not we find a */
    /*  feasible mu... Either just iterates off mu0, or to muStar */
    st.site = &f_emlrtRSI;
    a = 0.0;
    for (b_i = 0; b_i < 20; b_i++) {
      absx = muDoubleScalarAbs(b_d[b_i]);
      if (muDoubleScalarIsNaN(absx) || (absx > a)) {
        a = absx;
      }
    }

    b_st.site = &yb_emlrtRSI;
    a = muDoubleScalarMin(1.0, 1.0 / (a * a));
    for (b_i = 0; b_i < 20; b_i++) {
      v[b_i] = v0[b_i] + a * b_d[b_i];
    }

    *numIter = 1.0;
  } else {
    /*  Otherwise, truly give  up and cold start */
    mu = mu_0;

    /*  under the update at the end */
    for (b_i = 0; b_i < 20; b_i++) {
      b_d[b_i] *= rtInf;
      v[b_i] = 0.0;
    }

    startFlag = 0;
  }

  /*  If we have a pair (mu,d) with d > 1 and we iterate the Newton */
  /*  algorithm until we obtain convergence. If we found a muStar from our */
  /*  initialization, this loop will be skipped. */
  if (!(maxIter >= 0.0)) {
    emlrtNonNegativeCheckR2012b(maxIter, &b_emlrtDCI, sp);
  }

  d = (int32_T)muDoubleScalarFloor(maxIter);
  if (maxIter != d) {
    emlrtIntegerCheckR2012b(maxIter, &emlrtDCI, sp);
  }

  i = xError_vec->size[0];
  xError_vec->size[0] = (int32_T)maxIter;
  emxEnsureCapacity_real_T(sp, xError_vec, i, &g_emlrtRTEI);
  if (maxIter != d) {
    emlrtIntegerCheckR2012b(maxIter, &emlrtDCI, sp);
  }

  init = (int32_T)maxIter;
  for (i = 0; i < init; i++) {
    xError_vec->data[i] = 0.0;
  }

  for (b_i = 0; b_i < 20; b_i++) {
    zPrev[b_i] = muDoubleScalarExp(v[b_i]);
  }

  st.site = &g_emlrtRSI;
  if (mu < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  alpha = muDoubleScalarSqrt(mu);
  for (i = 0; i < 20; i++) {
    for (i1 = 0; i1 < 10; i1++) {
      y_tmp[i1 + 10 * i] = Aineq[i + 20 * i1];
    }

    d = zPrev[i];
    w[i] = d + d * b_d[i];
  }

  for (i = 0; i < 10; i++) {
    d = 0.0;
    for (i1 = 0; i1 < 20; i1++) {
      d += alpha * y_tmp[i + 10 * i1] * w[i1];
    }

    y[i] = d - c[i];
  }

  for (i = 0; i < 10; i++) {
    d = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      d += invW[i + 10 * i1] * y[i1];
    }

    x[i] = d;
  }

  for (i = 0; i < 10; i++) {
    y[i] = x[i] - xStar[i];
  }

  xError = c_norm(y);
  if (1 > (int32_T)maxIter) {
    emlrtDynamicBoundsCheckR2012b(1, 1, (int32_T)maxIter, &c_emlrtBCI, sp);
  }

  xError_vec->data[0] = xError;
  if (startFlag < 1) {
    init = 1;
    do {
      exitg1 = 0;
      alpha = 0.0;
      for (b_i = 0; b_i < 20; b_i++) {
        absx = muDoubleScalarAbs(b_d[b_i]);
        if (muDoubleScalarIsNaN(absx) || (absx > alpha)) {
          alpha = absx;
        }
      }

      if ((alpha > 1.0) || (init == 1)) {
        if (init == 1) {
          init = 0;
          memset(&b_d[0], 0, 20U * sizeof(real_T));
        }

        /*  Solve for d */
        if (*numIter >= maxIter) {
          /*  We want at least aa feasible soluton */
          exitg1 = 1;
        } else {
          /*          [d,cg,res,CGpcflag] = solveNewtonStep(mu,v,const,zeros(m,1)); */
          st.site = &h_emlrtRSI;

          /*  W = const.W; */
          /*  invW = const.invW; */
          if ((maxCGIter > 20.0) || muDoubleScalarIsNaN(maxCGIter)) {
            b_maxIter = 20.0;
          } else {
            b_maxIter = maxCGIter;
          }

          /*  Define the preconditioner MTilde */
          for (b_i = 0; b_i < 20; b_i++) {
            d = muDoubleScalarExp(v[b_i]);
            zPrev[b_i] = d;
            MTilde[b_i] = d * GDiag[b_i] * d + 1.0;
          }

          /*  Define the RHS vector b */
          b_st.site = &fc_emlrtRSI;
          if (mu < 0.0) {
            emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
              "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError",
              3, 4, 4, "sqrt");
          }

          a = 1.0 / muDoubleScalarSqrt(mu);
          b_st.site = &fc_emlrtRSI;
          if (mu < 0.0) {
            emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
              "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError",
              3, 4, 4, "sqrt");
          }

          alpha = muDoubleScalarSqrt(mu);

          /*  This is a placeholder function for when we eventually use Riccatti */
          /*  Returns zOut = invW*zIn */
          /*  First, determine whether or not to apply the diagonal preconditioner. Use */
          /*  a criterion than at least 1/4 of the variables have dropped below */
          /*  vTresh... This means that many elements of exp(v) will be near zero */
          b_st.site = &ob_emlrtRSI;
          for (i = 0; i < 20; i++) {
            b_v[i] = (v[i] < -4.0);
          }

          c_st.site = &ub_emlrtRSI;
          eml_find(&c_st, b_v, tmp_data, tmp_size);
          if (preCondFlag == 1.0) {
            startFlag = ((int8_T)tmp_size[0] > 5);
          } else {
            startFlag = 0;
          }

          /*  --------------- CONJUGATE GRADIENT --------------- */
          /*  Initialize and redefine the problem such that x0 = 0. */
          /*  correct at the end by d = x + d0 */
          /*  Function to evaluate M(v)*x */
          /*  Returns zOut = M(v)*zIn */
          for (b_i = 0; b_i < 20; b_i++) {
            d2Hat[b_i] = muDoubleScalarExp(v[b_i]);
          }

          /*  This is a placeholder function for when we eventually use Riccatti */
          /*  Returns zOut = invW*zIn */
          for (i = 0; i < 10; i++) {
            d = 0.0;
            for (i1 = 0; i1 < 20; i1++) {
              d += alpha * Aineq[i1 + 20 * i] * zPrev[i1];
            }

            y[i] = d - c[i];
          }

          for (i = 0; i < 10; i++) {
            d = 0.0;
            for (i1 = 0; i1 < 10; i1++) {
              d += invW[i + 10 * i1] * y[i1];
            }

            b_invW[i] = d;
          }

          for (i = 0; i < 20; i++) {
            d = 0.0;
            for (i1 = 0; i1 < 10; i1++) {
              d += Aineq[i + 20 * i1] * b_invW[i1];
            }

            b[i] = d + bineq[i];
            w[i] = d2Hat[i] * b_d[i];
          }

          for (i = 0; i < 10; i++) {
            d = 0.0;
            for (i1 = 0; i1 < 20; i1++) {
              d += Aineq[i1 + 20 * i] * w[i1];
            }

            y[i] = d;
          }

          for (i = 0; i < 10; i++) {
            d = 0.0;
            for (i1 = 0; i1 < 10; i1++) {
              d += invW[i + 10 * i1] * y[i1];
            }

            b_invW[i] = d;
          }

          /*  Run the first iteration of CG and iniialize iteration variables */
          for (b_i = 0; b_i < 20; b_i++) {
            d = 0.0;
            for (i = 0; i < 10; i++) {
              d += Aineq[b_i + 20 * i] * b_invW[i];
            }

            b[b_i] = (1.0 - a * zPrev[b_i] * b[b_i]) - (b_d[b_i] + d2Hat[b_i] *
              d);
            b_xStar[b_i] = 0.0;
          }

          /*  initialize */
          resStar = b_norm(b);

          /*  initialize */
          /*  Calculate iteration constants */
          b_st.site = &pb_emlrtRSI;
          if (startFlag != 0) {
            for (b_i = 0; b_i < 20; b_i++) {
              z[b_i] = b[b_i] / MTilde[b_i];
            }

            /*  preconditioner step */
          } else {
            memcpy(&z[0], &b[0], 20U * sizeof(real_T));
          }

          /*  Function to evaluate M(v)*x */
          /*  Returns zOut = M(v)*zIn */
          /*  This is a placeholder function for when we eventually use Riccatti */
          /*  Returns zOut = invW*zIn */
          for (b_i = 0; b_i < 20; b_i++) {
            d = z[b_i];
            p[b_i] = d;
            d1 = muDoubleScalarExp(v[b_i]);
            d2Hat[b_i] = d1;
            w[b_i] = d1 * d;
          }

          for (i = 0; i < 10; i++) {
            d = 0.0;
            for (i1 = 0; i1 < 20; i1++) {
              d += Aineq[i1 + 20 * i] * w[i1];
            }

            y[i] = d;
          }

          for (i = 0; i < 10; i++) {
            d = 0.0;
            for (i1 = 0; i1 < 10; i1++) {
              d += invW[i + 10 * i1] * y[i1];
            }

            b_invW[i] = d;
          }

          xError = 0.0;
          absx = 0.0;
          for (i = 0; i < 20; i++) {
            d = 0.0;
            for (i1 = 0; i1 < 10; i1++) {
              d += Aineq[i + 20 * i1] * b_invW[i1];
            }

            d1 = z[i];
            d = d1 + d2Hat[i] * d;
            w[i] = d;
            xError += b[i] * d1;
            absx += d1 * d;
          }

          alpha = xError / absx;

          /*  Define/initialize xNext and the associated b vectors */
          /*  Update x and r */
          for (b_i = 0; b_i < 20; b_i++) {
            b_x[b_i] = alpha * z[b_i];
            r[b_i] = b[b_i] - alpha * w[b_i];
          }

          /*  r1 */
          absx = b_norm(r);
          if (absx < resStar) {
            /*  Store d as minimal residual solution */
            memcpy(&b_xStar[0], &b_x[0], 20U * sizeof(real_T));
            resStar = absx;
          }

          b_numIter = 1;

          /*  Iterate...  */
          bNorm = b_norm(b);
          while ((b_numIter <= b_maxIter) && (absx / bNorm > 1.0E-12)) {
            memcpy(&zPrev[0], &z[0], 20U * sizeof(real_T));
            b_st.site = &rb_emlrtRSI;
            if (startFlag != 0) {
              for (b_i = 0; b_i < 20; b_i++) {
                z[b_i] = r[b_i] / MTilde[b_i];
              }
            } else {
              memcpy(&z[0], &r[0], 20U * sizeof(real_T));
            }

            alpha = 0.0;
            xError = 0.0;
            for (i = 0; i < 20; i++) {
              alpha += r[i] * z[i];
              xError += b[i] * zPrev[i];
            }

            a = alpha / xError;

            /*  Function to evaluate M(v)*x */
            /*  Returns zOut = M(v)*zIn */
            /*  This is a placeholder function for when we eventually use Riccatti */
            /*  Returns zOut = invW*zIn */
            for (b_i = 0; b_i < 20; b_i++) {
              d = z[b_i] + a * p[b_i];
              p[b_i] = d;
              d1 = muDoubleScalarExp(v[b_i]);
              d2Hat[b_i] = d1;
              w[b_i] = d1 * d;
            }

            for (i = 0; i < 10; i++) {
              d = 0.0;
              for (i1 = 0; i1 < 20; i1++) {
                d += Aineq[i1 + 20 * i] * w[i1];
              }

              y[i] = d;
            }

            for (i = 0; i < 10; i++) {
              d = 0.0;
              for (i1 = 0; i1 < 10; i1++) {
                d += invW[i + 10 * i1] * y[i1];
              }

              b_invW[i] = d;
            }

            absx = 0.0;
            for (i = 0; i < 20; i++) {
              d = 0.0;
              for (i1 = 0; i1 < 10; i1++) {
                d += Aineq[i + 20 * i1] * b_invW[i1];
              }

              d1 = p[i];
              d = d1 + d2Hat[i] * d;
              w[i] = d;
              absx += d1 * d;
            }

            alpha /= absx;

            /*  Update xNext and bNext */
            /*  Update x */
            for (b_i = 0; b_i < 20; b_i++) {
              b_x[b_i] += alpha * p[b_i];
              d = r[b_i];
              b[b_i] = d;
              d -= alpha * w[b_i];
              r[b_i] = d;
            }

            absx = b_norm(r);
            if (absx < resStar) {
              memcpy(&b_xStar[0], &b_x[0], 20U * sizeof(real_T));
              resStar = absx;
            }

            /*  i++ */
            b_numIter++;
            if (*emlrtBreakCheckR2012bFlagVar != 0) {
              emlrtBreakCheckR2012b(&st);
            }
          }

          /*  Undo the change of variables */
          for (i = 0; i < 20; i++) {
            b_d[i] += b_xStar[i];
          }

          for (b_i = 0; b_i < 20; b_i++) {
            zPrev[b_i] = muDoubleScalarExp(v[b_i]);
          }

          st.site = &i_emlrtRSI;
          if (mu < 0.0) {
            emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
              "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError",
              3, 4, 4, "sqrt");
          }

          alpha = muDoubleScalarSqrt(mu);
          for (i = 0; i < 20; i++) {
            d = zPrev[i];
            w[i] = d + d * b_d[i];
          }

          for (i = 0; i < 10; i++) {
            d = 0.0;
            for (i1 = 0; i1 < 20; i1++) {
              d += alpha * y_tmp[i + 10 * i1] * w[i1];
            }

            y[i] = d - c[i];
          }

          for (i = 0; i < 10; i++) {
            d = 0.0;
            for (i1 = 0; i1 < 10; i1++) {
              d += invW[i + 10 * i1] * y[i1];
            }

            x[i] = d;
          }

          for (i = 0; i < 10; i++) {
            y[i] = x[i] - xStar[i];
          }

          xError = c_norm(y);
          if (((int32_T)(*numIter + 1.0) < 1) || ((int32_T)(*numIter + 1.0) >
               xError_vec->size[0])) {
            emlrtDynamicBoundsCheckR2012b((int32_T)(*numIter + 1.0), 1,
              xError_vec->size[0], &e_emlrtBCI, sp);
          }

          xError_vec->data[(int32_T)(*numIter + 1.0) - 1] = xError;

          /*  Update v */
          st.site = &j_emlrtRSI;
          a = 0.0;
          for (b_i = 0; b_i < 20; b_i++) {
            absx = muDoubleScalarAbs(b_d[b_i]);
            if (muDoubleScalarIsNaN(absx) || (absx > a)) {
              a = absx;
            }
          }

          b_st.site = &yb_emlrtRSI;
          a = muDoubleScalarMin(1.0, 1.0 / (a * a));
          for (i = 0; i < 20; i++) {
            v[i] += a * b_d[i];
          }

          (*numIter)++;
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b(sp);
          }
        }
      } else {
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  /*  --------------------- MAIN NEWTON ITERATION LOOP --------------------- */
  /*  Then, we finally run the main loop which selects muStar */
  memcpy(&b_xStar[0], &b_d[0], 20U * sizeof(real_T));

  /*  initialize for use in WS the first step */
  alpha = 1.0;
  while (((xError > xTol) || (alpha > 1.0)) && (*numIter < maxIter)) {
    /*  First, we sample two points of (mu,d) and solve the linear system */
    /*  generated by d = d0 + k*d1, where k = 1/sqrt(mu) */
    /*  Run the first Newton system */
    st.site = &k_emlrtRSI;
    b_solveNewtonStep(&st, v, invW, c, Aineq, bineq, maxCGIter, GDiag, b_xStar,
                      preCondFlag, w, d2Hat, zPrev, MTilde, &absx, &alpha,
                      &xError);
    memcpy(&b_xStar[0], &w[0], 20U * sizeof(real_T));

    /*  Run the second */
    st.site = &l_emlrtRSI;
    solveNewtonStep_warmStart(&st, zPrev, v, invW, Aineq, maxCGIter, d2Hat,
      MTilde, xError);

    /*  Obtain affine representation of d = d0 + k*d1 */
    for (b_i = 0; b_i < 20; b_i++) {
      d = w[b_i];
      d1 = d2Hat[b_i];
      zPrev[b_i] = 3.4142135623730945 * d + -2.4142135623730945 * d1;
      d = -3.4142135623730945 * (d - d1);
      w[b_i] = d;
    }

    /*  Solve for muStar using bisection */
    st.site = &m_emlrtRSI;
    muStarSolve(&st, zPrev, w, mu_f, &absx, b_d);

    /*  Make sure that muStar is finite and catch it if not */
    if (muDoubleScalarIsInf(absx)) {
      st.site = &n_emlrtRSI;
      if (mu < 0.0) {
        emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
          "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3,
          4, 4, "sqrt");
      }

      a = 1.0 / muDoubleScalarSqrt(mu);
      for (b_i = 0; b_i < 20; b_i++) {
        b_d[b_i] = zPrev[b_i] + a * w[b_i];
      }

      /*  mu is unchanged, change d */
    } else {
      mu = absx;

      /*  mu is changed, we use the calculated d */
    }

    /*  Get xError */
    for (b_i = 0; b_i < 20; b_i++) {
      zPrev[b_i] = muDoubleScalarExp(v[b_i]);
    }

    st.site = &o_emlrtRSI;
    if (mu < 0.0) {
      emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    alpha = muDoubleScalarSqrt(mu);
    for (i = 0; i < 20; i++) {
      d = zPrev[i];
      d += d * b_d[i];
      zPrev[i] = d;
    }

    for (i = 0; i < 10; i++) {
      d = 0.0;
      for (i1 = 0; i1 < 20; i1++) {
        d += alpha * y_tmp[i + 10 * i1] * zPrev[i1];
      }

      y[i] = d - c[i];
    }

    for (i = 0; i < 10; i++) {
      d = 0.0;
      for (i1 = 0; i1 < 10; i1++) {
        d += invW[i + 10 * i1] * y[i1];
      }

      x[i] = d;
    }

    for (i = 0; i < 10; i++) {
      y[i] = x[i] - xStar[i];
    }

    xError = c_norm(y);
    if (((int32_T)(*numIter + 1.0) < 1) || ((int32_T)(*numIter + 1.0) >
         xError_vec->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(*numIter + 1.0), 1,
        xError_vec->size[0], &d_emlrtBCI, sp);
    }

    xError_vec->data[(int32_T)(*numIter + 1.0) - 1] = xError;

    /*  Update x, v, d */
    alpha = 0.0;
    for (b_i = 0; b_i < 20; b_i++) {
      absx = muDoubleScalarAbs(b_d[b_i]);
      if (muDoubleScalarIsNaN(absx) || (absx > alpha)) {
        alpha = absx;
      }
    }

    st.site = &p_emlrtRSI;
    b_st.site = &yb_emlrtRSI;
    a = muDoubleScalarMin(1.0, 1.0 / (alpha * alpha));
    for (i = 0; i < 20; i++) {
      v[i] += a * b_d[i];
    }

    (*numIter)++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  st.site = &q_emlrtRSI;
  *execTime = toc(&st);
  if (1.0 > *numIter) {
    i = 0;
  } else {
    if (1 > xError_vec->size[0]) {
      emlrtDynamicBoundsCheckR2012b(1, 1, xError_vec->size[0], &b_emlrtBCI, sp);
    }

    if (((int32_T)*numIter < 1) || ((int32_T)*numIter > xError_vec->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)*numIter, 1, xError_vec->size[0],
        &emlrtBCI, sp);
    }

    i = (int32_T)*numIter;
  }

  iv[0] = 1;
  iv[1] = i;
  st.site = &r_emlrtRSI;
  indexShapeCheck(&st, xError_vec->size[0], iv);
  i1 = xError_vec->size[0];
  xError_vec->size[0] = i;
  emxEnsureCapacity_real_T(sp, xError_vec, i1, &h_emlrtRTEI);
}

/* End of code generation (logInteriorPoint_conjgrad_rt.c) */
