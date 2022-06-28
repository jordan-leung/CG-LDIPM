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

static emlrtRSInfo b_emlrtRSI = { 27,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 36,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 54,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 57,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 66,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 76,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 90,  /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 106, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 107, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 112, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 128, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 131, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 140, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 144, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 150, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 156, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 161, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo s_emlrtRSI = { 162, /* lineNo */
  "logInteriorPoint_conjgrad_rt",      /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo eb_emlrtRSI = { 169,/* lineNo */
  "invWTimes",                         /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo fb_emlrtRSI = { 79, /* lineNo */
  "eml_mtimes_helper",                 /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pathName */
};

static emlrtRSInfo nb_emlrtRSI = { 195,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo ob_emlrtRSI = { 197,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo pb_emlrtRSI = { 205,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo qb_emlrtRSI = { 219,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo rb_emlrtRSI = { 228,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo sb_emlrtRSI = { 234,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo tb_emlrtRSI = { 240,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo ub_emlrtRSI = { 264,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo vb_emlrtRSI = { 271,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo wb_emlrtRSI = { 276,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo xb_emlrtRSI = { 301,/* lineNo */
  "solveNewtonStep",                   /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo yb_emlrtRSI = { 39, /* lineNo */
  "find",                              /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/elmat/find.m"/* pathName */
};

static emlrtRSInfo ec_emlrtRSI = { 45, /* lineNo */
  "mpower",                            /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/ops/mpower.m"/* pathName */
};

static emlrtRSInfo hc_emlrtRSI = { 325,/* lineNo */
  "solveNewtonStep_warmStart",         /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo jc_emlrtRSI = { 348,/* lineNo */
  "solveNewtonStep_warmStart",         /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pathName */
};

static emlrtRSInfo mc_emlrtRSI = { 413,/* lineNo */
  "muStarSolve",                       /* fcnName */
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

static emlrtRTEInfo e_emlrtRTEI = { 13,/* lineNo */
  9,                                   /* colNo */
  "sqrt",                              /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/elfun/sqrt.m"/* pName */
};

static emlrtRTEInfo g_emlrtRTEI = { 13,/* lineNo */
  13,                                  /* colNo */
  "toLogicalCheck",                    /* fName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/toLogicalCheck.m"/* pName */
};

static emlrtRTEInfo h_emlrtRTEI = { 89,/* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_conjgrad_rt",      /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pName */
};

static emlrtRTEInfo i_emlrtRTEI = { 162,/* lineNo */
  1,                                   /* colNo */
  "logInteriorPoint_conjgrad_rt",      /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/logInteriorPoint_conjgrad_rt.m"/* pName */
};

/* Function Declarations */
static void b_solveNewtonStep(const emlrtStack *sp, const real_T v[50], const
  real_T const_invW[10000], const real_T const_c[100], const real_T const_A[5000],
  const real_T const_b[50], real_T const_maxCGIter, const real_T const_GDiag[50],
  const real_T d0[50], real_T preCondFlag, real_T d[50], real_T dNext[50],
  real_T fNext[50], real_T MTilde[50], real_T *numIter, real_T *resStar, real_T *
  applyPreCond);
static void muStarSolve(const emlrtStack *sp, const real_T d0[50], const real_T
  d1[50], real_T mu_f, real_T *muStar, real_T d[50]);
static void solveNewtonStep(const emlrtStack *sp, const real_T v[50], const
  real_T const_invW[10000], const real_T const_c[100], const real_T const_A[5000],
  const real_T const_b[50], real_T const_maxCGIter, const real_T const_GDiag[50],
  real_T preCondFlag, real_T d[50], real_T dNext[50], real_T fNext[50], real_T
  MTilde[50], real_T *numIter, real_T *resStar, real_T *applyPreCond);
static void solveNewtonStep_warmStart(const emlrtStack *sp, const real_T b[50],
  const real_T v[50], const real_T const_invW[10000], const real_T const_A[5000],
  real_T const_maxCGIter, real_T d0[50], const real_T MTilde[50], real_T
  applyPreCond);

/* Function Definitions */
static void b_solveNewtonStep(const emlrtStack *sp, const real_T v[50], const
  real_T const_invW[10000], const real_T const_c[100], const real_T const_A[5000],
  const real_T const_b[50], real_T const_maxCGIter, const real_T const_GDiag[50],
  const real_T d0[50], real_T preCondFlag, real_T d[50], real_T dNext[50],
  real_T fNext[50], real_T MTilde[50], real_T *numIter, real_T *resStar, real_T *
  applyPreCond)
{
  emlrtStack b_st;
  emlrtStack st;
  real_T dv[5000];
  real_T y_tmp[5000];
  real_T b_b[100];
  real_T y[100];
  real_T D[50];
  real_T b[50];
  real_T b_D[50];
  real_T p[50];
  real_T r[50];
  real_T rPrev[50];
  real_T x[50];
  real_T z[50];
  real_T zPrev[50];
  real_T alpha;
  real_T bHat_i;
  real_T bNorm;
  real_T bSum;
  real_T b_d;
  real_T b_r;
  real_T beta;
  real_T maxIter;
  real_T res;
  int32_T tmp_data[50];
  int32_T tmp_size[1];
  int32_T i;
  int32_T k;
  boolean_T b_v[50];
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;

  /*  W = const.W; */
  /*  invW = const.invW; */
  if ((const_maxCGIter > 50.0) || muDoubleScalarIsNaN(const_maxCGIter)) {
    maxIter = 50.0;
  } else {
    maxIter = const_maxCGIter;
  }

  /*  Define the preconditioner MTilde */
  /*  Define the RHS vector b */
  for (k = 0; k < 50; k++) {
    b_d = muDoubleScalarExp(v[k]);
    b[k] = b_d;
    MTilde[k] = b_d * const_GDiag[k] * b_d + 1.0;
    for (i = 0; i < 100; i++) {
      y_tmp[i + 100 * k] = const_A[k + 50 * i];
    }
  }

  st.site = &nb_emlrtRSI;
  for (k = 0; k < 50; k++) {
    zPrev[k] = muDoubleScalarExp(v[k]);
  }

  for (i = 0; i < 5000; i++) {
    dv[i] = 1.4142135623730951 * y_tmp[i];
  }

  b_mtimes(dv, zPrev, y);
  st.site = &nb_emlrtRSI;

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (i = 0; i < 100; i++) {
    y[i] -= const_c[i];
  }

  mtimes(const_invW, y, b_b);
  c_mtimes(const_A, b_b, rPrev);
  st.site = &ob_emlrtRSI;
  for (k = 0; k < 50; k++) {
    zPrev[k] = muDoubleScalarExp(v[k]);
  }

  b_mtimes(y_tmp, zPrev, y);
  st.site = &ob_emlrtRSI;

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (i = 0; i < 100; i++) {
    y[i] -= const_c[i];
  }

  mtimes(const_invW, y, b_b);
  c_mtimes(const_A, b_b, fNext);
  for (i = 0; i < 50; i++) {
    fNext[i] = 1.0 - b[i] * (fNext[i] + const_b[i]);
  }

  /*  First, determine whether or not to apply the diagonal preconditioner. Use */
  /*  a criterion than at least 1/4 of the variables have dropped below */
  /*  vTresh... This means that many elements of exp(v) will be near zero */
  st.site = &pb_emlrtRSI;
  for (i = 0; i < 50; i++) {
    b_v[i] = (v[i] < -4.0);
  }

  b_st.site = &yb_emlrtRSI;
  eml_find(&b_st, b_v, tmp_data, tmp_size);
  if (preCondFlag == 1.0) {
    *applyPreCond = ((int8_T)tmp_size[0] > 12.5);
  } else {
    *applyPreCond = 0.0;
  }

  /*  --------------- CONJUGATE GRADIENT --------------- */
  /*  Initialize and redefine the problem such that x0 = 0. */
  /*  correct at the end by d = x + d0 */
  st.site = &qb_emlrtRSI;

  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  for (k = 0; k < 50; k++) {
    D[k] = muDoubleScalarExp(v[k]);
  }

  for (i = 0; i < 50; i++) {
    b_D[i] = D[i] * d0[i];
  }

  d_mtimes(const_A, b_D, y);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  mtimes(const_invW, y, b_b);
  c_mtimes(const_A, b_b, zPrev);
  for (i = 0; i < 50; i++) {
    D[i] = d0[i] + D[i] * zPrev[i];
  }

  /*  Run the first iteration of CG and iniialize iteration variables */
  for (k = 0; k < 50; k++) {
    b[k] = (1.0 - 0.70710678118654746 * b[k] * (rPrev[k] + const_b[k])) - D[k];
    d[k] = 0.0;
  }

  /*  initialize */
  *resStar = b_norm(b);

  /*  initialize */
  /*  Calculate iteration constants */
  st.site = &rb_emlrtRSI;
  if (*applyPreCond != 0.0) {
    for (k = 0; k < 50; k++) {
      z[k] = b[k] / MTilde[k];
    }

    /*  preconditioner step */
  } else {
    memcpy(&z[0], &b[0], 50U * sizeof(real_T));
  }

  memcpy(&p[0], &z[0], 50U * sizeof(real_T));
  st.site = &sb_emlrtRSI;

  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  for (k = 0; k < 50; k++) {
    zPrev[k] = muDoubleScalarExp(v[k]);
  }

  for (i = 0; i < 50; i++) {
    b_D[i] = zPrev[i] * z[i];
  }

  d_mtimes(const_A, b_D, y);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  mtimes(const_invW, y, b_b);
  c_mtimes(const_A, b_b, rPrev);
  for (i = 0; i < 50; i++) {
    zPrev[i] = z[i] + zPrev[i] * rPrev[i];
  }

  res = 0.0;
  bHat_i = 0.0;

  /*  Define/initialize xNext and the associated b vectors */
  for (i = 0; i < 50; i++) {
    b_d = z[i];
    res += b[i] * b_d;
    bHat_i += b_d * zPrev[i];
    D[i] = fNext[i] - D[i];
  }

  alpha = res / bHat_i;
  st.site = &tb_emlrtRSI;
  res = 0.0;
  for (i = 0; i < 50; i++) {
    res += b[i] * D[i];
  }

  bHat_i = res / (*resStar * *resStar);

  /*  update */
  bSum = bHat_i;

  /*  initialize */
  res = alpha * bHat_i;

  /*  Update x and r */
  for (k = 0; k < 50; k++) {
    b_d = b[k];
    D[k] -= bHat_i * b_d;
    dNext[k] = res * b_d;
    x[k] = alpha * z[k];
    rPrev[k] = b_d;
    r[k] = b_d - alpha * zPrev[k];
  }

  /*  r1 */
  res = b_norm(r);
  if (res < *resStar) {
    /*  Store d as minimal residual solution */
    memcpy(&d[0], &x[0], 50U * sizeof(real_T));
    *resStar = res;
  }

  *numIter = 1.0;

  /*  Iterate...  */
  bNorm = b_norm(b);
  while ((*numIter <= maxIter) && (res / bNorm > 1.0E-12)) {
    memcpy(&zPrev[0], &z[0], 50U * sizeof(real_T));
    st.site = &ub_emlrtRSI;
    if (*applyPreCond != 0.0) {
      for (k = 0; k < 50; k++) {
        z[k] = r[k] / MTilde[k];
      }
    } else {
      memcpy(&z[0], &r[0], 50U * sizeof(real_T));
    }

    b_r = 0.0;
    bHat_i = 0.0;
    for (i = 0; i < 50; i++) {
      b_r += r[i] * z[i];
      bHat_i += rPrev[i] * zPrev[i];
    }

    beta = b_r / bHat_i;
    for (i = 0; i < 50; i++) {
      p[i] = z[i] + beta * p[i];
    }

    st.site = &vb_emlrtRSI;

    /*  Function to evaluate M(v)*x */
    /*  Returns zOut = M(v)*zIn */
    for (k = 0; k < 50; k++) {
      zPrev[k] = muDoubleScalarExp(v[k]);
    }

    for (i = 0; i < 50; i++) {
      b_D[i] = zPrev[i] * p[i];
    }

    d_mtimes(const_A, b_D, y);

    /*  This is a placeholder function for when we eventually use Riccatti */
    /*  Returns zOut = invW*zIn */
    mtimes(const_invW, y, b_b);
    c_mtimes(const_A, b_b, rPrev);
    for (i = 0; i < 50; i++) {
      zPrev[i] = p[i] + zPrev[i] * rPrev[i];
    }

    bHat_i = 0.0;
    for (i = 0; i < 50; i++) {
      bHat_i += p[i] * zPrev[i];
    }

    alpha = b_r / bHat_i;

    /*  Update xNext and bNext */
    st.site = &wb_emlrtRSI;
    b_r = 0.0;
    for (i = 0; i < 50; i++) {
      b_r += r[i] * D[i];
    }

    bHat_i = b_r / (res * res);
    bSum += bHat_i;
    res = alpha * bSum;

    /*  Update x */
    for (k = 0; k < 50; k++) {
      b_d = r[k];
      D[k] -= bHat_i * b_d;
      b_r = beta * b[k] + b_d;
      b[k] = b_r;
      dNext[k] += res * b_r;
      x[k] += alpha * p[k];
      rPrev[k] = b_d;
      b_d -= alpha * zPrev[k];
      r[k] = b_d;
    }

    res = b_norm(r);
    if (res < *resStar) {
      memcpy(&d[0], &x[0], 50U * sizeof(real_T));
      *resStar = res;
    }

    /*  i++ */
    (*numIter)++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  /*  Undo the change of variables */
  for (i = 0; i < 50; i++) {
    d[i] += d0[i];
  }

  st.site = &xb_emlrtRSI;
  if (*applyPreCond != 0.0) {
    for (i = 0; i < 50; i++) {
      dNext[i] = dNext[i] / MTilde[i] + d0[i];
    }

    /*  note the rescaling for the precond case */
  } else {
    for (i = 0; i < 50; i++) {
      dNext[i] += d0[i];
    }
  }
}

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

  st.site = &mc_emlrtRSI;
  if (*muStar < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  upper_bound_i = 1.0 / muDoubleScalarSqrt(*muStar);
  for (i = 0; i < 50; i++) {
    d[i] = d0[i] + upper_bound_i * d1[i];
  }
}

static void solveNewtonStep(const emlrtStack *sp, const real_T v[50], const
  real_T const_invW[10000], const real_T const_c[100], const real_T const_A[5000],
  const real_T const_b[50], real_T const_maxCGIter, const real_T const_GDiag[50],
  real_T preCondFlag, real_T d[50], real_T dNext[50], real_T fNext[50], real_T
  MTilde[50], real_T *numIter, real_T *resStar, real_T *applyPreCond)
{
  emlrtStack b_st;
  emlrtStack st;
  real_T dv[5000];
  real_T y_tmp[5000];
  real_T b_b[100];
  real_T y[100];
  real_T D[50];
  real_T b[50];
  real_T b_D[50];
  real_T p[50];
  real_T r[50];
  real_T rPrev[50];
  real_T x[50];
  real_T z[50];
  real_T zPrev[50];
  real_T alpha;
  real_T bHat_i;
  real_T bNorm;
  real_T bSum;
  real_T b_d;
  real_T b_r;
  real_T beta;
  real_T maxIter;
  real_T res;
  int32_T tmp_data[50];
  int32_T tmp_size[1];
  int32_T i;
  int32_T k;
  boolean_T b_v[50];
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;

  /*  W = const.W; */
  /*  invW = const.invW; */
  if ((const_maxCGIter > 50.0) || muDoubleScalarIsNaN(const_maxCGIter)) {
    maxIter = 50.0;
  } else {
    maxIter = const_maxCGIter;
  }

  /*  Define the preconditioner MTilde */
  /*  Define the RHS vector b */
  for (k = 0; k < 50; k++) {
    b_d = muDoubleScalarExp(v[k]);
    b[k] = b_d;
    MTilde[k] = b_d * const_GDiag[k] * b_d + 1.0;
    for (i = 0; i < 100; i++) {
      y_tmp[i + 100 * k] = const_A[k + 50 * i];
    }
  }

  st.site = &nb_emlrtRSI;
  for (k = 0; k < 50; k++) {
    zPrev[k] = muDoubleScalarExp(v[k]);
  }

  for (i = 0; i < 5000; i++) {
    dv[i] = 1.4142135623730951 * y_tmp[i];
  }

  b_mtimes(dv, zPrev, y);
  st.site = &nb_emlrtRSI;

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (i = 0; i < 100; i++) {
    y[i] -= const_c[i];
  }

  mtimes(const_invW, y, b_b);
  c_mtimes(const_A, b_b, rPrev);
  st.site = &ob_emlrtRSI;
  for (k = 0; k < 50; k++) {
    zPrev[k] = muDoubleScalarExp(v[k]);
  }

  b_mtimes(y_tmp, zPrev, y);
  st.site = &ob_emlrtRSI;

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  for (i = 0; i < 100; i++) {
    y[i] -= const_c[i];
  }

  mtimes(const_invW, y, b_b);
  c_mtimes(const_A, b_b, fNext);
  for (i = 0; i < 50; i++) {
    fNext[i] = 1.0 - b[i] * (fNext[i] + const_b[i]);
  }

  /*  First, determine whether or not to apply the diagonal preconditioner. Use */
  /*  a criterion than at least 1/4 of the variables have dropped below */
  /*  vTresh... This means that many elements of exp(v) will be near zero */
  st.site = &pb_emlrtRSI;
  for (i = 0; i < 50; i++) {
    b_v[i] = (v[i] < -4.0);
  }

  b_st.site = &yb_emlrtRSI;
  eml_find(&b_st, b_v, tmp_data, tmp_size);
  if (preCondFlag == 1.0) {
    *applyPreCond = ((int8_T)tmp_size[0] > 12.5);
  } else {
    *applyPreCond = 0.0;
  }

  /*  --------------- CONJUGATE GRADIENT --------------- */
  /*  Initialize and redefine the problem such that x0 = 0. */
  /*  correct at the end by d = x + d0 */
  st.site = &qb_emlrtRSI;

  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  for (k = 0; k < 50; k++) {
    D[k] = muDoubleScalarExp(v[k]);
  }

  for (i = 0; i < 50; i++) {
    b_D[i] = D[i] * 0.0;
  }

  d_mtimes(const_A, b_D, y);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  mtimes(const_invW, y, b_b);
  c_mtimes(const_A, b_b, zPrev);
  for (i = 0; i < 50; i++) {
    D[i] *= zPrev[i];
  }

  /*  Run the first iteration of CG and iniialize iteration variables */
  for (k = 0; k < 50; k++) {
    b[k] = (1.0 - 0.70710678118654746 * b[k] * (rPrev[k] + const_b[k])) - D[k];
    d[k] = 0.0;
  }

  /*  initialize */
  *resStar = b_norm(b);

  /*  initialize */
  /*  Calculate iteration constants */
  st.site = &rb_emlrtRSI;
  if (*applyPreCond != 0.0) {
    for (k = 0; k < 50; k++) {
      z[k] = b[k] / MTilde[k];
    }

    /*  preconditioner step */
  } else {
    memcpy(&z[0], &b[0], 50U * sizeof(real_T));
  }

  memcpy(&p[0], &z[0], 50U * sizeof(real_T));
  st.site = &sb_emlrtRSI;

  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  for (k = 0; k < 50; k++) {
    zPrev[k] = muDoubleScalarExp(v[k]);
  }

  for (i = 0; i < 50; i++) {
    b_D[i] = zPrev[i] * z[i];
  }

  d_mtimes(const_A, b_D, y);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  mtimes(const_invW, y, b_b);
  c_mtimes(const_A, b_b, rPrev);
  for (i = 0; i < 50; i++) {
    zPrev[i] = z[i] + zPrev[i] * rPrev[i];
  }

  res = 0.0;
  bHat_i = 0.0;

  /*  Define/initialize xNext and the associated b vectors */
  for (i = 0; i < 50; i++) {
    b_d = z[i];
    res += b[i] * b_d;
    bHat_i += b_d * zPrev[i];
    D[i] = fNext[i] - D[i];
  }

  alpha = res / bHat_i;
  st.site = &tb_emlrtRSI;
  res = 0.0;
  for (i = 0; i < 50; i++) {
    res += b[i] * D[i];
  }

  bHat_i = res / (*resStar * *resStar);

  /*  update */
  bSum = bHat_i;

  /*  initialize */
  res = alpha * bHat_i;

  /*  Update x and r */
  for (k = 0; k < 50; k++) {
    b_d = b[k];
    D[k] -= bHat_i * b_d;
    dNext[k] = res * b_d;
    x[k] = alpha * z[k];
    rPrev[k] = b_d;
    r[k] = b_d - alpha * zPrev[k];
  }

  /*  r1 */
  res = b_norm(r);
  if (res < *resStar) {
    /*  Store d as minimal residual solution */
    memcpy(&d[0], &x[0], 50U * sizeof(real_T));
    *resStar = res;
  }

  *numIter = 1.0;

  /*  Iterate...  */
  bNorm = b_norm(b);
  while ((*numIter <= maxIter) && (res / bNorm > 1.0E-12)) {
    memcpy(&zPrev[0], &z[0], 50U * sizeof(real_T));
    st.site = &ub_emlrtRSI;
    if (*applyPreCond != 0.0) {
      for (k = 0; k < 50; k++) {
        z[k] = r[k] / MTilde[k];
      }
    } else {
      memcpy(&z[0], &r[0], 50U * sizeof(real_T));
    }

    b_r = 0.0;
    bHat_i = 0.0;
    for (i = 0; i < 50; i++) {
      b_r += r[i] * z[i];
      bHat_i += rPrev[i] * zPrev[i];
    }

    beta = b_r / bHat_i;
    for (i = 0; i < 50; i++) {
      p[i] = z[i] + beta * p[i];
    }

    st.site = &vb_emlrtRSI;

    /*  Function to evaluate M(v)*x */
    /*  Returns zOut = M(v)*zIn */
    for (k = 0; k < 50; k++) {
      zPrev[k] = muDoubleScalarExp(v[k]);
    }

    for (i = 0; i < 50; i++) {
      b_D[i] = zPrev[i] * p[i];
    }

    d_mtimes(const_A, b_D, y);

    /*  This is a placeholder function for when we eventually use Riccatti */
    /*  Returns zOut = invW*zIn */
    mtimes(const_invW, y, b_b);
    c_mtimes(const_A, b_b, rPrev);
    for (i = 0; i < 50; i++) {
      zPrev[i] = p[i] + zPrev[i] * rPrev[i];
    }

    bHat_i = 0.0;
    for (i = 0; i < 50; i++) {
      bHat_i += p[i] * zPrev[i];
    }

    alpha = b_r / bHat_i;

    /*  Update xNext and bNext */
    st.site = &wb_emlrtRSI;
    b_r = 0.0;
    for (i = 0; i < 50; i++) {
      b_r += r[i] * D[i];
    }

    bHat_i = b_r / (res * res);
    bSum += bHat_i;
    res = alpha * bSum;

    /*  Update x */
    for (k = 0; k < 50; k++) {
      b_d = r[k];
      D[k] -= bHat_i * b_d;
      b_r = beta * b[k] + b_d;
      b[k] = b_r;
      dNext[k] += res * b_r;
      x[k] += alpha * p[k];
      rPrev[k] = b_d;
      b_d -= alpha * zPrev[k];
      r[k] = b_d;
    }

    res = b_norm(r);
    if (res < *resStar) {
      memcpy(&d[0], &x[0], 50U * sizeof(real_T));
      *resStar = res;
    }

    /*  i++ */
    (*numIter)++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  /*  Undo the change of variables */
  st.site = &xb_emlrtRSI;
  if (*applyPreCond != 0.0) {
    for (i = 0; i < 50; i++) {
      dNext[i] /= MTilde[i];
    }

    /*  note the rescaling for the precond case */
  }
}

static void solveNewtonStep_warmStart(const emlrtStack *sp, const real_T b[50],
  const real_T v[50], const real_T const_invW[10000], const real_T const_A[5000],
  real_T const_maxCGIter, real_T d0[50], const real_T MTilde[50], real_T
  applyPreCond)
{
  emlrtStack st;
  real_T b_b[100];
  real_T y[100];
  real_T D[50];
  real_T b_D[50];
  real_T p[50];
  real_T rPrev[50];
  real_T x[50];
  real_T z[50];
  real_T zPrev[50];
  real_T absxk;
  real_T alpha;
  real_T bNorm;
  real_T b_scale;
  real_T d;
  real_T maxIter;
  real_T resStar;
  real_T scale;
  real_T t;
  int32_T k;
  int32_T numIter;
  st.prev = sp;
  st.tls = sp->tls;
  if ((const_maxCGIter > 50.0) || muDoubleScalarIsNaN(const_maxCGIter)) {
    maxIter = 50.0;
  } else {
    maxIter = const_maxCGIter;
  }

  /*  --------------- CONJUGATE GRADIENT --------------- */
  /*  Initialize */
  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  for (k = 0; k < 50; k++) {
    D[k] = muDoubleScalarExp(v[k]);
  }

  for (k = 0; k < 50; k++) {
    zPrev[k] = D[k] * d0[k];
  }

  d_mtimes(const_A, zPrev, y);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  mtimes(const_invW, y, b_b);
  c_mtimes(const_A, b_b, zPrev);
  resStar = 0.0;
  scale = 3.3121686421112381E-170;
  bNorm = 0.0;
  b_scale = 3.3121686421112381E-170;
  for (k = 0; k < 50; k++) {
    d = b[k];
    alpha = d - (d0[k] + D[k] * zPrev[k]);
    D[k] = alpha;
    absxk = muDoubleScalarAbs(alpha);
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
  st.site = &hc_emlrtRSI;
  if (muDoubleScalarIsNaN(applyPreCond)) {
    emlrtErrorWithMessageIdR2018a(&st, &g_emlrtRTEI, "MATLAB:nologicalnan",
      "MATLAB:nologicalnan", 0);
  }

  if (applyPreCond != 0.0) {
    for (k = 0; k < 50; k++) {
      z[k] = D[k] / MTilde[k];
    }

    /*  preconditioner step */
  } else {
    memcpy(&z[0], &D[0], 50U * sizeof(real_T));
  }

  memcpy(&p[0], &z[0], 50U * sizeof(real_T));

  /*  Function to evaluate M(v)*x */
  /*  Returns zOut = M(v)*zIn */
  for (k = 0; k < 50; k++) {
    b_D[k] = muDoubleScalarExp(v[k]);
  }

  for (k = 0; k < 50; k++) {
    zPrev[k] = b_D[k] * z[k];
  }

  d_mtimes(const_A, zPrev, y);

  /*  This is a placeholder function for when we eventually use Riccatti */
  /*  Returns zOut = invW*zIn */
  mtimes(const_invW, y, b_b);
  c_mtimes(const_A, b_b, zPrev);
  for (k = 0; k < 50; k++) {
    b_D[k] = z[k] + b_D[k] * zPrev[k];
  }

  b_scale = 0.0;
  alpha = 0.0;
  for (k = 0; k < 50; k++) {
    d = z[k];
    b_scale += D[k] * d;
    alpha += d * b_D[k];
  }

  alpha = b_scale / alpha;
  b_scale = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 50; k++) {
    x[k] = d0[k] + alpha * z[k];
    d = D[k];
    rPrev[k] = d;
    d -= alpha * b_D[k];
    D[k] = d;
    absxk = muDoubleScalarAbs(d);
    if (absxk > scale) {
      t = scale / absxk;
      b_scale = b_scale * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      b_scale += t * t;
    }
  }

  b_scale = scale * muDoubleScalarSqrt(b_scale);
  if (b_scale < resStar) {
    /*  Store d as minimal residual solution */
    memcpy(&d0[0], &x[0], 50U * sizeof(real_T));
    resStar = b_scale;
  }

  numIter = 1;

  /*  Iterate...  */
  while ((numIter < maxIter) && (b_scale / bNorm > 1.0E-12)) {
    memcpy(&zPrev[0], &z[0], 50U * sizeof(real_T));
    st.site = &jc_emlrtRSI;
    if (muDoubleScalarIsNaN(applyPreCond)) {
      emlrtErrorWithMessageIdR2018a(&st, &g_emlrtRTEI, "MATLAB:nologicalnan",
        "MATLAB:nologicalnan", 0);
    }

    if (applyPreCond != 0.0) {
      for (k = 0; k < 50; k++) {
        z[k] = D[k] / MTilde[k];
      }
    } else {
      memcpy(&z[0], &D[0], 50U * sizeof(real_T));
    }

    b_scale = 0.0;
    alpha = 0.0;
    for (k = 0; k < 50; k++) {
      b_scale += D[k] * z[k];
      alpha += rPrev[k] * zPrev[k];
    }

    alpha = b_scale / alpha;
    for (k = 0; k < 50; k++) {
      p[k] = z[k] + alpha * p[k];
    }

    /*  Function to evaluate M(v)*x */
    /*  Returns zOut = M(v)*zIn */
    for (k = 0; k < 50; k++) {
      b_D[k] = muDoubleScalarExp(v[k]);
    }

    for (k = 0; k < 50; k++) {
      zPrev[k] = b_D[k] * p[k];
    }

    d_mtimes(const_A, zPrev, y);

    /*  This is a placeholder function for when we eventually use Riccatti */
    /*  Returns zOut = invW*zIn */
    mtimes(const_invW, y, b_b);
    c_mtimes(const_A, b_b, zPrev);
    for (k = 0; k < 50; k++) {
      b_D[k] = p[k] + b_D[k] * zPrev[k];
    }

    alpha = 0.0;
    for (k = 0; k < 50; k++) {
      alpha += p[k] * b_D[k];
    }

    alpha = b_scale / alpha;
    b_scale = 0.0;
    scale = 3.3121686421112381E-170;
    for (k = 0; k < 50; k++) {
      x[k] += alpha * p[k];
      d = D[k];
      rPrev[k] = d;
      d -= alpha * b_D[k];
      D[k] = d;
      absxk = muDoubleScalarAbs(d);
      if (absxk > scale) {
        t = scale / absxk;
        b_scale = b_scale * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        b_scale += t * t;
      }
    }

    b_scale = scale * muDoubleScalarSqrt(b_scale);
    if (b_scale < resStar) {
      memcpy(&d0[0], &x[0], 50U * sizeof(real_T));
      resStar = b_scale;
    }

    /*  i++ */
    numIter++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }
}

void logInteriorPoint_conjgrad_rt(c_logInteriorPoint_conjgrad_rtS *SD, const
  emlrtStack *sp, const real_T W[10000], const real_T c[100], real_T Aineq[5000],
  const real_T bineq[50], real_T mu_f, real_T mu_0, const real_T v0[50], real_T
  maxIter, real_T maxCGIter, real_T preCondFlag, const real_T xStar[100], real_T
  xTol, real_T x[100], emxArray_real_T *xError_vec, real_T *execTime, real_T
  *numIter)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  real_T y[5000];
  real_T y_tmp[5000];
  real_T b[100];
  real_T b_Aineq[100];
  real_T GDiag[50];
  real_T MTilde[50];
  real_T b_b[50];
  real_T b_x[50];
  real_T b_xStar[50];
  real_T d[50];
  real_T d1[50];
  real_T d2Hat[50];
  real_T p[50];
  real_T r[50];
  real_T v[50];
  real_T z[50];
  real_T zPrev[50];
  real_T a;
  real_T appliedPreCond;
  real_T bNorm;
  real_T b_maxIter;
  real_T dNorm;
  real_T mu;
  real_T resStar;
  real_T xError;
  int32_T tmp_data[50];
  int32_T iv[2];
  int32_T tmp_size[1];
  int32_T b_numIter;
  int32_T i;
  int32_T init;
  int32_T k;
  int32_T startFlag;
  boolean_T b_v[50];
  boolean_T exitg1;
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
  for (i = 0; i < 5000; i++) {
    Aineq[i] = -Aineq[i];
  }

  /*  Pack */
  st.site = &emlrtRSI;
  inv(&st, W, SD->f0.invW);

  /*  Compute the diagonal values of G which we use for preconditioning in CG */
  for (k = 0; k < 50; k++) {
    st.site = &b_emlrtRSI;

    /*  This is a placeholder function for when we eventually use Riccatti */
    /*  Returns zOut = invW*zIn */
    b_st.site = &eb_emlrtRSI;
    for (i = 0; i < 100; i++) {
      b_Aineq[i] = Aineq[k + 50 * i];
    }

    mtimes(SD->f0.invW, b_Aineq, b);
    xError = 0.0;
    for (i = 0; i < 100; i++) {
      xError += Aineq[k + 50 * i] * b[i];
    }

    GDiag[k] = xError;
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
  solveNewtonStep(&st, v0, SD->f0.invW, c, Aineq, bineq, maxCGIter, GDiag,
                  preCondFlag, b_xStar, d2Hat, zPrev, MTilde, &xError, &dNorm,
                  &appliedPreCond);

  /*  Run seconds Newton system starting at warm-start d2Hat */
  st.site = &e_emlrtRSI;
  solveNewtonStep_warmStart(&st, zPrev, v0, SD->f0.invW, Aineq, maxCGIter, d2Hat,
    MTilde, appliedPreCond);

  /*  Obtain affine representation of d = d0 + k*d1 */
  /*  Solve for muStar */
  for (i = 0; i < 50; i++) {
    appliedPreCond = b_xStar[i];
    xError = d2Hat[i];
    d1[i] = 3.4142135623730945 * appliedPreCond + -2.4142135623730945 * xError;
    zPrev[i] = -3.4142135623730945 * (appliedPreCond - xError);
  }

  st.site = &f_emlrtRSI;
  muStarSolve(&st, d1, zPrev, mu_f, &mu, d);

  /*  If we found a muStar (or a feasible point), then use these values */
  if (!muDoubleScalarIsInf(mu)) {
    startFlag = 1;

    /*  Update v - this updates regardless of whether or not we find a */
    /*  feasible mu... Either just iterates off mu0, or to muStar */
    st.site = &g_emlrtRSI;
    a = c_norm(d);
    b_st.site = &ec_emlrtRSI;
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

  /*  If we have a pair (mu,d) with d > 1 and we iterate the Newton */
  /*  algorithm until we obtain convergence. If we found a muStar from our */
  /*  initialization, this loop will be skipped. */
  if (!(maxIter >= 0.0)) {
    emlrtNonNegativeCheckR2012b(maxIter, &b_emlrtDCI, sp);
  }

  appliedPreCond = (int32_T)muDoubleScalarFloor(maxIter);
  if (maxIter != appliedPreCond) {
    emlrtIntegerCheckR2012b(maxIter, &emlrtDCI, sp);
  }

  i = xError_vec->size[0];
  xError_vec->size[0] = (int32_T)maxIter;
  emxEnsureCapacity_real_T(sp, xError_vec, i, &h_emlrtRTEI);
  if (maxIter != appliedPreCond) {
    emlrtIntegerCheckR2012b(maxIter, &emlrtDCI, sp);
  }

  init = (int32_T)maxIter;
  for (i = 0; i < init; i++) {
    xError_vec->data[i] = 0.0;
  }

  st.site = &h_emlrtRSI;
  if (mu < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  dNorm = muDoubleScalarSqrt(mu);
  for (i = 0; i < 50; i++) {
    for (init = 0; init < 100; init++) {
      y_tmp[init + 100 * i] = Aineq[i + 50 * init];
    }
  }

  st.site = &h_emlrtRSI;
  for (k = 0; k < 50; k++) {
    xError = muDoubleScalarExp(v[k]);
    b_x[k] = xError;
    zPrev[k] = xError;
  }

  for (i = 0; i < 5000; i++) {
    y[i] = dNorm * y_tmp[i];
  }

  for (i = 0; i < 50; i++) {
    d1[i] = b_x[i] + zPrev[i] * d[i];
  }

  b_st.site = &fb_emlrtRSI;
  b_mtimes(y, d1, b_Aineq);
  st.site = &h_emlrtRSI;
  for (i = 0; i < 100; i++) {
    b_Aineq[i] -= c[i];
  }

  b_st.site = &fb_emlrtRSI;
  mtimes(SD->f0.invW, b_Aineq, x);
  for (i = 0; i < 100; i++) {
    b_Aineq[i] = x[i] - xStar[i];
  }

  xError = d_norm(b_Aineq);
  if (1 > (int32_T)maxIter) {
    emlrtDynamicBoundsCheckR2012b(1, 1, (int32_T)maxIter, &c_emlrtBCI, sp);
  }

  xError_vec->data[0] = xError;
  if (startFlag < 1) {
    init = 1;
    exitg1 = false;
    while ((!exitg1) && ((c_norm(d) > 1.0) || (init == 1))) {
      if (init == 1) {
        init = 0;
        memset(&d[0], 0, 50U * sizeof(real_T));
      }

      /*  Solve for d */
      if (*numIter >= maxIter) {
        /*  We want at least aa feasible soluton */
        exitg1 = true;
      } else {
        /*          [d,cg,res,CGpcflag] = solveNewtonStep(mu,v,const,zeros(m,1)); */
        st.site = &i_emlrtRSI;

        /*  W = const.W; */
        /*  invW = const.invW; */
        if ((maxCGIter > 50.0) || muDoubleScalarIsNaN(maxCGIter)) {
          b_maxIter = 50.0;
        } else {
          b_maxIter = maxCGIter;
        }

        /*  Define the preconditioner MTilde */
        for (k = 0; k < 50; k++) {
          appliedPreCond = muDoubleScalarExp(v[k]);
          b_b[k] = appliedPreCond;
          MTilde[k] = appliedPreCond * GDiag[k] * appliedPreCond + 1.0;
        }

        /*  Define the RHS vector b */
        b_st.site = &nb_emlrtRSI;
        if (mu < 0.0) {
          emlrtErrorWithMessageIdR2018a(&b_st, &e_emlrtRTEI,
            "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError",
            3, 4, 4, "sqrt");
        }

        a = 1.0 / muDoubleScalarSqrt(mu);
        b_st.site = &nb_emlrtRSI;
        if (mu < 0.0) {
          emlrtErrorWithMessageIdR2018a(&b_st, &e_emlrtRTEI,
            "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError",
            3, 4, 4, "sqrt");
        }

        dNorm = muDoubleScalarSqrt(mu);
        b_st.site = &nb_emlrtRSI;
        for (k = 0; k < 50; k++) {
          zPrev[k] = muDoubleScalarExp(v[k]);
          for (i = 0; i < 100; i++) {
            y[i + 100 * k] = dNorm * Aineq[k + 50 * i];
          }
        }

        b_mtimes(y, zPrev, b_Aineq);
        b_st.site = &nb_emlrtRSI;

        /*  This is a placeholder function for when we eventually use Riccatti */
        /*  Returns zOut = invW*zIn */
        for (i = 0; i < 100; i++) {
          b_Aineq[i] -= c[i];
        }

        mtimes(SD->f0.invW, b_Aineq, b);
        c_mtimes(Aineq, b, d1);

        /*  First, determine whether or not to apply the diagonal preconditioner. Use */
        /*  a criterion than at least 1/4 of the variables have dropped below */
        /*  vTresh... This means that many elements of exp(v) will be near zero */
        b_st.site = &pb_emlrtRSI;
        for (i = 0; i < 50; i++) {
          b_v[i] = (v[i] < -4.0);
        }

        c_st.site = &yb_emlrtRSI;
        eml_find(&c_st, b_v, tmp_data, tmp_size);
        if (preCondFlag == 1.0) {
          startFlag = ((int8_T)tmp_size[0] > 12.5);
        } else {
          startFlag = 0;
        }

        /*  --------------- CONJUGATE GRADIENT --------------- */
        /*  Initialize and redefine the problem such that x0 = 0. */
        /*  correct at the end by d = x + d0 */
        b_st.site = &qb_emlrtRSI;

        /*  Function to evaluate M(v)*x */
        /*  Returns zOut = M(v)*zIn */
        for (k = 0; k < 50; k++) {
          d2Hat[k] = muDoubleScalarExp(v[k]);
        }

        for (i = 0; i < 50; i++) {
          zPrev[i] = d2Hat[i] * d[i];
        }

        d_mtimes(Aineq, zPrev, b_Aineq);

        /*  This is a placeholder function for when we eventually use Riccatti */
        /*  Returns zOut = invW*zIn */
        mtimes(SD->f0.invW, b_Aineq, b);
        c_mtimes(Aineq, b, zPrev);

        /*  Run the first iteration of CG and iniialize iteration variables */
        for (k = 0; k < 50; k++) {
          b_b[k] = (1.0 - a * b_b[k] * (d1[k] + bineq[k])) - (d[k] + d2Hat[k] *
            zPrev[k]);
          b_xStar[k] = 0.0;
        }

        /*  initialize */
        resStar = b_norm(b_b);

        /*  initialize */
        /*  Calculate iteration constants */
        b_st.site = &rb_emlrtRSI;
        if (startFlag != 0) {
          for (k = 0; k < 50; k++) {
            z[k] = b_b[k] / MTilde[k];
          }

          /*  preconditioner step */
        } else {
          memcpy(&z[0], &b_b[0], 50U * sizeof(real_T));
        }

        memcpy(&p[0], &z[0], 50U * sizeof(real_T));
        b_st.site = &sb_emlrtRSI;

        /*  Function to evaluate M(v)*x */
        /*  Returns zOut = M(v)*zIn */
        for (k = 0; k < 50; k++) {
          d2Hat[k] = muDoubleScalarExp(v[k]);
        }

        for (i = 0; i < 50; i++) {
          zPrev[i] = d2Hat[i] * z[i];
        }

        d_mtimes(Aineq, zPrev, b_Aineq);

        /*  This is a placeholder function for when we eventually use Riccatti */
        /*  Returns zOut = invW*zIn */
        mtimes(SD->f0.invW, b_Aineq, b);
        c_mtimes(Aineq, b, d1);
        for (i = 0; i < 50; i++) {
          d2Hat[i] = z[i] + d2Hat[i] * d1[i];
        }

        dNorm = 0.0;
        xError = 0.0;
        for (i = 0; i < 50; i++) {
          appliedPreCond = z[i];
          dNorm += b_b[i] * appliedPreCond;
          xError += appliedPreCond * d2Hat[i];
        }

        xError = dNorm / xError;

        /*  Define/initialize xNext and the associated b vectors */
        /*  Update x and r */
        for (k = 0; k < 50; k++) {
          b_x[k] = xError * z[k];
          r[k] = b_b[k] - xError * d2Hat[k];
        }

        /*  r1 */
        xError = b_norm(r);
        if (xError < resStar) {
          /*  Store d as minimal residual solution */
          memcpy(&b_xStar[0], &b_x[0], 50U * sizeof(real_T));
          resStar = xError;
        }

        b_numIter = 1;

        /*  Iterate...  */
        bNorm = b_norm(b_b);
        while ((b_numIter <= b_maxIter) && (xError / bNorm > 1.0E-12)) {
          memcpy(&zPrev[0], &z[0], 50U * sizeof(real_T));
          b_st.site = &ub_emlrtRSI;
          if (startFlag != 0) {
            for (k = 0; k < 50; k++) {
              z[k] = r[k] / MTilde[k];
            }
          } else {
            memcpy(&z[0], &r[0], 50U * sizeof(real_T));
          }

          appliedPreCond = 0.0;
          dNorm = 0.0;
          for (i = 0; i < 50; i++) {
            appliedPreCond += r[i] * z[i];
            dNorm += b_b[i] * zPrev[i];
          }

          a = appliedPreCond / dNorm;
          for (i = 0; i < 50; i++) {
            p[i] = z[i] + a * p[i];
          }

          b_st.site = &vb_emlrtRSI;

          /*  Function to evaluate M(v)*x */
          /*  Returns zOut = M(v)*zIn */
          for (k = 0; k < 50; k++) {
            d2Hat[k] = muDoubleScalarExp(v[k]);
          }

          for (i = 0; i < 50; i++) {
            zPrev[i] = d2Hat[i] * p[i];
          }

          d_mtimes(Aineq, zPrev, b_Aineq);

          /*  This is a placeholder function for when we eventually use Riccatti */
          /*  Returns zOut = invW*zIn */
          mtimes(SD->f0.invW, b_Aineq, b);
          c_mtimes(Aineq, b, d1);
          for (i = 0; i < 50; i++) {
            d2Hat[i] = p[i] + d2Hat[i] * d1[i];
          }

          xError = 0.0;
          for (i = 0; i < 50; i++) {
            xError += p[i] * d2Hat[i];
          }

          xError = appliedPreCond / xError;

          /*  Update xNext and bNext */
          /*  Update x */
          for (k = 0; k < 50; k++) {
            b_x[k] += xError * p[k];
            appliedPreCond = r[k];
            b_b[k] = appliedPreCond;
            appliedPreCond -= xError * d2Hat[k];
            r[k] = appliedPreCond;
          }

          xError = b_norm(r);
          if (xError < resStar) {
            memcpy(&b_xStar[0], &b_x[0], 50U * sizeof(real_T));
            resStar = xError;
          }

          /*  i++ */
          b_numIter++;
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b(&st);
          }
        }

        /*  Undo the change of variables */
        for (i = 0; i < 50; i++) {
          d[i] += b_xStar[i];
        }

        st.site = &j_emlrtRSI;
        if (mu < 0.0) {
          emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
            "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError",
            3, 4, 4, "sqrt");
        }

        dNorm = muDoubleScalarSqrt(mu);
        st.site = &j_emlrtRSI;
        for (k = 0; k < 50; k++) {
          xError = muDoubleScalarExp(v[k]);
          b_x[k] = xError;
          zPrev[k] = xError;
        }

        for (i = 0; i < 5000; i++) {
          y[i] = dNorm * y_tmp[i];
        }

        for (i = 0; i < 50; i++) {
          d1[i] = b_x[i] + zPrev[i] * d[i];
        }

        b_st.site = &fb_emlrtRSI;
        b_mtimes(y, d1, b_Aineq);
        st.site = &j_emlrtRSI;
        for (i = 0; i < 100; i++) {
          b_Aineq[i] -= c[i];
        }

        b_st.site = &fb_emlrtRSI;
        mtimes(SD->f0.invW, b_Aineq, x);
        for (i = 0; i < 100; i++) {
          b_Aineq[i] = x[i] - xStar[i];
        }

        xError = d_norm(b_Aineq);
        if (((int32_T)(*numIter + 1.0) < 1) || ((int32_T)(*numIter + 1.0) >
             xError_vec->size[0])) {
          emlrtDynamicBoundsCheckR2012b((int32_T)(*numIter + 1.0), 1,
            xError_vec->size[0], &e_emlrtBCI, sp);
        }

        xError_vec->data[(int32_T)(*numIter + 1.0) - 1] = xError;

        /*  Update v */
        st.site = &k_emlrtRSI;
        a = c_norm(d);
        b_st.site = &ec_emlrtRSI;
        a = muDoubleScalarMin(1.0, 1.0 / (a * a));
        for (i = 0; i < 50; i++) {
          v[i] += a * d[i];
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
  memcpy(&b_xStar[0], &d[0], 50U * sizeof(real_T));

  /*  initialize for use in WS the first step */
  dNorm = 1.0;
  while (((xError > xTol) || (dNorm > 1.0)) && (*numIter < maxIter)) {
    /*  First, we sample two points of (mu,d) and solve the linear system */
    /*  generated by d = d0 + k*d1, where k = 1/sqrt(mu) */
    /*  Run the first Newton system */
    st.site = &l_emlrtRSI;
    b_solveNewtonStep(&st, v, SD->f0.invW, c, Aineq, bineq, maxCGIter, GDiag,
                      b_xStar, preCondFlag, d1, d2Hat, zPrev, MTilde, &xError,
                      &dNorm, &appliedPreCond);
    memcpy(&b_xStar[0], &d1[0], 50U * sizeof(real_T));

    /*  Run the second */
    st.site = &m_emlrtRSI;
    solveNewtonStep_warmStart(&st, zPrev, v, SD->f0.invW, Aineq, maxCGIter,
      d2Hat, MTilde, appliedPreCond);

    /*  Obtain affine representation of d = d0 + k*d1 */
    for (k = 0; k < 50; k++) {
      appliedPreCond = d1[k];
      xError = d2Hat[k];
      zPrev[k] = 3.4142135623730945 * appliedPreCond + -2.4142135623730945 *
        xError;
      appliedPreCond = -3.4142135623730945 * (appliedPreCond - xError);
      d1[k] = appliedPreCond;
    }

    /*  Solve for muStar using bisection */
    st.site = &n_emlrtRSI;
    muStarSolve(&st, zPrev, d1, mu_f, &xError, d);

    /*  Make sure that muStar is finite and catch it if not */
    if (muDoubleScalarIsInf(xError)) {
      st.site = &o_emlrtRSI;
      if (mu < 0.0) {
        emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
          "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3,
          4, 4, "sqrt");
      }

      a = 1.0 / muDoubleScalarSqrt(mu);
      for (k = 0; k < 50; k++) {
        d[k] = zPrev[k] + a * d1[k];
      }

      /*  mu is unchanged, change d */
    } else {
      mu = xError;

      /*  mu is changed, we use the calculated d */
    }

    /*  Get xError */
    st.site = &p_emlrtRSI;
    if (mu < 0.0) {
      emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    dNorm = muDoubleScalarSqrt(mu);
    st.site = &p_emlrtRSI;
    for (k = 0; k < 50; k++) {
      xError = muDoubleScalarExp(v[k]);
      b_x[k] = xError;
      zPrev[k] = xError;
    }

    for (i = 0; i < 5000; i++) {
      y[i] = dNorm * y_tmp[i];
    }

    for (i = 0; i < 50; i++) {
      zPrev[i] = b_x[i] + zPrev[i] * d[i];
    }

    b_st.site = &fb_emlrtRSI;
    b_mtimes(y, zPrev, b_Aineq);
    st.site = &p_emlrtRSI;
    for (i = 0; i < 100; i++) {
      b_Aineq[i] -= c[i];
    }

    b_st.site = &fb_emlrtRSI;
    mtimes(SD->f0.invW, b_Aineq, x);
    for (i = 0; i < 100; i++) {
      b_Aineq[i] = x[i] - xStar[i];
    }

    xError = d_norm(b_Aineq);
    if (((int32_T)(*numIter + 1.0) < 1) || ((int32_T)(*numIter + 1.0) >
         xError_vec->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(*numIter + 1.0), 1,
        xError_vec->size[0], &d_emlrtBCI, sp);
    }

    xError_vec->data[(int32_T)(*numIter + 1.0) - 1] = xError;

    /*  Update x, v, d */
    dNorm = c_norm(d);
    st.site = &q_emlrtRSI;
    b_st.site = &ec_emlrtRSI;
    a = muDoubleScalarMin(1.0, 1.0 / (dNorm * dNorm));
    for (i = 0; i < 50; i++) {
      v[i] += a * d[i];
    }

    (*numIter)++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  st.site = &r_emlrtRSI;
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
  st.site = &s_emlrtRSI;
  indexShapeCheck(&st, xError_vec->size[0], iv);
  init = xError_vec->size[0];
  xError_vec->size[0] = i;
  emxEnsureCapacity_real_T(sp, xError_vec, init, &i_emlrtRTEI);
}

/* End of code generation (logInteriorPoint_conjgrad_rt.c) */
