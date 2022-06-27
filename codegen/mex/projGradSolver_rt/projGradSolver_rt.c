/*
 * projGradSolver_rt.c
 *
 * Code generation for function 'projGradSolver_rt'
 *
 */

/* Include files */
#include "projGradSolver_rt.h"
#include "eig.h"
#include "indexShapeCheck.h"
#include "projGradSolver_rt_data.h"
#include "projGradSolver_rt_emxutil.h"
#include "projGradSolver_rt_types.h"
#include "relop.h"
#include "rt_nonfinite.h"
#include "tic.h"
#include "toc.h"
#include "mwmathutil.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 48,    /* lineNo */
  "projGradSolver_rt",                 /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/projGradSolver_rt.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 47,  /* lineNo */
  "projGradSolver_rt",                 /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/projGradSolver_rt.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 21,  /* lineNo */
  "projGradSolver_rt",                 /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/projGradSolver_rt.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 13,  /* lineNo */
  "projGradSolver_rt",                 /* fcnName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/projGradSolver_rt.m"/* pathName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  48,                                  /* lineNo */
  27,                                  /* colNo */
  "xError_vec",                        /* aName */
  "projGradSolver_rt",                 /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/projGradSolver_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  48,                                  /* lineNo */
  25,                                  /* colNo */
  "xError_vec",                        /* aName */
  "projGradSolver_rt",                 /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/projGradSolver_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  20,                                  /* lineNo */
  12,                                  /* colNo */
  "xError_vec",                        /* aName */
  "projGradSolver_rt",                 /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/projGradSolver_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo emlrtDCI = { 18,    /* lineNo */
  1,                                   /* colNo */
  "projGradSolver_rt",                 /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/projGradSolver_rt.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = { 18,  /* lineNo */
  1,                                   /* colNo */
  "projGradSolver_rt",                 /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/projGradSolver_rt.m",/* pName */
  4                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  45,                                  /* lineNo */
  5,                                   /* colNo */
  "xError_vec",                        /* aName */
  "projGradSolver_rt",                 /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/projGradSolver_rt.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo g_emlrtRTEI = { 18,/* lineNo */
  1,                                   /* colNo */
  "projGradSolver_rt",                 /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/projGradSolver_rt.m"/* pName */
};

static emlrtRTEInfo h_emlrtRTEI = { 48,/* lineNo */
  1,                                   /* colNo */
  "projGradSolver_rt",                 /* fName */
  "/home/jordanleung/Documents/GitHub/optimizationFunctions/projGradSolver_rt.m"/* pName */
};

/* Function Definitions */
void projGradSolver_rt(const emlrtStack *sp, const real_T H[100], const real_T
  f[10], const real_T x0[10], const real_T xl[10], const real_T xu[10], real_T
  MaxIter, const real_T xOpt[10], real_T xTol, real_T x[10], real_T *iterCount,
  emxArray_real_T *xError_vec, real_T *execTime)
{
  emlrtStack st;
  creal_T eigVec[10];
  creal_T L;
  creal_T ell;
  real_T xStar[10];
  real_T absxk;
  real_T re;
  real_T scale;
  real_T t;
  real_T xError;
  int32_T iv[2];
  int32_T i;
  int32_T k;
  boolean_T p;
  st.prev = sp;
  st.tls = sp->tls;

  /*  * May eventually be updated to solve for general linear constraints, but */
  /*  current this code solves for box constraints only */
  memcpy(&x[0], &x0[0], 10U * sizeof(real_T));

  /*  Initialize x */
  *iterCount = 0.0;

  /*  Find the Lipshitz constant L and set step-size */
  st.site = &f_emlrtRSI;
  eig(&st, H, eigVec);
  L = eigVec[0];
  for (k = 0; k < 9; k++) {
    if (muDoubleScalarIsNaN(eigVec[k + 1].re) || muDoubleScalarIsNaN(eigVec[k +
         1].im)) {
      p = false;
    } else if (muDoubleScalarIsNaN(L.re) || muDoubleScalarIsNaN(L.im)) {
      p = true;
    } else {
      absRelopProxies(L, eigVec[k + 1], &xError, &scale);
      p = (xError < scale);
    }

    if (p) {
      L = eigVec[k + 1];
    }
  }

  ell = eigVec[0];
  for (k = 0; k < 9; k++) {
    if (muDoubleScalarIsNaN(eigVec[k + 1].re) || muDoubleScalarIsNaN(eigVec[k +
         1].im)) {
      p = false;
    } else if (muDoubleScalarIsNaN(ell.re) || muDoubleScalarIsNaN(ell.im)) {
      p = true;
    } else {
      absRelopProxies(ell, eigVec[k + 1], &xError, &scale);
      p = (xError > scale);
    }

    if (p) {
      ell = eigVec[k + 1];
    }
  }

  absxk = L.re + ell.re;
  t = L.im + ell.im;
  if (t == 0.0) {
    re = 2.0 / absxk;
  } else if (absxk == 0.0) {
    re = 0.0;
  } else {
    xError = muDoubleScalarAbs(absxk);
    scale = muDoubleScalarAbs(t);
    if (xError > scale) {
      xError = t / absxk;
      re = (2.0 + xError * 0.0) / (absxk + xError * t);
    } else if (scale == xError) {
      if (absxk > 0.0) {
        absxk = 0.5;
      } else {
        absxk = -0.5;
      }

      if (t > 0.0) {
        t = 0.5;
      } else {
        t = -0.5;
      }

      re = (2.0 * absxk + 0.0 * t) / xError;
    } else {
      xError = absxk / t;
      re = xError * 2.0 / (t + xError * absxk);
    }
  }

  if (!(MaxIter >= 0.0)) {
    emlrtNonNegativeCheckR2012b(MaxIter, &b_emlrtDCI, sp);
  }

  xError = (int32_T)muDoubleScalarFloor(MaxIter);
  if (MaxIter != xError) {
    emlrtIntegerCheckR2012b(MaxIter, &emlrtDCI, sp);
  }

  i = xError_vec->size[0];
  xError_vec->size[0] = (int32_T)MaxIter;
  emxEnsureCapacity_real_T(sp, xError_vec, i, &g_emlrtRTEI);
  if (MaxIter != xError) {
    emlrtIntegerCheckR2012b(MaxIter, &emlrtDCI, sp);
  }

  k = (int32_T)MaxIter;
  for (i = 0; i < k; i++) {
    xError_vec->data[i] = 0.0;
  }

  xError = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 10; k++) {
    absxk = muDoubleScalarAbs(x0[k] - xOpt[k]);
    if (absxk > scale) {
      t = scale / absxk;
      xError = xError * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      xError += t * t;
    }
  }

  xError = scale * muDoubleScalarSqrt(xError);
  if (1 > (int32_T)MaxIter) {
    emlrtDynamicBoundsCheckR2012b(1, 1, (int32_T)MaxIter, &c_emlrtBCI, sp);
  }

  xError_vec->data[0] = xError;
  st.site = &c_emlrtRSI;
  tic(&st);
  while ((*iterCount < MaxIter) && (xError > xTol)) {
    (*iterCount)++;

    /*  Current gradient direction */
    /*  Compute temporary step */
    for (i = 0; i < 10; i++) {
      xError = 0.0;
      for (k = 0; k < 10; k++) {
        xError += H[i + 10 * k] * x[k];
      }

      xStar[i] = x[i] - re * (xError + f[i]);
    }

    /*  Project any components exceeding the constraints back onto the box */
    /*  constraint */
    for (k = 0; k < 10; k++) {
      xError = xStar[k];
      scale = xu[k];
      if (xError > scale) {
        x[k] = scale;
      } else if (xError < xl[k]) {
        x[k] = xl[k];
      } else {
        x[k] = xError;
      }

      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }

    /*  Check convergence criteria */
    xError = 0.0;
    scale = 3.3121686421112381E-170;
    for (k = 0; k < 10; k++) {
      absxk = muDoubleScalarAbs(xOpt[k] - x[k]);
      if (absxk > scale) {
        t = scale / absxk;
        xError = xError * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        xError += t * t;
      }
    }

    xError = scale * muDoubleScalarSqrt(xError);
    if (((int32_T)(*iterCount + 1.0) < 1) || ((int32_T)(*iterCount + 1.0) >
         xError_vec->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(*iterCount + 1.0), 1,
        xError_vec->size[0], &d_emlrtBCI, sp);
    }

    xError_vec->data[(int32_T)(*iterCount + 1.0) - 1] = xError;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  st.site = &b_emlrtRSI;
  *execTime = toc(&st);
  if (1 > xError_vec->size[0]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, xError_vec->size[0], &b_emlrtBCI, sp);
  }

  if (((int32_T)(*iterCount + 1.0) < 1) || ((int32_T)(*iterCount + 1.0) >
       xError_vec->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(*iterCount + 1.0), 1,
      xError_vec->size[0], &emlrtBCI, sp);
  }

  iv[0] = 1;
  iv[1] = (int32_T)(*iterCount + 1.0);
  st.site = &emlrtRSI;
  indexShapeCheck(&st, xError_vec->size[0], iv);
  i = xError_vec->size[0];
  xError_vec->size[0] = (int32_T)(*iterCount + 1.0);
  emxEnsureCapacity_real_T(sp, xError_vec, i, &h_emlrtRTEI);
}

/* End of code generation (projGradSolver_rt.c) */
