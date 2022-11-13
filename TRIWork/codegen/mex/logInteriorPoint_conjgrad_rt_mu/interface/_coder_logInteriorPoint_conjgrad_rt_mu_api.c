/*
 * _coder_logInteriorPoint_conjgrad_rt_mu_api.c
 *
 * Code generation for function '_coder_logInteriorPoint_conjgrad_rt_mu_api'
 *
 */

/* Include files */
#include "_coder_logInteriorPoint_conjgrad_rt_mu_api.h"
#include "logInteriorPoint_conjgrad_rt_mu.h"
#include "logInteriorPoint_conjgrad_rt_mu_data.h"
#include "logInteriorPoint_conjgrad_rt_mu_mexutil.h"
#include "logInteriorPoint_conjgrad_rt_mu_types.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static const mxArray *b_emlrt_marshallOut(const real_T u[120]);
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *W, const
  char_T *identifier))[14400];
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[14400];
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *c, const
  char_T *identifier))[120];
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[120];
static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *mu_f,
  const char_T *identifier);
static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[14400];
static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[120];
static real_T l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */
static const mxArray *b_emlrt_marshallOut(const real_T u[120])
{
  static const int32_T iv[1] = { 0 };

  static const int32_T iv1[1] = { 120 };

  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(1, &iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, iv1, 1);
  emlrtAssign(&y, m);
  return y;
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *W, const
  char_T *identifier))[14400]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[14400];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(W), &thisId);
  emlrtDestroyArray(&W);
  return y;
}
  static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[14400]
{
  real_T (*y)[14400];
  y = j_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *c, const
  char_T *identifier))[120]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[120];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(c), &thisId);
  emlrtDestroyArray(&c);
  return y;
}
  static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[120]
{
  real_T (*y)[120];
  y = k_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *mu_f,
  const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(mu_f), &thisId);
  emlrtDestroyArray(&mu_f);
  return y;
}

static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = l_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[14400]
{
  static const int32_T dims[2] = { 120, 120 };

  real_T (*ret)[14400];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[14400])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[120]
{
  static const int32_T dims[1] = { 120 };

  real_T (*ret)[120];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[120])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void d_logInteriorPoint_conjgrad_rt_(e_logInteriorPoint_conjgrad_rt_ *SD, const
  mxArray * const prhs[10], int32_T nlhs, const mxArray *plhs[5])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  const mxArray *prhs_copy_idx_2;
  real_T (*Aineq)[14400];
  real_T (*W)[14400];
  real_T (*bineq)[120];
  real_T (*c)[120];
  real_T (*v0)[120];
  real_T (*x)[120];
  real_T execTime;
  real_T maxCGIter;
  real_T maxIter;
  real_T mu;
  real_T mu_0;
  real_T mu_f;
  real_T numIter;
  real_T preCondFlag;
  real_T totalCGIters;
  st.tls = emlrtRootTLSGlobal;
  x = (real_T (*)[120])mxMalloc(sizeof(real_T [120]));
  prhs_copy_idx_2 = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  W = c_emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "W");
  c = e_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "c");
  Aineq = c_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_2), "Aineq");
  bineq = e_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "bineq");
  mu_f = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "mu_f");
  mu_0 = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "mu_0");
  v0 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[6]), "v0");
  maxIter = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "maxIter");
  maxCGIter = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "maxCGIter");
  preCondFlag = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "preCondFlag");

  /* Invoke the target function */
  logInteriorPoint_conjgrad_rt_mu(SD, &st, *W, *c, *Aineq, *bineq, mu_f, mu_0,
    *v0, maxIter, maxCGIter, preCondFlag, *x, &mu, &execTime, &numIter,
    &totalCGIters);

  /* Marshall function outputs */
  plhs[0] = b_emlrt_marshallOut(*x);
  if (nlhs > 1) {
    plhs[1] = emlrt_marshallOut(mu);
  }

  if (nlhs > 2) {
    plhs[2] = emlrt_marshallOut(execTime);
  }

  if (nlhs > 3) {
    plhs[3] = emlrt_marshallOut(numIter);
  }

  if (nlhs > 4) {
    plhs[4] = emlrt_marshallOut(totalCGIters);
  }
}

/* End of code generation (_coder_logInteriorPoint_conjgrad_rt_mu_api.c) */
