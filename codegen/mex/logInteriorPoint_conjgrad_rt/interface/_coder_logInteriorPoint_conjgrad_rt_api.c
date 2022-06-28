/*
 * _coder_logInteriorPoint_conjgrad_rt_api.c
 *
 * Code generation for function '_coder_logInteriorPoint_conjgrad_rt_api'
 *
 */

/* Include files */
#include "_coder_logInteriorPoint_conjgrad_rt_api.h"
#include "logInteriorPoint_conjgrad_rt.h"
#include "logInteriorPoint_conjgrad_rt_data.h"
#include "logInteriorPoint_conjgrad_rt_emxutil.h"
#include "logInteriorPoint_conjgrad_rt_mexutil.h"
#include "logInteriorPoint_conjgrad_rt_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo j_emlrtRTEI = { 1, /* lineNo */
  1,                                   /* colNo */
  "_coder_logInteriorPoint_conjgrad_rt_api",/* fName */
  ""                                   /* pName */
};

/* Function Declarations */
static const mxArray *b_emlrt_marshallOut(const real_T u[100]);
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *W, const
  char_T *identifier))[10000];
static const mxArray *c_emlrt_marshallOut(const emxArray_real_T *u);
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[10000];
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *c, const
  char_T *identifier))[100];
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[100];
static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Aineq,
  const char_T *identifier))[5000];
static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[5000];
static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *bineq,
  const char_T *identifier))[50];
static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[50];
static real_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *mu_f,
  const char_T *identifier);
static real_T l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[10000];
static real_T (*o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[100];
static real_T (*p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[5000];
static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[50];
static real_T r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */
static const mxArray *b_emlrt_marshallOut(const real_T u[100])
{
  static const int32_T iv[1] = { 0 };

  static const int32_T iv1[1] = { 100 };

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
  char_T *identifier))[10000]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[10000];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(W), &thisId);
  emlrtDestroyArray(&W);
  return y;
}
  static const mxArray *c_emlrt_marshallOut(const emxArray_real_T *u)
{
  static const int32_T iv[1] = { 0 };

  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(1, &iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u->data[0]);
  emlrtSetDimensions((mxArray *)m, u->size, 1);
  emlrtAssign(&y, m);
  return y;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[10000]
{
  real_T (*y)[10000];
  y = n_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *c,
  const char_T *identifier))[100]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[100];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(c), &thisId);
  emlrtDestroyArray(&c);
  return y;
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[100]
{
  real_T (*y)[100];
  y = o_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Aineq,
  const char_T *identifier))[5000]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[5000];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(Aineq), &thisId);
  emlrtDestroyArray(&Aineq);
  return y;
}

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[5000]
{
  real_T (*y)[5000];
  y = p_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *bineq,
  const char_T *identifier))[50]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[50];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(sp, emlrtAlias(bineq), &thisId);
  emlrtDestroyArray(&bineq);
  return y;
}

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[50]
{
  real_T (*y)[50];
  y = q_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *mu_f,
  const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = l_emlrt_marshallIn(sp, emlrtAlias(mu_f), &thisId);
  emlrtDestroyArray(&mu_f);
  return y;
}

static real_T l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = r_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[10000]
{
  static const int32_T dims[2] = { 100, 100 };

  real_T (*ret)[10000];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[10000])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[100]
{
  static const int32_T dims[1] = { 100 };

  real_T (*ret)[100];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[100])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[5000]
{
  static const int32_T dims[2] = { 50, 100 };

  real_T (*ret)[5000];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[5000])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[50]
{
  static const int32_T dims[1] = { 50 };

  real_T (*ret)[50];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[50])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void d_logInteriorPoint_conjgrad_rt_(c_logInteriorPoint_conjgrad_rtS *SD, const
  mxArray * const prhs[12], int32_T nlhs, const mxArray *plhs[4])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  emxArray_real_T *xError_vec;
  const mxArray *prhs_copy_idx_2;
  real_T (*W)[10000];
  real_T (*Aineq)[5000];
  real_T (*c)[100];
  real_T (*x)[100];
  real_T (*xStar)[100];
  real_T (*bineq)[50];
  real_T (*v0)[50];
  real_T execTime;
  real_T maxCGIter;
  real_T maxIter;
  real_T mu_0;
  real_T mu_f;
  real_T numIter;
  real_T preCondFlag;
  real_T xTol;
  st.tls = emlrtRootTLSGlobal;
  x = (real_T (*)[100])mxMalloc(sizeof(real_T [100]));
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &xError_vec, 1, &j_emlrtRTEI, true);
  prhs_copy_idx_2 = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  W = c_emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "W");
  c = e_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "c");
  Aineq = g_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_2), "Aineq");
  bineq = i_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "bineq");
  mu_f = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "mu_f");
  mu_0 = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "mu_0");
  v0 = i_emlrt_marshallIn(&st, emlrtAlias(prhs[6]), "v0");
  maxIter = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "maxIter");
  maxCGIter = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "maxCGIter");
  preCondFlag = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "preCondFlag");
  xStar = e_emlrt_marshallIn(&st, emlrtAlias(prhs[10]), "xStar");
  xTol = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "xTol");

  /* Invoke the target function */
  logInteriorPoint_conjgrad_rt(SD, &st, *W, *c, *Aineq, *bineq, mu_f, mu_0, *v0,
    maxIter, maxCGIter, preCondFlag, *xStar, xTol, *x, xError_vec, &execTime,
    &numIter);

  /* Marshall function outputs */
  plhs[0] = b_emlrt_marshallOut(*x);
  if (nlhs > 1) {
    xError_vec->canFreeData = false;
    plhs[1] = c_emlrt_marshallOut(xError_vec);
  }

  emxFree_real_T(&xError_vec);
  if (nlhs > 2) {
    plhs[2] = emlrt_marshallOut(execTime);
  }

  if (nlhs > 3) {
    plhs[3] = emlrt_marshallOut(numIter);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_logInteriorPoint_conjgrad_rt_api.c) */
