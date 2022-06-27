/*
 * _coder_projGradSolver_rt_api.c
 *
 * Code generation for function '_coder_projGradSolver_rt_api'
 *
 */

/* Include files */
#include "_coder_projGradSolver_rt_api.h"
#include "projGradSolver_rt.h"
#include "projGradSolver_rt_data.h"
#include "projGradSolver_rt_emxutil.h"
#include "projGradSolver_rt_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo i_emlrtRTEI = { 1, /* lineNo */
  1,                                   /* colNo */
  "_coder_projGradSolver_rt_api",      /* fName */
  ""                                   /* pName */
};

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[100];
static const mxArray *b_emlrt_marshallOut(const real_T u);
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *f, const
  char_T *identifier))[10];
static const mxArray *c_emlrt_marshallOut(const emxArray_real_T *u);
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[10];
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *MaxIter,
  const char_T *identifier);
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *H, const
  char_T *identifier))[100];
static const mxArray *emlrt_marshallOut(const real_T u[10]);
static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[100];
static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[10];
static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[100]
{
  real_T (*y)[100];
  y = g_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static const mxArray *b_emlrt_marshallOut(const real_T u)
{
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m);
  return y;
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *f, const
  char_T *identifier))[10]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[10];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(f), &thisId);
  emlrtDestroyArray(&f);
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
  emlrtMsgIdentifier *parentId))[10]
{
  real_T (*y)[10];
  y = h_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *MaxIter,
  const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(MaxIter), &thisId);
  emlrtDestroyArray(&MaxIter);
  return y;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *H, const
  char_T *identifier))[100]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[100];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(H), &thisId);
  emlrtDestroyArray(&H);
  return y;
}
  static const mxArray *emlrt_marshallOut(const real_T u[10])
{
  static const int32_T iv[1] = { 0 };

  static const int32_T iv1[1] = { 10 };

  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(1, &iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, iv1, 1);
  emlrtAssign(&y, m);
  return y;
}

static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[100]
{
  static const int32_T dims[2] = { 10, 10 };

  real_T (*ret)[100];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[100])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[10]
{
  static const int32_T dims[1] = { 10 };

  real_T (*ret)[10];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[10])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void projGradSolver_rt_api(const mxArray * const prhs[8], int32_T nlhs, const
  mxArray *plhs[4])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  emxArray_real_T *xError_vec;
  real_T (*H)[100];
  real_T (*f)[10];
  real_T (*x)[10];
  real_T (*x0)[10];
  real_T (*xOpt)[10];
  real_T (*xl)[10];
  real_T (*xu)[10];
  real_T MaxIter;
  real_T execTime;
  real_T iterCount;
  real_T xTol;
  st.tls = emlrtRootTLSGlobal;
  x = (real_T (*)[10])mxMalloc(sizeof(real_T [10]));
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &xError_vec, 1, &i_emlrtRTEI, true);

  /* Marshall function inputs */
  H = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "H");
  f = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "f");
  x0 = c_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "x0");
  xl = c_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "xl");
  xu = c_emlrt_marshallIn(&st, emlrtAlias(prhs[4]), "xu");
  MaxIter = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "MaxIter");
  xOpt = c_emlrt_marshallIn(&st, emlrtAlias(prhs[6]), "xOpt");
  xTol = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "xTol");

  /* Invoke the target function */
  projGradSolver_rt(&st, *H, *f, *x0, *xl, *xu, MaxIter, *xOpt, xTol, *x,
                    &iterCount, xError_vec, &execTime);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*x);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(iterCount);
  }

  if (nlhs > 2) {
    xError_vec->canFreeData = false;
    plhs[2] = c_emlrt_marshallOut(xError_vec);
  }

  emxFree_real_T(&xError_vec);
  if (nlhs > 3) {
    plhs[3] = b_emlrt_marshallOut(execTime);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_projGradSolver_rt_api.c) */
