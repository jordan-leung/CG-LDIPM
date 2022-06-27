/*
 * logInteriorPoint_rt_mexutil.c
 *
 * Code generation for function 'logInteriorPoint_rt_mexutil'
 *
 */

/* Include files */
#include "logInteriorPoint_rt_mexutil.h"
#include "rt_nonfinite.h"

/* Function Definitions */
const mxArray *emlrt_marshallOut(const real_T u)
{
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m);
  return y;
}

/* End of code generation (logInteriorPoint_rt_mexutil.c) */
