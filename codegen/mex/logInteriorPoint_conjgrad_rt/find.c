/*
 * find.c
 *
 * Code generation for function 'find'
 *
 */

/* Include files */
#include "find.h"
#include "indexShapeCheck.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo vb_emlrtRSI = { 144,/* lineNo */
  "eml_find",                          /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/elmat/find.m"/* pathName */
};

static emlrtRSInfo wb_emlrtRSI = { 402,/* lineNo */
  "find_first_indices",                /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/lib/matlab/elmat/find.m"/* pathName */
};

/* Function Definitions */
void eml_find(const emlrtStack *sp, const boolean_T x[20], int32_T i_data[],
              int32_T i_size[1])
{
  emlrtStack b_st;
  emlrtStack st;
  int32_T iv[2];
  int32_T idx;
  int32_T ii;
  boolean_T exitg1;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &vb_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  idx = 0;
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii < 20)) {
    if (x[ii]) {
      idx++;
      i_data[idx - 1] = ii + 1;
      if (idx >= 20) {
        exitg1 = true;
      } else {
        ii++;
      }
    } else {
      ii++;
    }
  }

  if (1 > idx) {
    idx = 0;
  }

  iv[0] = 1;
  iv[1] = idx;
  b_st.site = &wb_emlrtRSI;
  indexShapeCheck(&b_st, 20, iv);
  i_size[0] = idx;
}

/* End of code generation (find.c) */
