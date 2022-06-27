/*
 * xzlarfg.c
 *
 * Code generation for function 'xzlarfg'
 *
 */

/* Include files */
#include "xzlarfg.h"
#include "eml_int_forloop_overflow_check.h"
#include "projGradSolver_rt_data.h"
#include "rt_nonfinite.h"
#include "xnrm2.h"
#include "mwmathutil.h"

/* Function Definitions */
real_T xzlarfg(const emlrtStack *sp, int32_T n, real_T *alpha1, real_T x[3])
{
  emlrtStack b_st;
  emlrtStack st;
  real_T beta1;
  real_T tau;
  real_T xnorm;
  int32_T k;
  int32_T knt;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  tau = 0.0;
  if (n > 0) {
    xnorm = b_xnrm2(n - 1, x);
    if (xnorm != 0.0) {
      beta1 = muDoubleScalarHypot(*alpha1, xnorm);
      if (*alpha1 >= 0.0) {
        beta1 = -beta1;
      }

      if (muDoubleScalarAbs(beta1) < 1.0020841800044864E-292) {
        knt = 0;
        do {
          knt++;
          for (k = 2; k <= n; k++) {
            x[k - 1] *= 9.9792015476736E+291;
          }

          beta1 *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while (!(muDoubleScalarAbs(beta1) >= 1.0020841800044864E-292));

        beta1 = muDoubleScalarHypot(*alpha1, b_xnrm2(n - 1, x));
        if (*alpha1 >= 0.0) {
          beta1 = -beta1;
        }

        tau = (beta1 - *alpha1) / beta1;
        xnorm = 1.0 / (*alpha1 - beta1);
        for (k = 2; k <= n; k++) {
          x[k - 1] *= xnorm;
        }

        st.site = &bb_emlrtRSI;
        if ((1 <= knt) && (knt > 2147483646)) {
          b_st.site = &m_emlrtRSI;
          check_forloop_overflow_error(&b_st);
        }

        for (k = 0; k < knt; k++) {
          beta1 *= 1.0020841800044864E-292;
        }

        *alpha1 = beta1;
      } else {
        tau = (beta1 - *alpha1) / beta1;
        xnorm = 1.0 / (*alpha1 - beta1);
        for (k = 2; k <= n; k++) {
          x[k - 1] *= xnorm;
        }

        *alpha1 = beta1;
      }
    }
  }

  return tau;
}

/* End of code generation (xzlarfg.c) */
