/*
 * xdhseqr.c
 *
 * Code generation for function 'xdhseqr'
 *
 */

/* Include files */
#include "xdhseqr.h"
#include "rt_nonfinite.h"
#include "xdlanv2.h"
#include "xrot.h"
#include "xzlarfg.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo cc_emlrtRSI = { 263,/* lineNo */
  "eml_dlahqr",                        /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xdhseqr.m"/* pathName */
};

static emlrtRSInfo gc_emlrtRSI = { 359,/* lineNo */
  "eml_dlahqr",                        /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xdhseqr.m"/* pathName */
};

static emlrtRSInfo hc_emlrtRSI = { 365,/* lineNo */
  "eml_dlahqr",                        /* fcnName */
  "/home/jordanleung/MATLAB2020b/toolbox/eml/eml/+coder/+internal/+reflapack/xdhseqr.m"/* pathName */
};

/* Function Definitions */
int32_T eml_dlahqr(const emlrtStack *sp, real_T h[100])
{
  emlrtStack st;
  real_T v[3];
  real_T aa;
  real_T ab;
  real_T ba;
  real_T bb;
  real_T rt1r;
  real_T rt2r;
  real_T s;
  real_T s_tmp;
  real_T tst;
  int32_T L;
  int32_T b;
  int32_T b_k;
  int32_T hoffset;
  int32_T i;
  int32_T info;
  int32_T its;
  int32_T j;
  int32_T k;
  int32_T m;
  int32_T nr;
  int32_T sum1_tmp;
  boolean_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  boolean_T goto150;
  st.prev = sp;
  st.tls = sp->tls;
  info = 0;
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 0.0;
  for (j = 0; j < 7; j++) {
    b = j + 10 * j;
    h[b + 2] = 0.0;
    h[b + 3] = 0.0;
  }

  h[79] = 0.0;
  i = 9;
  exitg1 = false;
  while ((!exitg1) && (i + 1 >= 1)) {
    L = 1;
    goto150 = false;
    its = 0;
    exitg2 = false;
    while ((!exitg2) && (its < 301)) {
      k = i;
      exitg3 = false;
      while ((!exitg3) && (k + 1 > L)) {
        b = k + 10 * (k - 1);
        ba = muDoubleScalarAbs(h[b]);
        if (ba <= 1.0020841800044864E-291) {
          exitg3 = true;
        } else {
          hoffset = k + 10 * k;
          bb = muDoubleScalarAbs(h[hoffset]);
          tst = muDoubleScalarAbs(h[b - 1]) + bb;
          if (tst == 0.0) {
            if (k - 1 >= 1) {
              tst = muDoubleScalarAbs(h[(k + 10 * (k - 2)) - 1]);
            }

            if (k + 2 <= 10) {
              tst += muDoubleScalarAbs(h[hoffset + 1]);
            }
          }

          if (ba <= 2.2204460492503131E-16 * tst) {
            tst = muDoubleScalarAbs(h[hoffset - 1]);
            if (ba > tst) {
              ab = ba;
              ba = tst;
            } else {
              ab = tst;
            }

            tst = muDoubleScalarAbs(h[b - 1] - h[hoffset]);
            if (bb > tst) {
              aa = bb;
              bb = tst;
            } else {
              aa = tst;
            }

            s = aa + ab;
            if (ba * (ab / s) <= muDoubleScalarMax(1.0020841800044864E-291,
                 2.2204460492503131E-16 * (bb * (aa / s)))) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }
      }

      L = k + 1;
      if (k + 1 > 1) {
        h[k + 10 * (k - 1)] = 0.0;
      }

      if (k + 1 >= i) {
        goto150 = true;
        exitg2 = true;
      } else {
        if (its == 10) {
          hoffset = k + 10 * k;
          s = muDoubleScalarAbs(h[hoffset + 1]) + muDoubleScalarAbs(h[(k + 10 *
            (k + 1)) + 2]);
          tst = 0.75 * s + h[hoffset];
          aa = -0.4375 * s;
          ab = s;
          bb = tst;
        } else if (its == 20) {
          s = muDoubleScalarAbs(h[i + 10 * (i - 1)]) + muDoubleScalarAbs(h[(i +
            10 * (i - 2)) - 1]);
          tst = 0.75 * s + h[i + 10 * i];
          aa = -0.4375 * s;
          ab = s;
          bb = tst;
        } else {
          hoffset = i + 10 * (i - 1);
          tst = h[hoffset - 1];
          ab = h[hoffset];
          aa = h[(i + 10 * i) - 1];
          bb = h[i + 10 * i];
        }

        s = ((muDoubleScalarAbs(tst) + muDoubleScalarAbs(aa)) +
             muDoubleScalarAbs(ab)) + muDoubleScalarAbs(bb);
        if (s == 0.0) {
          rt1r = 0.0;
          ba = 0.0;
          rt2r = 0.0;
          bb = 0.0;
        } else {
          tst /= s;
          ab /= s;
          aa /= s;
          bb /= s;
          ba = (tst + bb) / 2.0;
          tst = (tst - ba) * (bb - ba) - aa * ab;
          ab = muDoubleScalarSqrt(muDoubleScalarAbs(tst));
          if (tst >= 0.0) {
            rt1r = ba * s;
            rt2r = rt1r;
            ba = ab * s;
            bb = -ba;
          } else {
            rt1r = ba + ab;
            rt2r = ba - ab;
            if (muDoubleScalarAbs(rt1r - bb) <= muDoubleScalarAbs(rt2r - bb)) {
              rt1r *= s;
              rt2r = rt1r;
            } else {
              rt2r *= s;
              rt1r = rt2r;
            }

            ba = 0.0;
            bb = 0.0;
          }
        }

        m = i - 1;
        exitg3 = false;
        while ((!exitg3) && (m >= k + 1)) {
          hoffset = m + 10 * (m - 1);
          tst = h[hoffset];
          s_tmp = h[hoffset - 1];
          ab = s_tmp - rt2r;
          s = (muDoubleScalarAbs(ab) + muDoubleScalarAbs(bb)) +
            muDoubleScalarAbs(tst);
          aa = tst / s;
          hoffset = m + 10 * m;
          v[0] = (aa * h[hoffset - 1] + (s_tmp - rt1r) * (ab / s)) - ba * (bb /
            s);
          tst = h[hoffset];
          v[1] = aa * (((s_tmp + tst) - rt1r) - rt2r);
          v[2] = aa * h[hoffset + 1];
          s = (muDoubleScalarAbs(v[0]) + muDoubleScalarAbs(v[1])) +
            muDoubleScalarAbs(v[2]);
          v[0] /= s;
          v[1] /= s;
          v[2] /= s;
          if (m == k + 1) {
            exitg3 = true;
          } else {
            b = m + 10 * (m - 2);
            if (muDoubleScalarAbs(h[b - 1]) * (muDoubleScalarAbs(v[1]) +
                 muDoubleScalarAbs(v[2])) <= 2.2204460492503131E-16 *
                muDoubleScalarAbs(v[0]) * ((muDoubleScalarAbs(h[b - 2]) +
                  muDoubleScalarAbs(s_tmp)) + muDoubleScalarAbs(tst))) {
              exitg3 = true;
            } else {
              m--;
            }
          }
        }

        for (b_k = m; b_k <= i; b_k++) {
          hoffset = (i - b_k) + 2;
          nr = muIntScalarMin_sint32(3, hoffset);
          if (b_k > m) {
            hoffset = (b_k + 10 * (b_k - 2)) - 1;
            for (j = 0; j < nr; j++) {
              v[j] = h[j + hoffset];
            }
          }

          tst = v[0];
          st.site = &cc_emlrtRSI;
          bb = xzlarfg(&st, nr, &tst, v);
          v[0] = tst;
          if (b_k > m) {
            h[(b_k + 10 * (b_k - 2)) - 1] = tst;
            b = b_k + 10 * (b_k - 2);
            h[b] = 0.0;
            if (b_k < i) {
              h[b + 1] = 0.0;
            }
          } else {
            if (m > k + 1) {
              h[(b_k + 10 * (b_k - 2)) - 1] *= 1.0 - bb;
            }
          }

          rt1r = v[1];
          tst = bb * v[1];
          if (nr == 3) {
            s = v[2];
            ba = bb * v[2];
            for (j = b_k; j < 11; j++) {
              hoffset = b_k + 10 * (j - 1);
              aa = (h[hoffset - 1] + rt1r * h[hoffset]) + s * h[hoffset + 1];
              h[hoffset - 1] -= aa * bb;
              h[hoffset] -= aa * tst;
              h[hoffset + 1] -= aa * ba;
            }

            hoffset = b_k + 3;
            b = i + 1;
            b = muIntScalarMin_sint32(hoffset, b);
            for (j = 0; j < b; j++) {
              hoffset = j + 10 * (b_k - 1);
              ab = h[hoffset];
              nr = j + 10 * b_k;
              sum1_tmp = j + 10 * (b_k + 1);
              aa = (ab + rt1r * h[nr]) + s * h[sum1_tmp];
              h[hoffset] = ab - aa * bb;
              h[nr] -= aa * tst;
              h[sum1_tmp] -= aa * ba;
            }
          } else {
            if (nr == 2) {
              for (j = b_k; j < 11; j++) {
                hoffset = b_k + 10 * (j - 1);
                ab = h[hoffset - 1];
                aa = ab + rt1r * h[hoffset];
                h[hoffset - 1] = ab - aa * bb;
                h[hoffset] -= aa * tst;
              }

              for (j = 0; j <= i; j++) {
                hoffset = j + 10 * (b_k - 1);
                nr = j + 10 * b_k;
                aa = h[hoffset] + rt1r * h[nr];
                h[hoffset] -= aa * bb;
                h[nr] -= aa * tst;
              }
            }
          }
        }

        its++;
      }
    }

    if (!goto150) {
      info = i + 1;
      exitg1 = true;
    } else {
      if ((L != i + 1) && (L == i)) {
        b = i + 10 * i;
        rt1r = h[b - 1];
        hoffset = i + 10 * (i - 1);
        s = h[hoffset];
        tst = h[b];
        xdlanv2(&h[(i + 10 * (i - 1)) - 1], &rt1r, &s, &tst, &ab, &aa, &ba, &bb,
                &s_tmp, &rt2r);
        h[b - 1] = rt1r;
        h[hoffset] = s;
        h[b] = tst;
        if (10 > i + 1) {
          st.site = &gc_emlrtRSI;
          xrot(&st, 9 - i, h, i + (i + 1) * 10, (i + (i + 1) * 10) + 1, s_tmp,
               rt2r);
        }

        st.site = &hc_emlrtRSI;
        b_xrot(&st, i - 1, h, (i - 1) * 10 + 1, i * 10 + 1, s_tmp, rt2r);
      }

      i = L - 2;
    }
  }

  return info;
}

/* End of code generation (xdhseqr.c) */
