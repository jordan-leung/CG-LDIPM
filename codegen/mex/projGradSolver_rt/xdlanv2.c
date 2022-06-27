/*
 * xdlanv2.c
 *
 * Code generation for function 'xdlanv2'
 *
 */

/* Include files */
#include "xdlanv2.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Function Definitions */
void xdlanv2(real_T *a, real_T *b, real_T *c, real_T *d, real_T *rt1r, real_T
             *rt1i, real_T *rt2r, real_T *rt2i, real_T *cs, real_T *sn)
{
  real_T bcmax;
  real_T bcmis;
  real_T p;
  real_T scale;
  real_T tau;
  real_T z;
  int32_T b_b;
  int32_T b_c;
  if (*c == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (*b == 0.0) {
    *cs = 0.0;
    *sn = 1.0;
    bcmax = *d;
    *d = *a;
    *a = bcmax;
    *b = -*c;
    *c = 0.0;
  } else {
    tau = *a - *d;
    if ((tau == 0.0) && ((*b < 0.0) != (*c < 0.0))) {
      *cs = 1.0;
      *sn = 0.0;
    } else {
      p = 0.5 * tau;
      bcmis = muDoubleScalarAbs(*b);
      scale = muDoubleScalarAbs(*c);
      bcmax = muDoubleScalarMax(bcmis, scale);
      if (!(*b < 0.0)) {
        b_b = 1;
      } else {
        b_b = -1;
      }

      if (!(*c < 0.0)) {
        b_c = 1;
      } else {
        b_c = -1;
      }

      bcmis = muDoubleScalarMin(bcmis, scale) * (real_T)b_b * (real_T)b_c;
      scale = muDoubleScalarMax(muDoubleScalarAbs(p), bcmax);
      z = p / scale * p + bcmax / scale * bcmis;
      if (z >= 8.8817841970012523E-16) {
        *a = muDoubleScalarSqrt(scale) * muDoubleScalarSqrt(z);
        if (p < 0.0) {
          *a = -*a;
        }

        z = p + *a;
        *a = *d + z;
        *d -= bcmax / z * bcmis;
        tau = muDoubleScalarHypot(*c, z);
        *cs = z / tau;
        *sn = *c / tau;
        *b -= *c;
        *c = 0.0;
      } else {
        bcmis = *b + *c;
        tau = muDoubleScalarHypot(bcmis, tau);
        *cs = muDoubleScalarSqrt(0.5 * (muDoubleScalarAbs(bcmis) / tau + 1.0));
        if (!(bcmis < 0.0)) {
          b_b = 1;
        } else {
          b_b = -1;
        }

        *sn = -(p / (tau * *cs)) * (real_T)b_b;
        bcmax = *a * *cs + *b * *sn;
        scale = -*a * *sn + *b * *cs;
        z = *c * *cs + *d * *sn;
        bcmis = -*c * *sn + *d * *cs;
        *b = scale * *cs + bcmis * *sn;
        *c = -bcmax * *sn + z * *cs;
        bcmax = 0.5 * ((bcmax * *cs + z * *sn) + (-scale * *sn + bcmis * *cs));
        *a = bcmax;
        *d = bcmax;
        if (*c != 0.0) {
          if (*b != 0.0) {
            if ((*b < 0.0) == (*c < 0.0)) {
              bcmis = muDoubleScalarSqrt(muDoubleScalarAbs(*b));
              z = muDoubleScalarSqrt(muDoubleScalarAbs(*c));
              *a = bcmis * z;
              if (!(*c < 0.0)) {
                p = *a;
              } else {
                p = -*a;
              }

              tau = 1.0 / muDoubleScalarSqrt(muDoubleScalarAbs(*b + *c));
              *a = bcmax + p;
              *d = bcmax - p;
              *b -= *c;
              *c = 0.0;
              scale = bcmis * tau;
              bcmis = z * tau;
              bcmax = *cs * scale - *sn * bcmis;
              *sn = *cs * bcmis + *sn * scale;
              *cs = bcmax;
            }
          } else {
            *b = -*c;
            *c = 0.0;
            bcmax = *cs;
            *cs = -*sn;
            *sn = bcmax;
          }
        }
      }
    }
  }

  *rt1r = *a;
  *rt2r = *d;
  if (*c == 0.0) {
    *rt1i = 0.0;
    *rt2i = 0.0;
  } else {
    *rt1i = muDoubleScalarSqrt(muDoubleScalarAbs(*b)) * muDoubleScalarSqrt
      (muDoubleScalarAbs(*c));
    *rt2i = -*rt1i;
  }
}

/* End of code generation (xdlanv2.c) */
