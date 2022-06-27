/*
 * relop.c
 *
 * Code generation for function 'relop'
 *
 */

/* Include files */
#include "relop.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"
#include <math.h>

/* Function Declarations */
static boolean_T iseq(real_T x, real_T y);

/* Function Definitions */
static boolean_T iseq(real_T x, real_T y)
{
  real_T absx;
  int32_T exponent;
  boolean_T p;
  absx = muDoubleScalarAbs(y / 2.0);
  if ((!muDoubleScalarIsInf(absx)) && (!muDoubleScalarIsNaN(absx))) {
    if (absx <= 2.2250738585072014E-308) {
      absx = 4.94065645841247E-324;
    } else {
      frexp(absx, &exponent);
      absx = ldexp(1.0, exponent - 53);
    }
  } else {
    absx = rtNaN;
  }

  if ((muDoubleScalarAbs(y - x) < absx) || (muDoubleScalarIsInf(x) &&
       muDoubleScalarIsInf(y) && ((x > 0.0) == (y > 0.0)))) {
    p = true;
  } else {
    p = false;
  }

  return p;
}

void absRelopProxies(const creal_T a, const creal_T b, real_T *x, real_T *y)
{
  real_T Ma;
  real_T absai;
  real_T absbi;
  real_T ma;
  real_T mb;
  boolean_T SCALEA;
  boolean_T SCALEB;
  ma = muDoubleScalarAbs(a.re);
  if ((ma > 8.9884656743115785E+307) || (muDoubleScalarAbs(a.im) >
       8.9884656743115785E+307)) {
    SCALEA = true;
  } else {
    SCALEA = false;
  }

  mb = muDoubleScalarAbs(b.re);
  if ((mb > 8.9884656743115785E+307) || (muDoubleScalarAbs(b.im) >
       8.9884656743115785E+307)) {
    SCALEB = true;
  } else {
    SCALEB = false;
  }

  if (SCALEA || SCALEB) {
    *x = muDoubleScalarHypot(a.re / 2.0, a.im / 2.0);
    *y = muDoubleScalarHypot(b.re / 2.0, b.im / 2.0);
  } else {
    *x = muDoubleScalarHypot(a.re, a.im);
    *y = muDoubleScalarHypot(b.re, b.im);
  }

  if (iseq(*x, *y)) {
    absai = muDoubleScalarAbs(a.im);
    absbi = muDoubleScalarAbs(b.im);
    if (ma > absai) {
      Ma = ma;
      ma = absai;
    } else {
      Ma = absai;
    }

    if (mb > absbi) {
      absai = mb;
      mb = absbi;
    } else {
      absai = absbi;
    }

    if (Ma > absai) {
      if (ma < mb) {
        *x = Ma - absai;
        *y = (ma / 2.0 + mb / 2.0) / (Ma / 2.0 + absai / 2.0) * (mb - ma);
      } else {
        *x = Ma;
        *y = absai;
      }
    } else if (Ma < absai) {
      if (ma > mb) {
        *y = absai - Ma;
        *x = (ma / 2.0 + mb / 2.0) / (Ma / 2.0 + absai / 2.0) * (ma - mb);
      } else {
        *x = Ma;
        *y = absai;
      }
    } else {
      *x = ma;
      *y = mb;
    }

    if (iseq(*x, *y)) {
      *x = muDoubleScalarAtan2(a.im, a.re);
      *y = muDoubleScalarAtan2(b.im, b.re);
      if (iseq(*x, *y)) {
        if (*x > 0.78539816339744828) {
          if (*x > 2.3561944901923448) {
            *x = -a.im;
            *y = -b.im;
          } else {
            *x = -a.re;
            *y = -b.re;
          }
        } else if (*x > -0.78539816339744828) {
          *x = a.im;
          *y = b.im;
        } else if (*x > -2.3561944901923448) {
          *x = a.re;
          *y = b.re;
        } else {
          *x = -a.im;
          *y = -b.im;
        }

        if (iseq(*x, *y)) {
          *x = 0.0;
          *y = 0.0;
        }
      }
    }
  }
}

/* End of code generation (relop.c) */
