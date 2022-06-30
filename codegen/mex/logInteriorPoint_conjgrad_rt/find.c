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

/* Function Definitions */
void eml_find(const boolean_T x[120], int32_T i_data[], int32_T i_size[1])
{
  int32_T idx;
  int32_T ii;
  boolean_T exitg1;
  idx = 0;
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii < 120)) {
    if (x[ii]) {
      idx++;
      i_data[idx - 1] = ii + 1;
      if (idx >= 120) {
        exitg1 = true;
      } else {
        ii++;
      }
    } else {
      ii++;
    }
  }

  indexShapeCheck();
  if (1 > idx) {
    i_size[0] = 0;
  } else {
    i_size[0] = idx;
  }
}

/* End of code generation (find.c) */
