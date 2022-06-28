/*
 * mtimes.h
 *
 * Code generation for function 'mtimes'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void b_mtimes(const real_T A[5000], const real_T B[50], real_T C[100]);
void c_mtimes(const real_T A[5000], const real_T B[100], real_T C[50]);
void d_mtimes(const real_T A[5000], const real_T B[50], real_T C[100]);
void mtimes(const real_T A[10000], const real_T B[100], real_T C[100]);

/* End of code generation (mtimes.h) */
