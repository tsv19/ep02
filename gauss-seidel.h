#ifndef __GAUSS_SEIDEL__
#define __GAUSS_SEIDEL__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <matheval.h>
#include <math.h>
#include "utils.h"


// 
void gauss_seidel(SisLinear_T *sl, unsigned int max_iter, double epsilon);

#endif