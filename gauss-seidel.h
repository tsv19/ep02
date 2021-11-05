#ifndef __GAUSS_SEIDEL__
#define __GAUSS_SEIDEL__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <matheval.h>
#include <math.h>
#include "utils.h"


// metodo de gauss-seidel
double *gauss_seidel(SisLinear_T *sl, unsigned int max_iter, double epsilon);
void gauss_seidel_otimizado(double **diagonais, int n_diagonais, unsigned int max_iter, double epsilon, int n_matriz, double *termos_indep);

// 
// double *gauss_seidel_sem_vetores(SisLinear_T *sl, unsigned int max_iter, double epsilon);

// 
double pega_coeficiente(SisLinear_T *sl, int linha, int coluna);

// 
double *gauss_seidel_teste(SisLinear_T *sl, unsigned int max_iter, double epsilon);

#endif