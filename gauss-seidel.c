
#include "gauss-seidel.h"

// // supondo NAO pivotamento total
// void retrosub(SisLinear_T *sl){
// 	double *incognitas = calloc(sl->n, sizeof(double));
// 	for (int i = sl->n - 1; i >= 0; i--){
// 		for (int j = sl->n -1; j >=0; j--){
// 			incognitas[i] += sl->termos_indep[i] / sl->coeficientes[i][j];
// 			for (int j = i; j < sl->n; j++){
// 				incognitas[i] += sl->coeficientes[k][j]
// 			}
// 		}
// 	}
// }


// void gauss_seidel(SisLinear_T *sl){
// 	/* para cada linha a partir da primeira */
// 	for (int i = 0; i < sl->n; i++){
// 		// uint iPivo = encontraMax(A, i);
// 		// if ( i != iPivo )
// 		// 	trocaLinha( A, b, i, iPivo );

// 		for (int k = i + 1; k < sl->n; k++){
// 			double m = sl->coeficientes[k][i] / sl->coeficientes[i][i];
// 			sl->coeficientes[k][i] = 0.0;
// 			for (int j = i + 1; j < sl->n; j++)
// 				sl->coeficientes[k][j] -= sl->coeficientes[i][j] * m;
// 			sl->termos_indep[k] -= sl->termos_indep[i] * m;
// 		}
// 	}
// }


void gauss_seidel(SisLinear_T *sl, unsigned int max_iter, double epsilon){ //retorna int que indica erro
    double itime, etime, time;
    time = 0;

	// Contem aproximacao inicial
    double *incognitas = malloc(sl->n * sizeof(double));
	for (int l = 0; l < sl->n; l++)
		incognitas[l] = 1;

    for (int iter = 0; iter < max_iter; iter++){
        itime = timestamp(); //tempo inicial
		for (int i = 0; i < sl->n; i++){
	        double resultado = sl->termos_indep[i];
			for (int j = 0; j < sl->n; j++)
				// somar ou subtrair??
				resultado -= sl->coeficientes[i][j] * incognitas[i];
			incognitas[i] = resultado / sl->coeficientes[i][i];
        }

        etime = timestamp(); //tempo final
        time += etime - itime;
        // sl->residuos[iter] = Residuo(ls, cf, res);
    }
    // *gs = (time-(*res))/maxIter;
    // *res = (*res)/maxIter;
}