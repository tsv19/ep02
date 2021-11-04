
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <matheval.h>
#include <math.h>
#include <sys/time.h>
#include "gauss-seidel.h"
#include "utils.h"


int main(){
	double tempo;
	double epsilon;					// tolerância no método de Gauss-Seidel.
	unsigned int max_iter;			// número máximo de iterações.
	SisLinear_T *sistema_linear = aloca_sistema_linear();

	read_input(sistema_linear, &epsilon, &max_iter);

	ajusta_sistema_linear(sistema_linear);
	gera_sistema_linear(sistema_linear);
	gera_termos_indep(sistema_linear);


    tempo = timestamp();
	
	// double *d = {1, 4, 7, 10, 13, 16};
	// double *a0 = {0, 0, 0, 0, 0};
	// double *a1 = {1, 1, 1, 1};
	// double *c0 = {0, 0, 0, 0, 0};
	// double *c1 = {0, 1, 2, 3};

	gauss_seidel(sistema_linear, max_iter, epsilon);
	tempo = timestamp() - tempo;
	imprime_sistema_linear(sistema_linear);
	// imprime_termos_indep(sistema_linear);
	// imprime_solucao(sistema_linear);
	// printf("tempo: %lf\n", tempo);

	destroi_sistema_linear(sistema_linear);
	return EXIT_SUCCESS;
}