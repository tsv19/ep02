
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
	double * incognitas;
	unsigned int max_iter;			// número máximo de iterações.
	SisLinear_T *sistema_linear;
	
	sistema_linear = read_input(&epsilon, &max_iter);

	gera_sistema_linear(sistema_linear);
	gera_termos_indep(sistema_linear);

    tempo = timestamp();
	
	double c1[4] = {0, 1, 2, 3};
	double c0[5] = {0, 0, 0, 0, 0};
	double d[6] = {1, 4, 7, 10, 13, 16};
	double a0[5] = {0, 2, 4, 6, 8};
	double a1[4] = {1, 1, 1, 1};
      
	double* diagonais[5];
	diagonais[0] = c1;
	diagonais[1] = c0;
	diagonais[2] = d;
	diagonais[3] = a0;
	diagonais[4] = a1;

	gauss_seidel_otimizado(diagonais, 5, max_iter, epsilon, 6, sistema_linear->termos_indep);

	// incognitas = gauss_seidel(sistema_linear, max_iter, epsilon);
	// tempo = timestamp() - tempo;
	// imprime_diagonais(sistema_linear);
	// imprime_termos_indep(sistema_linear);
	// imprime_solucao(incognitas, sistema_linear->n);
	// imprime_tempo(tempo);

	destroi_sistema_linear(sistema_linear);
	return EXIT_SUCCESS;
}