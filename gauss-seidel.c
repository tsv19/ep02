
#include "gauss-seidel.h"


double *calcula_residuo(SisLinear_T *sl, double *incognitas){
	double *residuos = malloc(sl->n * sizeof(double));

	for (int i = 0; i < sl->n; i++){
		double res = 0.0;
		for (int j = 0; j < sl->n; j++)
			res += sl->coeficientes[i][j] * incognitas[j];
		residuos[i] = res - sl->termos_indep[i];
	}
	return residuos;
}


double *gauss_seidel(SisLinear_T *sl, unsigned int max_iter, double epsilon){
	// Contem aproximacao inicial
    double *incognitas = calloc(sl->n, sizeof(double));
	// 	criterio de parada quando a norma das incognitas eh menor que epsilon 
	double norma = 1.0 + epsilon;

    for (int iter = 0; norma > epsilon && iter < max_iter; iter++){
		norma = 0.0;
		for (int i = 0; i < sl->n; i++){
	        double resultado = 0;
			for (int j = 0; j < sl->n; j++)
				if (i != j)
					resultado += sl->coeficientes[i][j] * incognitas[j];
			
			// calcula termo independente para com o valor de i
			double termo_indep = evaluator_evaluate_x(sl->fx[sl->k], (double)(i));
			// calcula o xi da próxima iteração
			resultado = termo_indep - resultado;
			resultado = resultado / sl->coeficientes[i][i];
        
			// Diferenca entre o xi (próximo) e o X[i] (atual)
			double diferenca = fabs(resultado - incognitas[i]);
			incognitas[i] = resultado;

			if (diferenca > norma)
				norma = diferenca;
		}
    }

	return incognitas;
}


// double *calcula_residuo_otimizado(double **diagonais, double *termos_indep, int n_diagonais, double *incognitas, int n_matriz){
// 	double *residuos = malloc(n_matriz * sizeof(double));
// 	int metade_diagonal = (n_diagonais - 1) / 2;

// 	for (int linha = 0; linha < n_matriz; linha++){
// 		double resultado = 0.0;

// 		// esquerda (abaixo da diagonal)
// 		for (int coluna = linha - 1; ((coluna >= 0) && (coluna >= linha - metade_diagonal)); coluna--){
// 			// na esquerda, o indice que identifica o vetor eh 
// 			int indice_diagonal = metade_diagonal + linha;

// 			// na esquerda, o indice dos valores do vetor vai diminuindo de um em um conforme se afasta da diagonal
// 			int indice_valores = coluna;

// 			if (indice_diagonal < n_diagonais)
// 				resultado += diagonais[indice_diagonal][indice_valores] * incognitas[coluna];
// 		}

// 		// direita (acima da diagonal)
// 		for (int coluna = linha + 1; coluna < metade_diagonal; coluna++){
// 			// na direita, o indice dos valores do vetor eh igual ao numero da linha
// 			int indice_valores = linha;

// 			// na direita, o indice que identifica o vetor eh igual ao numero da coluna
// 			int indice_diagonal = metade_diagonal - coluna;

// 			if (indice_diagonal < n_diagonais)
// 				resultado += diagonais[indice_diagonal][indice_valores] * incognitas[coluna];
// 		}
// 		resultado += diagonais[metade_diagonal][linha] * incognitas[linha];
// 		residuos[linha] = resultado - termos_indep[linha];
// 	}

// 	return residuos;
// }


void gauss_seidel_otimizado(double **diagonais, int n_diagonais, unsigned int max_iter, double epsilon, int n_matriz, double *termos_indep){
    double *incognitas = calloc(n_matriz, sizeof(double));
	int metade_diagonal = (n_diagonais - 1) / 2;
	double norma = 1.0 + epsilon;

    for (int iter = 0; norma > epsilon && iter < max_iter; iter++){
		norma = 0.0;
		for (int linha = 0; linha < n_matriz; linha++){
			double resultado = 0.0;

			// direita (acima da diagonal)
			printf("dir\n");
			int indice_diagonal = 0;
			for (int coluna = 0; coluna < n_matriz; coluna++){
				// na direita, o indice dos valores do vetor eh igual ao numero da linha
				int indice_valores = linha;
				int limite = n_matriz - metade_diagonal + indice_diagonal;
				int coluna_incognitas = metade_diagonal - indice_diagonal;
				// int coluna_incognitas = metade_diagonal - indice_diagonal + linha;
				if ((indice_diagonal < metade_diagonal) && (indice_valores < limite)){
					printf("DIA[%d][%d] %lf INC[%d] %lf\n", indice_diagonal, indice_valores, diagonais[indice_diagonal][indice_valores], coluna_incognitas, incognitas[coluna_incognitas]);
					resultado += diagonais[indice_diagonal][indice_valores] * incognitas[coluna_incognitas];
				}
				indice_diagonal += 1;
			}

			// esquerda (abaixo da diagonal)
			printf("esq\n");
			indice_diagonal = metade_diagonal + 1;
			for (int coluna = linha; (coluna >= 0); coluna--){
				// na esquerda, o indice dos valores do vetor vai diminuindo de um em um conforme se afasta da diagonal
				int indice_valores = coluna - 1;
				if ((indice_diagonal < n_diagonais) && (indice_valores >= 0)){
					printf("DIA[%d][%d] %lf INC[%d] %lf\n", indice_diagonal, indice_valores, diagonais[indice_diagonal][indice_valores], coluna-1, incognitas[coluna-1]);
					resultado += diagonais[indice_diagonal][indice_valores] * incognitas[coluna-1];
				}
				indice_diagonal += 1;
			}

			resultado = (termos_indep[linha] - resultado) / diagonais[metade_diagonal][linha];
			double diferenca = fabs(resultado - incognitas[linha]);
			incognitas[linha] = resultado;

			if (diferenca > norma)
				norma = diferenca;

			// for (int l = 0; l < n_matriz; l++)
			// 	printf("x: %lf\n", incognitas[l]);
			// printf("\n");
		}
		printf("\n");
	}
}


// double *gauss_seidel_sem_vetores(SisLinear_T *sl, unsigned int max_iter, double epsilon){
// 	// Contem aproximacao inicial
//     double *incognitas = calloc(sl->n, sizeof(double));
// 	// 	criterio de parada quando a norma das incognitas eh menor que epsilon 
// 	double norma = 1.0 + epsilon;

//     for (int iter = 0; norma > epsilon && iter < max_iter; iter++){
// 		norma = 0.0;
// 		for (int i = 0; i < sl->n; i++){
// 	        double resultado = 0;
// 			for (int j = 0; j < sl->n; j++)
// 				if (i != j)
// 					resultado += sl->coeficientes[i][j] * incognitas[j];
			
// 			// calcula termo independente para com o valor de i
// 			double termo_indep = evaluator_evaluate_x(sl->fx[sl->k], (double)(i));
// 			// calcula o xi da próxima iteração
// 			resultado = termo_indep - resultado;
// 			resultado = resultado / sl->coeficientes[i][i];
        
// 			// Diferenca entre o xi (próximo) e o X[i] (atual)
// 			double diferenca = fabs(resultado - incognitas[i]);
// 			incognitas[i] = resultado;

// 			if (diferenca > norma)
// 				norma = diferenca;
// 		}
//     }

// 	return incognitas;
// }



double pega_coeficiente(SisLinear_T *sl, int linha, int coluna){
	int qual_vetor = 0, coeficiente = 0; 

	qual_vetor = (sl->k-1)/2 + linha - coluna;

	// diagonal superior
	if ( linha - coluna <= 0 )
		coeficiente = linha;
	// inferior
	else if ( linha - coluna > 0 )
		coeficiente = coluna;

	// sl->coeficientes[qual_vetor][coeficiente];
	return ;
}


double *gauss_seidel_teste(SisLinear_T *sl, unsigned int max_iter, double epsilon){
	// Contem aproximacao inicial
    double *incognitas = calloc(sl->n, sizeof(double));
	// 	criterio de parada quando a norma das incognitas eh menor que epsilon 
	double norma = 1.0 + epsilon;

    for (int iter = 0; norma > epsilon && iter < max_iter; iter++){
		norma = 0.0;
		
		// metade de cima da matriz
		for (int i = 0; i <= (sl->k-1)/2; i++){
	        double resultado = 0;
			for (int j = 0; j <= (sl->k-1)/2 + i; j++)
				if (i != j)
					resultado += sl->coeficientes[i][j] * incognitas[j];
			
			// calcula termo independente para com o valor de i
			double termo_indep = evaluator_evaluate_x(sl->fx[sl->k], (double)(i));
			// calcula o xi da próxima iteração
			resultado = termo_indep - resultado;
			resultado = resultado / sl->coeficientes[i][i];
        
			// Diferenca entre o xi (próximo) e o X[i] (atual)
			double diferenca = fabs(resultado - incognitas[i]);
			incognitas[i] = resultado;

			if (diferenca > norma)
				norma = diferenca;
		}
		
		// metade de baixo da matriz
		for (int i = (sl->k-1)/2+1; i < sl->n; i++){
	        double resultado = 0;
			for (int j = 1 + i; j < sl->n; j++)
				if (i != j)
					resultado += sl->coeficientes[i][j] * incognitas[j];
			
			// calcula termo independente para com o valor de i
			double termo_indep = evaluator_evaluate_x(sl->fx[sl->k], (double)(i));
			// calcula o xi da próxima iteração
			resultado = termo_indep - resultado;
			resultado = resultado / sl->coeficientes[i][i];
        
			// Diferenca entre o xi (próximo) e o X[i] (atual)
			double diferenca = fabs(resultado - incognitas[i]);
			incognitas[i] = resultado;

			if (diferenca > norma)
				norma = diferenca;
		}
    }

	return incognitas;
}