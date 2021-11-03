
#include "utils.h"



// 
void imprime_solucao(SisLinear_T *sl){
	printf("solucao: ");
	
	// for (int i = 0; i < sl->n; i++){
	// 	printf("%lf ", );
	// }

	printf("\n");
}


// 
void imprime_termos_indep(SisLinear_T *sl){
	for (int i = 0; i < sl->n; i++)
		printf("%1.0f ", sl->termos_indep[i]);
	printf("\n");
}


// 
void gera_termos_indep(SisLinear_T *sl){
	for (int i = 0; i < sl->n; i++)
		sl->termos_indep[i] = evaluator_evaluate_x(sl->fx[sl->k], (double)(i));
}


// 
void imprime_sistema_linear(SisLinear_T *sl){
	for (int i = 0; i < sl->n; i++){
		for (int j = 0; j < sl->n; j++)
			printf("%1.0f ", sl->coeficientes[i][j]);
		printf("\n");
	}
	printf("\n");
}


// 
void gera_sistema_linear(SisLinear_T *sl){
	int w = 0, l = sl->k - 1;

	for (int i = 0; i < sl->n; i++)
		sl->coeficientes[i][i] = evaluator_evaluate_x(sl->fx[l/2], (double)i);


	for (int indice = (sl->k - 1)/2; indice >= 1; indice--){
		for (int col = indice, lin = 0; col < sl->n; col++, lin++){
			sl->coeficientes[lin][col] = evaluator_evaluate_x(sl->fx[w], (double)lin);
			sl->coeficientes[col][lin] = evaluator_evaluate_x(sl->fx[l], (double)lin);
		}
		w++;
		l--;
	}
}


// 
void ajusta_sistema_linear(SisLinear_T *sl){
	// aloca um vetor de n ponteiros para linhas
	sl->coeficientes = malloc(sl->n * sizeof (double *));
	
	// aloca um vetor com todos os elementos da matriz
	sl->coeficientes[0] = sl->matriz = calloc(sl->n * sl->n, sizeof(double));
	
	// ajusta os demais ponteiros de linhas (i > 0)
	for (int i = 1; i < sl->n; i++)
		sl->coeficientes[i] = sl->coeficientes[0] + i * sl->n;
}


// destroi o sl
void destroi_sistema_linear(SisLinear_T *sl){
	for (int i = 0; i <= sl->k; i++)
		free(sl->fx[i]);

	free(sl->matriz);
	free(sl->coeficientes);
	free(sl->termos_indep);
	free(sl);
}


// aloca uma regiao de memoria para o sl
SisLinear_T *aloca_sistema_linear(){
	return (SisLinear_T *)malloc(sizeof(SisLinear_T));
}


// le a os dados de entrada
void read_input(SisLinear_T *sl, double *epsilon, unsigned int *max_iter){
	// dimensao do SL
	scanf(" %u", &sl->n);
	
	// número de diagonais (sempre um número ímpar, k < n)
    scanf(" %u", &sl->k);

	sl->termos_indep = (double *)malloc(sl->n * sizeof(double));
    
	// string da expressão da função fi(x), sem espaços, na variável x.
	char fx[SIZE_BUFFER + 1];

	sl->fx = (int **)malloc((sl->k + 1) * sizeof(void *));
	
	for (int i = 0; i <= sl->k; i++){
		scanf(" %s", fx);
		sl->fx[i] = evaluator_create(fx);
	}

	// tolerância no método de Gauss-Seidel
    scanf(" %lf", epsilon);
	
	// número máximo de iterações
    scanf(" %u", max_iter);
}


// Retorna tempo em milisegundos
double timestamp(void){
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}