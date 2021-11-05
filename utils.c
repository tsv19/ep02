
#include "utils.h"


// Imprime o tempo que o metodo de gauss-seidel demorou para solucionar o 
// sistema linear
void imprime_tempo(double tempo){
	printf("tempo: %lf\n", tempo);
}

// Imprime as k diagonais do sistema linear
void imprime_diagonais(SisLinear_T *sl){
	int fx = 0;
	double value = 0.0;

	// Imprime as diagonais superiores da mais externa ate a
	// diagonal principal
	for (int indice = (sl->k - 1)/2; indice >= 0; indice--){
		for (int k = indice, i = 0; k < sl->n; k++, i++){
			value = evaluator_evaluate_x(sl->fx[fx], (double)i);
			printf("%1.0f ", value);
		}
		printf("\n");
		fx++;
	}

	// Imprime as diagonais inferiores da mais interna para a
	// mais externa  
	for (int indice = 1; indice <= (sl->k - 1)/2; indice++){
		for (int k = indice, i = 0; k < sl->n; k++, i++){
			value = evaluator_evaluate_x(sl->fx[fx], (double)i);
			printf("%1.0f ", value);
		}
		printf("\n");
		fx++;
	}
}


// Imprime a solucao encontrada por gauss-seidel
void imprime_solucao(double *vet, int n){
	printf("----\nsolucao: ");
	for (int i = 0; i < n; i++)
		printf("%lf ", vet[i]);
	printf("\n");
}


// Imprime os termos independentes gerados para o sistema linear
void imprime_termos_indep(SisLinear_T *sl){
	for (int i = 0; i < sl->n; i++)
		printf("%1.0f ", evaluator_evaluate_x(sl->fx[sl->k], (double)(i)));
	printf("\n");
}


// apagar(?)
void gera_termos_indep(SisLinear_T *sl){
	for (int i = 0; i < sl->n; i++)
		sl->termos_indep[i] = evaluator_evaluate_x(sl->fx[sl->k], (double)(i));
}


// Gera os valores das k diagonais
void gera_sistema_linear(SisLinear_T *sl){
	int fx_sup = 0, fx_inf = sl->k - 1;

	// gera os valores para a diagonal principal
	for (int i = 0; i < sl->n; i++)
		sl->coeficientes[i][i] = evaluator_evaluate_x(sl->fx[fx_inf/2], (double)i);

	// gera os valores para as diagonais superiores e inferiores simultaneamente
	// das mais distantes da diagonal principal ate as mais proximas
	for (int indice = (sl->k - 1)/2; indice >= 1; indice--){
		for (int col = indice, lin = 0; col < sl->n; col++, lin++){
			sl->coeficientes[lin][col] = evaluator_evaluate_x(sl->fx[fx_sup], (double)lin);
			sl->coeficientes[col][lin] = evaluator_evaluate_x(sl->fx[fx_inf], (double)lin);
		}
		fx_sup++;
		fx_inf--;
	}
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
SisLinear_T *aloca_sistema_linear(unsigned int n, unsigned int k, int **fx){
	SisLinear_T *sl = (SisLinear_T *)malloc(sizeof(SisLinear_T));
	sl->n = n;
	sl->k = k;
	sl->termos_indep = (double *)malloc(sl->n * sizeof(double));
	sl->fx = fx;

	// aloca um vetor de n ponteiros para linhas
	sl->coeficientes = malloc(sl->n * sizeof (double *));
	
	// aloca um vetor com todos os elementos da matriz
	sl->coeficientes[0] = sl->matriz = calloc(sl->n * sl->n, sizeof(double));
	
	// ajusta os demais ponteiros de linhas (i > 0)
	for (int i = 1; i < sl->n; i++)
		sl->coeficientes[i] = sl->coeficientes[0] + i * sl->n;
	
	return sl;
}


// le a os dados de entrada
SisLinear_T *read_input(double *epsilon, unsigned int *max_iter){
	unsigned int n, k;
	int **temp;

	// dimensao do SL
	scanf(" %u", &n);
	// número de diagonais (sempre um número ímpar, k < n)
    scanf(" %u", &k);
	// string da expressão da função fi(x), sem espaços, na variável x.
	char fx[SIZE_BUFFER + 1];

	temp = (int **)malloc((k + 1) * sizeof(void *));
	
	for (int i = 0; i <= k; i++){
		scanf(" %s", fx);
		temp[i] = evaluator_create(fx);
	}

	// tolerância no método de Gauss-Seidel
    scanf(" %lf", epsilon);    
	// número máximo de iterações
    scanf(" %u", max_iter);

	return aloca_sistema_linear(n, k, temp);
}


// Retorna tempo em milisegundos
double timestamp(void){
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}