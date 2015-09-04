#include <math.h>
#include <stdio.h>
#include <time.h>
#define nmax 700
#define E 0.0001 /*erro*/

void swap(double* a, double* b) {
	double temp;
	temp = *a;
	*a = *b;
	*b = temp;
}

int cholrow(int n, double A[][nmax]) {
	int i, j, k;
	for (j = 0; j < n; j ++) {
		for (k = 0; k < j; k ++)
			A[j][j] -= A[j][k]*A[j][k];
		if (A[j][j] < E) /* <= 0 */
			return -1;
		A[j][j] = sqrt(A[j][j]);
		for (i = j + 1; i < n; i ++) {
			for (k = 0; k < j; k ++) {
				A[i][j] -= A[i][k]*A[j][k];
			}
			A[i][j] = A[i][j]/A[j][j];
		}
	}
	return 0;
}

int forwrow(int n, double A[][nmax], double b[]) {
	int i, j;
	for (i = 0; i < n; i ++) {
		if (fabs(A[i][i]) < E) /* == 0 */
			return -1;
		for (j = 0; j < i; j ++)
			b[i] -= A[i][j]*b[j];
		b[i] = b[i]/A[i][i];
	}
	return 0;
}

int backrow(int n, double A[][nmax], double b[], int trans) {
	int i, j;
	if (trans) { /* trans == 1 */
		for (j = n - 1; j >= 0; j --) {
			if (fabs(A[j][j]) < E) /* == 0 */
				return -1;
			b[j] = b[j]/A[j][j];
			for (i = 0; i < j; i ++)
				b[i] -= A[j][i]*b[j];
		}
	}
	else { /* trans == 0 */
		for (i = n - 1; i >= 0; i --) {
			if (fabs(A[i][i]) < E) /* == 0 */
				return -1;
			for (j = i + 1; j < n; j ++)
				b[i] -= A[i][j]*b[j];
			b[i] = b[i]/A[i][i];
		}
	}
	return 0;
}

int main() {
	char file_name[100];
	FILE *file;
	double A[nmax][nmax], b[nmax], duration;
	int n, i, j, k;
	clock_t start, end;

	printf("Nome do Arquivo: ");
	scanf("%s", file_name);
	file = fopen(file_name, "r");
	if (file == NULL) {
		fprintf(stderr, "Não foi possível abrir o arquivo!\n");
		return -1;
	}

	fscanf(file, "%d", &n);
	for (k = 0; k < n*n; k ++) {
		fscanf(file, "%d %d", &i, &j);
		fscanf(file, "%lf", &A[i][j]);
	}
	for (k= 0; k < n; k ++) {
		fscanf(file, "%d", &i);
		fscanf(file, "%lf", &b[i]);
	}

	start = clock();
	if (cholrow(n, A) == -1) {
		printf("Matriz não é definida positiva.\n");
		return -1;
	}
	end = clock();
	duration = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Choleskyrow tempo %e segundos\n", duration);
	start = clock();
	if (forwrow(n, A, b) == -1) {
		printf("Matriz é singular.\n");
		return -1;
	}
	end = clock();
	duration = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Forwrow tempo %e segundos\n", duration);
	/*for (i = 0; i < n; i++)
		for (j = 0; j < i; j++)
			swap(&A[i][j], &A[j][i]);*/
	start = clock();
	if (backrow(n, A, b, 1) == -1) {
		printf("Matriz é singular.\n");
		return -1;
	}
	end = clock();
	duration = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Backrow tempo %e segundos\n", duration);
	for (i = 0; i < n; i ++) {
		if (fabs(b[i] - (1 + i%(n/100))) > E)
			printf("Erro! %e  %d %d\n", b[i],-(1 + i%(n/100)), i);
	}
	printf("Fim da Análise!\n");
	return 0;
}
