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

int lucol(int n, double A[][nmax], int p[]) {
	int i, j, k, imax;
	for (k = 0; k < n - 1; k ++) {
		imax = k;
		for (i = k + 1; i < n; i ++)
			if (fabs(A[i][k]) > fabs(A[imax][k]))
				imax = i;
		p[k] = imax;
		if (imax != k)
			for (j = 0; j < n; j ++)
				swap(&A[k][j], &A[imax][j]);
		if (fabs(A[k][k]) < E)
			return -1;
		for (i = k + 1; i < n; i ++)
			A[i][k] = A[i][k]/A[k][k];
		for (j = k + 1; j < n; j ++)
			for (i = k + 1; i < n; i ++)
				A[i][j] -= A[k][j]*A[i][k];
	}
	if (fabs(A[n-1][n-1]) < E)
		return -1;
	return 0;
}

int sscol(int n, double A[][nmax], int p[], double b[]) {
	int i, j;
	for (i = 0; i < n - 1; i ++)
		swap(&b[i], &b[p[i]]);

	for (j = 0; j < n; j ++) {
		for (i = j + 1; i < n; i ++)
			b[i] -= A[i][j]*b[j];
	}

	for (j = n - 1; j >= 0; j --) {
		b[j] = b[j]/A[j][j];
		for (i = j - 1; i >= 0; i --)
			b[i] -= A[i][j]*b[j];
	}
	return 0;
}

int main() {
	char file_name[100];
	FILE *file;
	double A[nmax][nmax], b[nmax], duration;
	int n, i, j, k, p[nmax];
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
	if (lucol(n, A, p) == -1) {
		printf("Matriz é singular.\n");
		return -1;
	}
	end = clock();
	duration = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Lucol tempo %e segundos\n", duration);
	start = clock();
	sscol(n, A, p, b);
	end = clock();
	duration = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Sscol tempo %e segundos\n", duration);
	for (i = 0; i < n; i ++) {
		if (b[i] - (1 + i%(n/100)) > E || b[i] - (1 + i%(n/100)) < -E)
			printf("Erro! %e  %d %d\n", b[i],-(1 + i%(n/100)), i);
	}
	printf("Fim da Análise!\n");
	return 0;
}
