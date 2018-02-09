#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "Lab3IO.h"
#include "timer.h"

double *X;
double **G;
int N, p;

void swapRow(int row1, int row2);
int findAbsoluteMax(double **U, int k, int n);
void Gaussian(void * rank);
void Jordan(void * rank);


int main(int argc, char const *argv[])
{
	if(argc != 2){
		printf("Please pass in number of threads\n");
		exit(1);
	}

	p = atoi(argv[1]);

	double start, end;

	Lab3LoadInput(&G, &N);
	
	GET_TIME(start);
	X = calloc(N, sizeof(double));

	Gaussian(NULL);
	Jordan(NULL);

	
	for(int i = 0; i < N; i++){
		//open MP here
		X[i] = G[i][N] / G[i][i];
	}
	

	GET_TIME(end);

	double duration = end - start;
	Lab3SaveOutput(X, N, duration);

	printf("%f\n", duration);

	free(X);
	free(G);
	return 0;
}

void Jordan(void * rank){
	for(int k = N-1; k > 0 ; k--){
		//Eliminate elements to zero for each column one after another
		for(int i = 0; i < k; i++){
			//Row replacement one row after another

			//Open MP HERE
			G[i][N] = G[i][N] - G[i][k]/G[k][k] * G[k][N];
			G[i][k] = 0;
		}
	}
}

void Gaussian(void * rank){
	for(int k = 0; k< N-1; k++){
		// eliminate elements below the diagonal to zero one column after another
		// Pivoting
		// In U, from row k to row n, find the row kp, that has the maximum absolute value
		int kp = findAbsoluteMax(G, k, N);

		//Swap row k and row kp
		swapRow(k, kp);

		/*Elimination*/

		//OPEN MP HERE

		for(int i = k+1; i< N; i++){
			double temp = G[i][k] / G[k][k];
			for(int j = k; j<N+1; j++){
				G[i][j] = G[i][j] - temp * G[k][j]; /* Row replacement */
			}
		}
	}
}

void swapRow(int row1, int row2){

	//OPEN MP HERE
# pragma omp parallel num_threads(p){
		for(int i = 0; i< N+1; i++){
			double temp = G[row1][i];
			G[row1][i] = G[row2][i];
			G[row2][i] = temp;
		}
	}
}

int findAbsoluteMax(double **U, int k, int n){

	//OPEN MP HERE
	int kp = 0;
	int maxAbsolute = 0;
	for(int i = k; i<n; i++){
		if(fabs(U[i][k]) >= maxAbsolute){
			maxAbsolute = fabs(U[i][k]);
			kp = i;
		}
	}
	return kp;
}
