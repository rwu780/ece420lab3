#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

	Lab3LoadInput(&G, &N);
	X = calloc(N, sizeof(double));

	Gaussian(NULL);
	Jordan(NULL);


	for(int i = 0; i < N; i++){
		X[i] = G[i][i];
	}

	Lab3SaveOutput(X, N, 0.0);

	/* code */
	return 0;
}

void Jordan(void * rank){
	for(int k = N-1; k > 2 ; k--){
		//Eliminate elements to zero for each column one after another
		for(int i = 0; i < k-1; i++){
			//Row replacement one row after another
			G[i][N+1] = G[i][N+1] - G[i][k]/G[k][k] * G[k][N+1];
			G[i][k] = 0;
		}
	}
}

void Gaussian(void * rank){
	// double **U = G;
	for(int k = 0; k< N-1; k++){
		/* eliminate elements below the diagonal to zero one column after another */
		//Pivoting
		// In U, from row k to row n, find the row kp, that has the maximum absolute value
		int kp = findAbsoluteMax(G, k, N);

		//Swap row k and row kp
		swapRow(k, kp);

		/*Elimination*/
		for(int i = k+1; i< N; i++){
			double temp = G[i][k] / G[k][k];
			for(int j = k; j<N+1; j++){
				G[i][j] = G[i][j] - temp * G[k][j]; /* Row replacement */
			}
		}
	}
}

void swapRow(int row1, int row2){

	double temp[N+1];

	for(int i = 0; i<N+1; i++){
		temp[i] = G[row1][i];
		G[row1][i] = G[row2][i];
		G[row2][i] = temp[i];
	}
}

int findAbsoluteMax(double **U, int k, int n){
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
