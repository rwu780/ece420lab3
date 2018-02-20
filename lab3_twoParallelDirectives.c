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

	for(int k = 0; k< N-1; k++){
		// eliminate elements below the diagonal to zero one column after another
		// Pivoting
		int kp = 0;
		double maxAbsolute = 0;

		# pragma omp parallel num_threads(p)
		{
			# pragma omp for
			for(int i = k; i<N; i++){
				#pragma omp critical
				{
					if(fabs(G[i][k]) >= maxAbsolute){
						maxAbsolute = fabs(G[i][k]);
						kp = i;
					}
				}
			}
		
			//Swap row k and row kp
			# pragma omp critical
			{
				double *temp = G[k];
				G[k] = G[kp];
				G[kp] = temp;
			}
			/*Elimination*/
	
			//OPEN MP HERE
			# pragma omp barrier
			# pragma omp for
			for(int i = k+1; i< N; i++){
				double temp = G[i][k] / G[k][k];
				for(int j = k; j<N+1; j++){
					G[i][j] = G[i][j] - temp * G[k][j]; /* Row replacement */
				}
			}
		}
	}

	# pragma omp parallel num_threads(p)
	{

		for(int k = N-1; k > 0 ; k--){
			# pragma omp for
		//Eliminate elements to zero for each column one after another
			for(int i = 0; i < k; i++){
				//Row replacement one row after another
				//Open MP HERE
				G[i][N] = G[i][N] - G[i][k]/G[k][k] * G[k][N];
				G[i][k] = 0;
	
			}
		}
	


		# pragma omp for
		for(int i = 0; i < N; i++){
			//open MP here
			X[i] = G[i][N] / G[i][i];
		}	
	}
	GET_TIME(end);

	double duration = end - start;
	Lab3SaveOutput(X, N, duration);

	printf("%f\n", duration);

	free(X);
	free(G);
	return 0;
}

