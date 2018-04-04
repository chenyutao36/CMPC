#include <stdio.h>
#include "common.h"
#include <sys/time.h>

void Block_Fill(int m, int n, double *Gi, double *G, int idm, int idn, int ldG){
       
    int i,j;
    int s;
    for (j=0;j<n;j++){
        s = idn*ldG + idm + j*ldG;
        for (i=0;i<m;i++){
            G[s+i] = Gi[j*m+i];
        }
    }
       
}

void Block_Fill_Trans(int m, int n, double *Gi, double *G, int idm, int idn, int ldG){
       
    int i,j;
    int s;
    for (j=0;j<m;j++){
        s = idn*ldG + idm + j*ldG;
        for (i=0;i<n;i++){
            G[s+i] = Gi[i*m+j];
        }
    }
       
}

void set_zeros(int dim, double *A){
    int i;
    for (i=0;i<dim;i++)
        A[i] = 0.0;
}

void print_matrix(double *A, int m, int n){
    int i,j;
    for(i=0;i<m;i++){
        for(j=0;j<n;j++)
            printf("%6.4f ", A[j*m+i]);
        printf("\n");
    }
}

void print_vector(double *x, int m){
    int i;
    for(i=0;i<m;i++){
        printf("%6.4e ", x[i]);
        printf("\n");
    }
}

void print_vector_trans(double *x, int m){
    int i;
    for(i=0;i<m;i++){
        printf("%6.4e ", x[i]);      
    }
    printf("\n");
}

void print_vector_int(int *x, int m){
    int i;
    for(i=0;i<m;i++){
        printf("%4d ", x[i]);
        printf("\n");
    }
}

void print_vector_int_trans(int *x, int m){
    int i;
    for(i=0;i<m;i++){
        printf("%4d ", x[i]);
    }
    printf("\n");
}

void regularization(int n, double *A, double reg){
    int i;
    for (i=0;i<n;i++)
        if (A[i*n+i]<reg)
            A[i*n+i] = reg;
}

void dgemv_n_3l(int m, int n, double alpha, double *A, int lda, double *x, double beta, double *y, double *z)
	{

	int ii, jj;

	double tmp;

	for(ii=0; ii<m; ii++)
		z[ii] = beta * y[ii];

	for(jj=0; jj<n; jj++)
		{
		tmp = alpha * x[jj];
		for(ii=0; ii<m; ii++)
			{
			z[ii] += A[ii+lda*jj] * tmp;
			}
		}
	
	}

void Tic(Timer *timer)
{
    gettimeofday(&timer->t1, 0);
}

double Toc(Timer *timer)
{
    gettimeofday(&timer->t2, 0);

    double elapsedTime;

    elapsedTime = (timer->t2.tv_sec - timer->t1.tv_sec)*1000;      // sec to ms
    elapsedTime += (timer->t2.tv_usec - timer->t1.tv_usec) / 1000.0;   // us to ms

    return elapsedTime;
}