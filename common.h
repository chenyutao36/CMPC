#ifndef COMMON_H_
#define COMMON_H_

#include <sys/time.h>

typedef struct{
    int nx;
    int nu;
    int ny;
    int nyN;
    int np;
    int nbx;
    int nbu;
    int nbg;
    int nbgN;
    int N;
    int N2;

    int *nbx_idx;
    int *nbu_idx;
}model_size;

typedef struct{
    int qpsolver; // 0:qore, 1:hpipm_ocp, 2:hpipm_pcond
    int shifting; // 0:no, 1: yes
}rti_opt;

typedef struct{
    struct timeval t1;
    struct timeval t2;
}Timer;

void Block_Fill(int m, int n, double *Gi, double *G,
     int idm, int idn, int ldG);

void Block_Fill_Trans(int m, int n, double *Gi, double *G,
     int idm, int idn, int ldG);

void set_zeros(int dim, double *A);

void print_matrix(double *A, int m, int n);

void print_vector(double *x, int m);
void print_vector_trans(double *x, int m);

void print_vector_int(int *x, int m);
void print_vector_int_trans(int *x, int m);

void regularization(int n, double *A, double reg);

void dgemv_n_3l(int m, int n, double alpha, double *A, int lda, double *x, double beta, double *y, double *z);

void Tic(Timer *timer);

double Toc(Timer *timer);


#endif