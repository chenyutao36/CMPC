#include <stdlib.h>
#include <stdio.h>
#include "string.h"

#include "common.h"
#include "full_condensing.h"

extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern void dsymm_(char*, char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern void dsymv_(char*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern void dsyrk_(char*, char*, int*, int*, double*, double*, int*, double*, double*, int*);

full_condensing_workspace* full_condensing_workspace_create(model_size *size)
{
    int nx=size->nx;
    int nu=size->nu;
    // int ny=size->ny;
    // int nyN=size->nyN;
    // int np=size->np;
    int nbg = size->nbg;
    int nbgN = size->nbgN;
    int N = size->N;

    full_condensing_workspace *work = (full_condensing_workspace*)malloc(sizeof(full_condensing_workspace));

    work->G = (double*)malloc(nx*N*N*nu * sizeof(double));
    work->W = (double *)malloc(nx*N*N*nu * sizeof(double));       
    work->w = (double *)malloc(nx*N * sizeof(double));             
    work->L = (double *)malloc((N+1)*nx * sizeof(double));         
    work->Hi = (double *)malloc(nu*nu * sizeof(double));        
    work->Cci = (double *)malloc(nbg*nu * sizeof(double));   
    work->CcN = (double *)malloc(nbgN*nu * sizeof(double)); 

    work->Hc = (double*)malloc(nu*N*N*nu * sizeof(double));
    work->gc = (double*)malloc(nu*N * sizeof(double));
    work->Cc = (double*)malloc(N*nu*(N*nbg+nbgN) * sizeof(double));
    work->lcc = (double*)malloc((N*nbg+nbgN) * sizeof(double));
    work->ucc = (double*)malloc((N*nbg+nbgN) * sizeof(double));

    return work;
}

int full_condensing(double *dx0, model_size *size, qp_problem *qp, full_condensing_workspace *work)
{
    int nx=size->nx;
    int nu=size->nu;
    // int ny=size->ny;
    // int nyN=size->nyN;
    // int np=size->np;
    int nbg = size->nbg;
    int nbgN = size->nbgN;
    int N = size->N;

    double *Q = qp->Q;
    double *S = qp->S;
    double *R = qp->R;
    double *A = qp->A;
    double *B = qp->B;
    double *C = qp->C;
    double *D = qp->D;
    double *CN = qp->CN;
    double *gx = qp->gx;
    double *gu = qp->gu;   
    double *b = qp->b;
    // double *lb_u = qp->lb_u;
    // double *ub_u = qp->ub_u;
    double *lb_g = qp->lb_g;
    double *ub_g = qp->ub_g;

    double *G = work->G;    
    double *W = work->W;
    double *w = work->w;
    double *L = work->L;
    double *Hi = work->Hi;
    double *Cci = work->Cci;
    double *CcN = work->CcN;
    double *Hc = work->Hc;
    double *gc = work->gc;
    double *Cc = work->Cc;
    double *lcc = work->lcc;
    double *ucc = work->ucc;     

    char *nTrans = "N", *Trans="T", *SIDE = "L", *UPLO = "L";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    int one_i = 1; 

    int i,j;
    /* compute G */
    for(i=0;i<N;i++){
        memcpy(G+(i*N+i)*nx*nu, B+i*nx*nu, nx*nu*sizeof(double));
        for (j=i+1;j<N;j++){
            dgemm_(nTrans, nTrans, &nx, &nu, &nx, &one_d, A+j*nx*nx, &nx, G+(i*N+j-1)*nx*nu, &nx, &zero, G+(i*N+j)*nx*nu, &nx);
        }
    }

    /* Compute Hc */
    for(i=0;i<N;i++){
        dsymm_(SIDE, UPLO, &nx, &nu, &one_d, Q+N*nx*nx, &nx, G+(i*N+N-1)*nx*nu, &nx, &zero, W+(i*N+N-1)*nx*nu, &nx);
        for(j=N-1;j>i;j--){        
            dgemm_(Trans, nTrans, &nu, &nu, &nx, &one_d, S+j*nx*nu, &nx, G+(i*N+j-1)*nx*nu, &nx, &zero, Hi, &nu);                     
            dgemm_(Trans, nTrans, &nu, &nu, &nx, &one_d, B+j*nx*nu, &nx, W+(i*N+j)*nx*nu, &nx, &one_d, Hi, &nu);
            Block_Fill(nu, nu, Hi, Hc, j*nu, i*nu, N*nu);
            Block_Fill_Trans(nu, nu, Hi, Hc, i*nu, j*nu, N*nu);
            dgemm_(Trans, nTrans, &nx, &nu, &nx, &one_d, A+j*nx*nx, &nx, W+(i*N+j)*nx*nu, &nx, &zero, W+(i*N+j-1)*nx*nu, &nx); 
            dsymm_(SIDE, UPLO, &nx, &nu, &one_d, Q+j*nx*nx, &nx, G+(i*N+j-1)*nx*nu, &nx, &one_d, W+(i*N+j-1)*nx*nu, &nx);
        }
        memcpy(Hi,R+i*nu*nu,nu*nu*sizeof(double));
        dgemm_(Trans, nTrans, &nu, &nu, &nx, &one_d, B+i*nx*nu, &nx, W+(i*N+i)*nx*nu, &nx, &one_d, Hi, &nu);
        Block_Fill(nu, nu, Hi, Hc, i*nu, i*nu, N*nu);      
    }

    /* Compute Cc */
    if (nbg>0){         
        for(i=0;i<N;i++){             
            Block_Fill_Trans(nbg,nu, D+i*nbg*nu, Cc, i*nu, i*nbg, N*nu);
        for(j=i+1;j<N;j++){   
            dgemm_(nTrans, nTrans, &nbg, &nu, &nx, &one_d, C+j*nbg*nx, &nbg, G+(i*N+j-1)*nx*nu, &nx, &zero, Cci, &nbg);                   
            Block_Fill_Trans(nbg, nu, Cci, Cc, i*nu, j*nbg, N*nu);
        }    
        }  
    }

    /* Compute CcN */
    if (nbgN>0){          
        for(i=0;i<N;i++){                 
            dgemm_(nTrans, nTrans, &nbgN, &nu, &nx, &one_d, CN, &nbgN, G+(i*N+N-1)*nx*nu, &nx, &zero, CcN, &nbgN);
            Block_Fill_Trans(nbgN, nu, CcN, Cc, i*nu, N*nbg, N*nu);
        }
    }

    /* compute L */
    memcpy(L,dx0, nx*sizeof(double)); 
    for(i=0;i<N;i++){
        memcpy(L+(i+1)*nx, b+i*nx, nx*sizeof(double)); 
        dgemv_(nTrans,&nx,&nx,&one_d,A+i*nx*nx,&nx,L+i*nx,&one_i,&one_d,L+(i+1)*nx,&one_i);
    }

    /* compute gc */
    memcpy(w+(N-1)*nx,gx+N*nx,nx*sizeof(double));
    dsymv_(UPLO, &nx, &one_d, Q+N*nx*nx, &nx, L+N*nx, &one_i, &one_d, w+(N-1)*nx, &one_i);
    for(i=N-1;i>0;i--){
        memcpy(gc+i*nu,gu+i*nu,nu*sizeof(double));
        dgemv_(Trans,&nx,&nu,&one_d,S+i*nx*nu,&nx,L+i*nx,&one_i,&one_d,gc+i*nu,&one_i);
        dgemv_(Trans,&nx,&nu,&one_d,B+i*nx*nu,&nx,w+i*nx,&one_i,&one_d,gc+i*nu,&one_i);
         
        memcpy(w+(i-1)*nx, gx+i*nx, nx*sizeof(double));
        dsymv_(UPLO, &nx, &one_d, Q+i*nx*nx, &nx, L+i*nx, &one_i, &one_d, w+(i-1)*nx, &one_i);
        dgemv_(Trans,&nx,&nx,&one_d,A+i*nx*nx,&nx,w+i*nx,&one_i,&one_d,w+(i-1)*nx,&one_i);
    }   
    memcpy(gc,gu,nu*sizeof(double));
    dgemv_(Trans,&nx,&nu,&one_d,S,&nx,L,&one_i,&one_d,gc,&one_i);
    dgemv_(Trans,&nx,&nu,&one_d,B,&nx,w,&one_i,&one_d,gc,&one_i);

    /* Compute cc */
    if (nbg>0){                    
        for(i=0;i<N;i++){
            memcpy(lcc+i*nbg,lb_g+i*nbg,nbg*sizeof(double));
            dgemv_(nTrans,&nbg,&nx,&minus_one_d,C+i*nbg*nx,&nbg,L+i*nx,&one_i,&one_d,lcc+i*nbg,&one_i); 

            memcpy(ucc+i*nbg,ub_g+i*nbg,nbg*sizeof(double));
            dgemv_(nTrans,&nbg,&nx,&minus_one_d,C+i*nbg*nx,&nbg,L+i*nx,&one_i,&one_d,ucc+i*nbg,&one_i);
        }        
    }   
    
    /* Compute ccN */
    if (nbgN>0){          
        memcpy(lcc+N*nbg,lb_g+N*nbg,nbgN*sizeof(double));
        dgemv_(nTrans,&nbgN,&nx,&minus_one_d,CN,&nbgN,L+N*nx,&one_i,&one_d,lcc+N*nbg,&one_i);
        memcpy(ucc+N*nbg,ub_g+N*nbg,nbgN*sizeof(double));
        dgemv_(nTrans,&nbgN,&nx,&minus_one_d,CN,&nbgN,L+N*nx,&one_i,&one_d,ucc+N*nbg,&one_i);
    }

    return 0;
}

void full_condensing_workspace_free(full_condensing_workspace* work)
{
    free(work->G);
    free(work->W);
    free(work->w);
    free(work->L);
    free(work->Hi);
    free(work->Cci);
    free(work->CcN);
    free(work->Hc);
    free(work->gc);
    free(work->Cc);
    free(work->lcc);
    free(work->ucc);

    free(work);
}