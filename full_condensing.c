#include <stdlib.h>
#include <stdio.h>
#include "string.h"

#include "common.h"
#include "full_condensing.h"

#include <blasfeo_target.h>
#include "blasfeo_common.h"
#include "blasfeo_i_aux_ext_dep.h"
#include "blasfeo_d_aux_ext_dep.h"
#include "blasfeo_v_aux_ext_dep.h"
#include "blasfeo_d_aux.h"
#include "blasfeo_d_kernel.h"
#include "blasfeo_d_blas.h"

extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern void dsymm_(char*, char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern void dsymv_(char*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern void dsyrk_(char*, char*, int*, int*, double*, double*, int*, double*, double*, int*);

full_condensing_workspace* full_condensing_workspace_create(model_size *size)
{
    int nx=size->nx;
    int nu=size->nu;
    int nbg = size->nbg;
    int nbgN = size->nbgN;
    int N = size->N;

    full_condensing_workspace *work = (full_condensing_workspace*)malloc(sizeof(full_condensing_workspace));

    work->sStQ = malloc(sizeof(struct blasfeo_dmat));
    work->sRt = malloc(sizeof(struct blasfeo_dmat));
    work->sBtAt = malloc(sizeof(struct blasfeo_dmat));
    work->sCt = malloc(sizeof(struct blasfeo_dmat));
    blasfeo_allocate_dmat(nx+nu, (N+1)*nx, work->sStQ);
    blasfeo_allocate_dmat(nu, N*nu, work->sRt);
    blasfeo_allocate_dmat(nx+nu, N*nx, work->sBtAt);
    blasfeo_allocate_dmat(nx, N*nbg+nbgN, work->sCt);

    work->srq = malloc(sizeof(struct blasfeo_dvec)); 
    blasfeo_allocate_dvec((N+1)*(nx+nu), work->srq);

    work->sGt = malloc(sizeof(struct blasfeo_dmat));
    work->sHtWt = malloc(sizeof(struct blasfeo_dmat));
    work->sWttmp = malloc(sizeof(struct blasfeo_dmat));
    work->shw = malloc(sizeof(struct blasfeo_dvec));
    work->sL = malloc(sizeof(struct blasfeo_dvec));
    blasfeo_allocate_dmat(N*nu, N*nx, work->sGt);
    blasfeo_allocate_dmat(nu, nx+nu, work->sHtWt);
    blasfeo_allocate_dmat(nu, nx, work->sWttmp);
    blasfeo_allocate_dvec(N*nx, work->shw);
    blasfeo_allocate_dvec((N+1)*nx, work->sL);

    work->sHc = malloc(sizeof(struct blasfeo_dmat));
    work->sCct = malloc(sizeof(struct blasfeo_dmat));
    work->sgc = malloc(sizeof(struct blasfeo_dvec));
    blasfeo_allocate_dmat(N*nu, N*nu, work->sHc);
    blasfeo_allocate_dmat(N*nu, N*nbg+nbgN, work->sCct);
    blasfeo_allocate_dvec(N*nu, work->sgc);
    
    work->G = (double*)calloc(nx*N*N*nu , sizeof(double));
    work->W = (double *)calloc(nx*N*N*nu , sizeof(double));       
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
    double *lb_g = qp->lb_g;
    double *ub_g = qp->ub_g;

    struct blasfeo_dmat *sStQ = work->sStQ;
    struct blasfeo_dmat *sRt = work->sRt;
    struct blasfeo_dmat *sBtAt = work->sBtAt;
    struct blasfeo_dmat *sCt = work->sCt;
    struct blasfeo_dmat *sGt = work->sGt;
    struct blasfeo_dmat *sHtWt = work->sHtWt;
    struct blasfeo_dmat *sWttmp = work->sWttmp;
    struct blasfeo_dmat *sHc = work->sHc;
    // struct blasfeo_dmat *sCct = work->sCct;

    // struct blasfeo_dvec *srq = work->srq;
    // struct blasfeo_dvec *shw = work->shw;
    // struct blasfeo_dvec *sL = work->sL;

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

    /* converte col maj to blasfeo struct */
    for(i=0;i<N;i++){
        blasfeo_pack_dmat(nu, nx, S+i*nx*nu, nu, sStQ, 0, i*nx);
        blasfeo_pack_dmat(nx, nx, Q+i*nx*nx, nx, sStQ, nu, i*nx);            
        blasfeo_pack_tran_dmat(nu, nu, R+i*nu*nu, nu, sRt, 0, i*nu);
        blasfeo_pack_tran_dmat(nx, nu, B+i*nx*nu, nx, sBtAt, 0, i*nx);
        blasfeo_pack_tran_dmat(nx, nx, A+i*nx*nx, nx, sBtAt, nu, i*nx);
        blasfeo_pack_tran_dmat(nbg, nx, C+i*nbg*nx, nbg, sCt, 0, i*nbg);

        // blasfeo_pack_dvec(nx, gx, srq, i*(nx+nu)+nu);
        // blasfeo_pack_dvec(nu, gu, srq, i*(nx+nu));
    }
    // blasfeo_pack_dmat(nu, nx, S+N*nx*nu, nu, sStQ, 0, N*nx);
    blasfeo_pack_dmat(nx, nx, Q+N*nx*nx, nx, sStQ, nu, N*nx);   

    /* compute G */
    for(i=0;i<N;i++){
        memcpy(G+(i*N+i)*nx*nu, B+i*nx*nu, nx*nu*sizeof(double));
        blasfeo_dgecp(nu, nx, sBtAt, 0, i*nx, sGt, i*nu, i*nx);
        for (j=i+1;j<N;j++){
            dgemm_(nTrans, nTrans, &nx, &nu, &nx, &one_d, A+j*nx*nx, &nx, G+(i*N+j-1)*nx*nu, &nx, &zero, G+(i*N+j)*nx*nu, &nx);
            blasfeo_dgemm_nn(nu, nx, nx, 1.0, sGt, (j-1)*nu, i*nx, sBtAt, nu, j*nx, 0.0, sGt, j*nu, i*nx, sGt, j*nu, i*nx);
        }
    }

    /* Compute Hc */
    for(i=0;i<N;i++){
        dsymm_(SIDE, UPLO, &nx, &nu, &one_d, Q+N*nx*nx, &nx, G+(i*N+N-1)*nx*nu, &nx, &zero, W+(i*N+N-1)*nx*nu, &nx);
        blasfeo_dgemm_nn(nu, nx, nx, 1.0, sGt, (N-1)*nu, i*nx, sStQ, nu, N*nx, 0.0, sWttmp, 0, 0, sWttmp, 0, 0);
        for(j=N-1;j>i;j--){        
            dgemm_(nTrans, nTrans, &nu, &nu, &nx, &one_d, S+j*nx*nu, &nu, G+(i*N+j-1)*nx*nu, &nx, &zero, Hi, &nu);                     
            dgemm_(Trans, nTrans, &nu, &nu, &nx, &one_d, B+j*nx*nu, &nx, W+(i*N+j)*nx*nu, &nx, &one_d, Hi, &nu);
            Block_Fill(nu, nu, Hi, Hc, j*nu, i*nu, N*nu);
            Block_Fill_Trans(nu, nu, Hi, Hc, i*nu, j*nu, N*nu);
            dgemm_(Trans, nTrans, &nx, &nu, &nx, &one_d, A+j*nx*nx, &nx, W+(i*N+j)*nx*nu, &nx, &zero, W+(i*N+j-1)*nx*nu, &nx); 
            dsymm_(SIDE, UPLO, &nx, &nu, &one_d, Q+j*nx*nx, &nx, G+(i*N+j-1)*nx*nu, &nx, &one_d, W+(i*N+j-1)*nx*nu, &nx);

            blasfeo_dgemm_nt(nu,nu+nx,nx,1.0,sGt,(j-1)*nu,i*nx,sStQ,0,i*nx,0.0,sHtWt,0,0,sHtWt,0,0);
            blasfeo_dgemm_nt(nu,nu+nx,nx,1.0,sWttmp,0,0,sBtAt,0,i*nx,1.0,sHtWt,0,0,sHtWt,0,0);
            blasfeo_dgecp(nu, nu, sHtWt, 0, 0, sHc, i*nu, j*nu);
            blasfeo_dgecp(nu, nx, sHtWt, 0, nu, sWttmp, 0, 0);
            blasfeo_dgetr(nu, nu, sHtWt, 0, 0, sHc, j*nu, i*nu);          
        }
        memcpy(Hi,R+i*nu*nu,nu*nu*sizeof(double));
        dgemm_(Trans, nTrans, &nu, &nu, &nx, &one_d, B+i*nx*nu, &nx, W+(i*N+i)*nx*nu, &nx, &one_d, Hi, &nu);
        Block_Fill(nu, nu, Hi, Hc, i*nu, i*nu, N*nu); 

        blasfeo_dgemm_nt(nu,nu,nx,1.0,sWttmp,0,0,sBtAt,0,i*nx,1.0,sRt,0,i*nu,sHc,i*nu,i*nu);     
    }

    FILE *fHc = fopen("Hc.dat","w");
    d_print_to_file_mat(fHc,N*nu,N*nu,Hc,N*nu);
    fclose(fHc);
    FILE *fsHc = fopen("sHc.dat","w");
    blasfeo_print_to_file_dmat(fsHc, N*nu, N*nu, sHc, 0, 0);
    fclose(fsHc);
    // blasfeo_unpack_dmat(N*nu,N*nu,sHc,0,0,Hc,N*nu);

    /* Compute Cc */
    if (nbg>0){         
        for(i=0;i<N;i++){             
            Block_Fill_Trans(nbg,nu, D+i*nbg*nu, Cc, i*nu, i*nbg, N*nu);
            // blasfeo_pack_tran_dmat(nbg,nu,D+i*nbg*nu,nbg,sCct,i*nu,i*nbg);
            for(j=i+1;j<N;j++){   
                dgemm_(nTrans, nTrans, &nbg, &nu, &nx, &one_d, C+j*nbg*nx, &nbg, G+(i*N+j-1)*nx*nu, &nx, &zero, Cci, &nbg);                   
                Block_Fill_Trans(nbg, nu, Cci, Cc, i*nu, j*nbg, N*nu);

                // blasfeo_dgemm_nt(nu,nu,nx,1.0,sWttmp,0,0,sBtAt,0,i*nx,1.0,sRt,0,i*nu,sHc,i*nu,i*nu); 
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
        dgemv_(nTrans,&nu,&nx,&one_d,S+i*nx*nu,&nu,L+i*nx,&one_i,&one_d,gc+i*nu,&one_i);
        dgemv_(Trans,&nx,&nu,&one_d,B+i*nx*nu,&nx,w+i*nx,&one_i,&one_d,gc+i*nu,&one_i);
         
        memcpy(w+(i-1)*nx, gx+i*nx, nx*sizeof(double));
        dsymv_(UPLO, &nx, &one_d, Q+i*nx*nx, &nx, L+i*nx, &one_i, &one_d, w+(i-1)*nx, &one_i);
        dgemv_(Trans,&nx,&nx,&one_d,A+i*nx*nx,&nx,w+i*nx,&one_i,&one_d,w+(i-1)*nx,&one_i);
    }   
    memcpy(gc,gu,nu*sizeof(double));
    dgemv_(Trans,&nu,&nx,&one_d,S,&nu,L,&one_i,&one_d,gc,&one_i);
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