
#include <stdio.h>
#include <stdlib.h>
#include "string.h"

#include "qp_generation.h"
#include "casadi_wrapper.h"
#include "casadi_src.h"
#include "common.h"

extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern void dsyrk_(char*, char*, int*, int*, double*, double*, int*, double*, double*, int*);
extern void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

qp_in* qp_in_create(model_size *size)
{
    int nx=size->nx;
    int nu=size->nu;
    int ny=size->ny;
    int nyN=size->nyN;
    int np=size->np;
    int nbx = size->nbx;
    int nbu = size->nbu;
    int nbg = size->nbg;
    int nbgN = size->nbgN;
    int N = size->N;

    qp_in *in = (qp_in*)malloc(sizeof(qp_in));

    in->x = (double*)calloc((N+1)*nx,sizeof(double));
    in->u = (double*)calloc(N*nu,sizeof(double));
    in->y = (double*)calloc(N*ny,sizeof(double));
    in->yN = (double*)calloc(nyN,sizeof(double));
    in->W = (double*)calloc(ny*ny,sizeof(double));
    in->WN = (double*)calloc(nyN*nyN,sizeof(double));
    in->p = (double*)calloc((N+1)*np,sizeof(double));
    in->lb = (double*)calloc(nbu+nbx,sizeof(double));
    in->ub = (double*)calloc(nbu+nbx,sizeof(double));
    in->lbg = (double*)calloc(nbg,sizeof(double));
    in->ubg = (double*)calloc(nbg,sizeof(double));
    in->lbgN = (double*)calloc(nbgN,sizeof(double));
    in->ubgN = (double*)calloc(nbgN,sizeof(double));

    return in;
}

qp_problem* qp_problem_create(model_size *size)
{
    int nx=size->nx;
    int nu=size->nu;
    int nbx = size->nbx;
    int nbu = size->nbu;
    int nbg = size->nbg;
    int nbgN = size->nbgN;
    int N = size->N;

    qp_problem *qp = (qp_problem*)malloc(sizeof(qp_problem));

    qp->Q = (double*)calloc((N+1)*nx*nx,sizeof(double));
    qp->S = (double*)calloc((N+1)*nx*nu,sizeof(double));
    qp->R = (double*)calloc((N+1)*nu*nu,sizeof(double));
    qp->A = (double*)calloc(N*nx*nx,sizeof(double));
    qp->B = (double*)calloc(N*nx*nu,sizeof(double));
    qp->b = (double*)calloc(N*nx,sizeof(double));
    qp->C = (double*)calloc(N*nbg*nx,sizeof(double));
    qp->CN = (double*)calloc(nbgN*nx,sizeof(double));
    qp->D = (double*)calloc(N*nbg*nu,sizeof(double));
    qp->DN = (double*)calloc(nbgN*nu,sizeof(double));
    qp->gx = (double*)calloc(nx*(N+1),sizeof(double));
    qp->gu = (double*)calloc(nu*(N),sizeof(double));
    qp->lb = (double*)calloc((nbu+nbx)*N+nbx,sizeof(double));
    qp->ub = (double*)calloc((nbu+nbx)*N+nbx,sizeof(double));
    qp->lb_g = (double*)calloc(nbg*N+nbgN,sizeof(double));
    qp->ub_g = (double*)calloc(nbg*N+nbgN,sizeof(double));

    return qp;
}

qp_generation_workspace* qp_generation_workspace_create(model_size *size)
{
    int nx=size->nx;
    int nu=size->nu;
    int ny=size->ny;
    int nyN=size->nyN;

    qp_generation_workspace *work = (qp_generation_workspace*)malloc(sizeof(qp_generation_workspace));

    work->Jac = (double**)malloc(2*sizeof(double*));
    work->Jac[0] = (double*)calloc(ny*nx,sizeof(double));
    work->Jac[1] = (double*)calloc(ny*nu,sizeof(double));
    work->Jac_N = (double*)calloc(nyN*nx,sizeof(double));

    return work;
}

qp_out* qp_out_create(model_size *size)
{
    int nx=size->nx;
    int nu=size->nu; 
    int nbg=size->nbg;
    int nbgN=size->nbgN;  
    int N = size->N;

    qp_out *out = (qp_out*)malloc(sizeof(qp_out));

    out->dx = (double*)calloc((N+1)*nx,sizeof(double));
    out->du = (double*)calloc(N*nu,sizeof(double));
    out->lam = (double*)calloc(nx*(N+1),sizeof(double));
    out->mu = (double*)calloc((N*nbg+nbgN),sizeof(double));

    return out;
}

int qp_generation(qp_in *in, model_size *size, 
    qp_problem *qp, qp_out *out, qp_generation_workspace *qp_work)
{
    int nx = size->nx;
    int nu = size->nu;
    int ny = size->ny;
    int nyN = size->nyN;
    int np = size->np;
    int nbx = size->nbx;
    int nbu = size->nbu;
    int nbg = size->nbg;
    int nbgN = size->nbgN;
    int N = size->N;
    int *nbx_idx = size->nbx_idx;
    int *nbu_idx = size->nbu_idx;

    double *x = in->x;
    double *u = in->u;
    double *y = in->y;
    double *yN = in->yN; 
    double *W = in->W;
    double *WN = in->WN;
    double *p = in->p;
    double *lb = in->lb;
    double *ub = in->ub;
    double *lbg = in->lbg;
    double *ubg = in->ubg;
    double *lbgN = in->lbgN;
    double *ubgN = in->ubgN;
       
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
    double *lb_dux = qp->lb;
    double *ub_dux = qp->ub;
    double *lb_g = qp->lb_g;
    double *ub_g = qp->ub_g;
   
    double **Jac = qp_work->Jac;
    double *Jac_N = qp_work->Jac_N; 
        
    int i=0,j=0;
    char *nTrans = "N", *Trans="T", *UPLO="L";
    double one_d = 1.0, zero = 0.0;
        
    // allocate array of pointers
    double *Sens[2];    
    double *Cons[2];
      
    double *casadi_in[6];
    double *casadi_out[2];    
    casadi_in[4] = W;
                      
    // start loop
    for(i=0;i<N;i++){
        casadi_in[0] = x+i*nx;
        casadi_in[1] = u+i*nu;
        casadi_in[2] = p+i*np;
        casadi_in[3] = y+i*ny;
        
        // control bounds
        for (j=0;j<nbu;j++){
            lb_dux[i*(nbu+nbx)+j] = lb[j]-u[i*nu+nbu_idx[j]];
            ub_dux[i*(nbu+nbx)+j] = ub[j]-u[i*nu+nbu_idx[j]];
        }

        // state bounds
        for (j=0;j<nbx;j++){
            lb_dux[i*(nbu+nbx)+nbu+j] = lb[nbu+j]-x[i*nx+nbx_idx[j]];
            ub_dux[i*(nbu+nbx)+nbu+j] = ub[nbu+j]-x[i*nx+nbx_idx[j]];
        }
        
        // integration                             
        casadi_out[0] = b+i*nx;
        F_Fun(casadi_in, casadi_out);
      
        // sensitivity computation
        Sens[0] = A + i*nx*nx;
        Sens[1] = B + i*nx*nu;
        D_Fun(casadi_in, Sens);	
                           
        // equality residual        
        for (j=0;j<nx;j++)
            b[i*nx+j] -= x[(i+1)*nx+j];
       
        // Hessian
        Ji_Fun(casadi_in, Jac);
        dsyrk_(UPLO, Trans, &nx, &ny, &one_d, Jac[0], &ny, &zero, Q+i*nx*nx, &nx);
        dgemm_(Trans, nTrans, &nu, &nx, &ny, &one_d, Jac[1], &ny, Jac[0], &ny, &zero, S+i*nx*nu, &nu); // S^T
        dsyrk_(UPLO, Trans, &nu, &ny, &one_d, Jac[1], &ny, &zero, R+i*nu*nu, &nu);
                
        // gradient
        casadi_out[0] = gx+i*nx;
        casadi_out[1] = gu+i*nu;
        gi_Fun(casadi_in, casadi_out);
                
        // constraint residual
        if (nbg>0){
            casadi_out[0] = lb_g + i*nbg;
            path_con_Fun(casadi_in, casadi_out);
            for (j=0;j<nbg;j++){
                ub_g[i*nbg+j] = ubg[j] - casadi_out[0][j];
                casadi_out[0][j] = lbg[j] - casadi_out[0][j];            
            }
        
            // constraint Jacobian
            Cons[0] = C+i*nbg*nx;
            Cons[1] = D+i*nbg*nu;
            Ci_Fun(casadi_in, Cons);
        }
    }
    
    // terminal data
    casadi_in[0] = x+N*nx;
    casadi_in[1] = p+N*np;
    casadi_in[2] = yN;
    casadi_in[3] = WN;

    for (j=0;j<nbx;j++){
        lb_dux[N*(nbu+nbx)+j] = lb[nbu+j]-x[N*nx+nbx_idx[j]];
        ub_dux[N*(nbu+nbx)+j] = ub[nbu+j]-x[N*nx+nbx_idx[j]];
    }
    
    JN_Fun(casadi_in, Jac_N);
    dsyrk_(UPLO, Trans, &nx, &nyN, &one_d, Jac_N, &nyN, &zero, Q+N*nx*nx, &nx);
    
    casadi_out[0] = gx+N*nx;
    gN_Fun(casadi_in, casadi_out);

    if (nbgN>0){
        casadi_out[0] = lb_g + N*nbg;
        path_con_N_Fun(casadi_in, casadi_out);
        for (j=0;j<nbgN;j++){
            ub_g[i*nbg+j] = ubgN[j] - casadi_out[0][j];
            casadi_out[0][j] = lbgN[j] - casadi_out[0][j];            
        }

        CN_Fun(casadi_in, &CN);
    }
    
    
    // printf("QN=\n");
    // print_matrix(Q+N*nx*nx, nx, nx);

    // printf("S1=\n");
    // print_matrix(S+nx*nu, nu, nx);

    // printf("R1=\n");
    // print_matrix(R+nu*nu, nu, nu);

    // printf("A1=\n");
    // print_matrix(A+nx*nx, nx, nx);

    // printf("B1=\n");
    // print_matrix(B+nx*nu, nx, nu);

    // printf("b1=\n");
    // print_vector(b+nx, nx);

    // printf("gxN=\n");
    // print_vector_trans(gx+N*nx, nx);

    // printf("C1=\n");
    // print_matrix(C+nbg*nx, nbg, nx);

    // printf("D1=\n");
    // print_matrix(D+nbg*nu, nbg, nu);

    // printf("CN=\n");
    // print_matrix(CN, nbgN, nx);

    // printf("lb=\n");
    // print_vector_trans(lb_dux, N*(nbu+nbx)+nbx);

    // printf("ub=\n");
    // print_vector_trans(ub_dux, N*(nbu+nbx)+nbx);

    // printf("lb_g=\n");
    // print_vector_trans(lb_g, N*nbg+nbgN);

    // printf("ub_g=\n");
    // print_vector_trans(ub_g, N*nbg+nbgN);
	
    return 0;
    
}

int expand(model_size *size, qp_problem *qp, qp_out *out)
 {
    int nx = size->nx;
    int nu = size->nu;
    int N = size->N;

    double *dx = out->dx;
    double *du = out->du;

    double *A = qp->A;
    double *B = qp->B;
    double *b = qp->b;

    int i;
    char *nTrans = "N";
    double one_d = 1.0;
    int one_i = 1;

    for(i=0;i<N;i++){
        memcpy(dx+(i+1)*nx, b+i*nx, nx*sizeof(double));          
        dgemv_(nTrans,&nx,&nx,&one_d,A+i*nx*nx,&nx,dx+i*nx,&one_i,&one_d,dx+(i+1)*nx,&one_i);
        dgemv_(nTrans,&nx,&nu,&one_d,B+i*nx*nu,&nx,du+i*nu,&one_i,&one_d,dx+(i+1)*nx,&one_i);
    }

    return 0;
 }

void qp_in_free(qp_in *in){
    free(in->x);
    free(in->u);
    free(in->y);
    free(in->yN);
    free(in->W);
    free(in->WN);
    free(in->p);
    free(in->lb);
    free(in->ub);
    free(in->lbg);
    free(in->ubg);
    free(in->lbgN);
    free(in->ubgN);

    free(in);
}

void qp_problem_free(qp_problem *qp){
    free(qp->Q);
    free(qp->S);
    free(qp->R);
    free(qp->A);
    free(qp->B);
    free(qp->b);
    free(qp->gx);
    free(qp->gu);
    free(qp->C);
    free(qp->CN);
    free(qp->D);
    free(qp->lb);
    free(qp->ub);
    free(qp->lb_g);
    free(qp->ub_g);

    free(qp);
}

void qp_generation_workspace_free(qp_generation_workspace *work){
    free(work->Jac[0]);
    free(work->Jac[1]);
    free(work->Jac);
    free(work->Jac_N);

    free(work);
}

void qp_out_free(qp_out *out){
    free(out->dx);
    free(out->du);
    free(out->lam);
    free(out->mu);

    free(out);
}