#include <stdlib.h>
#include <stdio.h>

#include "rti_step.h"
#include "common.h"

#define PI 3.141592653589793

int main()
{
    int i,j;

    // Chain of Masses
    // int nx = 27;
    // int nu = 3;
    // int ny = 18;
    // int nyN = 15;
    // int np = 1;
    // int nbu = 3;
    // int nbg = 0;
    // int nbgN = 0;
    // int N = 40;
    // int Ns = 20;
    
    // InvertedPendulum
    int nx = 4;
    int nu = 1;
    int ny = 5;
    int nyN = 4;
    int np = 1;
    int nbg = 1;
    int nbgN = 1;
    int N = 40;

    /* Initialize initial condition */

    model_size size;
    size.nx = nx;
    size.nu = nu;
    size.ny = ny;
    size.nyN = nyN;
    size.np = np;
    size.nbg = nbg;
    size.nbgN = nbgN;
    size.N = N;
    
    double *x0 = calloc(nx, sizeof(double));

    // int n=5; // fixed
    // for(j=0;j<n;j++)
    //     x0[j] = 7.5*(j+1)/n;

    // for(i=0;i<N;i++){
    //     for(j=0;j<n;j++)
    //         in->x[nx*i+j] = x0[j];

    //     for (j=0;j<nu;j++)
    //         in->u[nu*i+j] = 0;

    //     in->y[ny*i+0] = 7.5;
    // }

    // for (i=0;i<3;i++){
    //     in->W[i*ny+i] = 25;
    //     in->WN[i*nyN+i] = 25;
    // }
    // for (i=3;i<3+3*(n-1);i++){
    //     in->W[i*ny+i] = 0.25;
    //     in->WN[i*nyN+i] = 0.25;
    // }
    // for (i=3+3*(n-1);i<ny;i++)
    //     in->W[i*ny+i] = 0.1;

    // in->yN[0]=7.5;

    // for (i=0;i<nbu;i++){
    //     in->lbu[i] = -1;
    //     in->ubu[i] = 1;
    // }

    // for (i=0;i<nbg;i++){
    //     in->lbg[i] = 0;
    //     in->ubg[i] = 0;
    // }

    // for (i=0;i<nbgN;i++){
    //     in->lbgN[i] = 0;
    //     in->ubgN[i] = 0;
    // }

    x0[1] = PI;

    rti_step_workspace *rti_work= rti_step_workspace_create(&size);

    qp_in *in = rti_work->in;

    for(i=0;i<N;i++){
        for(j=0;j<nx;j++)
            in->x[nx*i+j] = x0[j];

        for (j=0;j<nu;j++)
            in->u[nu*i+j] = 0;
    }
    for(j=0;j<nx;j++)
        in->x[nx*N+j] = x0[j];

    in->W[0] = 10; in->W[1*ny+1]=10; in->W[2*ny+2]=0.1; in->W[3*ny+3]=0.1; in->W[4*ny+4]=0.01; 
    in->WN[0] = 10; in->WN[1*nyN+1]=10; in->WN[2*nyN+2]=0.1; in->WN[3*nyN+3]=0.1; 

    for (i=0;i<nu;i++){
        in->lbu[i] = -20;
        in->ubu[i] = 20;
    }

    for (i=0;i<nbg;i++){
        in->lbg[i] = -2;
        in->ubg[i] = 2;
    }

    for (i=0;i<nbgN;i++){
        in->lbgN[i] = -2;
        in->ubgN[i] = 2;
    }

    /* call rti step */

    
    rti_work->sample = 0;

    rti_step_init(&size, rti_work);
    
    int Nrep = 20;
    double cpt_rti=0;
    double cpt_prep=0, cpt_cond=0, cpt_qp_dense=0, cpt_qp_sparse=0;
    CMPC_timer t;
    
    t=CMPC_tic(t);
    for(i=0;i<Nrep;i++){
        rti_step(x0, &size, rti_work);
        cpt_prep += rti_work->cpt_prep;
        cpt_cond += rti_work->cpt_cond;
        cpt_qp_dense += rti_work->cpt_qp_dense;
        cpt_qp_sparse += rti_work->cpt_qp_sparse;

        rti_work->sample ++;
    }
    t=CMPC_toc(t);

    cpt_rti = t.t;
    
    printf("CPU Time:");
    printf("%8.6f ms\n", cpt_rti*1E3/Nrep);

    printf("QP preparation time:");
    printf("%8.6f ms\n", cpt_prep*1E3/Nrep);

    printf("Condensing time:");
    printf("%8.6f ms\n", cpt_cond*1E3/Nrep);

    printf("QORE time:");
    printf("%8.6f ms\n", cpt_qp_dense*1E3/Nrep);

    printf("HPIPM time:");
    printf("%8.6f ms\n", cpt_qp_sparse*1E3/Nrep);

    // print_vector(rti_work->qore_work->sol_qore, size.N*size.nu);

    /* free */
    free(x0);    
    rti_step_workspace_free(rti_work);
  
    return 0;
    
}