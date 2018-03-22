#include <stdlib.h>
#include <stdio.h>
#include "string.h"

#include "rti_step.h"
#include "common.h"
#include "qpsolver_hpipm_ocp.h"
#include "qp_generation.h"
#include "full_condensing.h"

extern void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

rti_step_workspace* rti_step_workspace_create(model_size *size)
{

    rti_step_workspace *work = (rti_step_workspace*)malloc(sizeof(rti_step_workspace));

    work->in = qp_in_create(size);
    work->qp = qp_problem_create(size);
    work->qp_work = qp_generation_workspace_create(size);
    work->out = qp_out_create(size);
    work->hpipm_ocp_work = qpsolver_hpipm_ocp_workspace_create(size);

    work->full_condensing_work = full_condensing_workspace_create(size);
    work->qore_work = qpsolver_qore_workspace_create(size);


    work->dx0 = (double*)malloc(size->nx * sizeof(double));
    work->u_opt = (double*)calloc(size->nu , sizeof(double));

    return work;
}

void rti_step_init(model_size *size, rti_step_workspace *rti_work)
{
    qpsolver_hpipm_ocp_workspace_init(size, rti_work->qp,
        rti_work->out, rti_work->hpipm_ocp_work);
    
}

int rti_step(double *x0, model_size *size, rti_opt *opt, rti_step_workspace *rti_work)
{
    CMPC_timer tprep;
    CMPC_timer tcond;
    CMPC_timer tqore;
    CMPC_timer thpipm;

    int jj;

    int nx = size->nx;
    int nu = size->nu;
    int N = size->N;

    qp_in *in=rti_work->in;
    qp_problem *qp=rti_work->qp;
    qp_generation_workspace *qp_work=rti_work->qp_work;
    qp_out *out=rti_work->out;
    qpsolver_hpipm_ocp_workspace *hpipm_ocp_work = rti_work->hpipm_ocp_work;
    full_condensing_workspace *full_condensing_work = rti_work->full_condensing_work;
    qpsolver_qore_workspace *qore_work = rti_work->qore_work;

    qp_work->sample = rti_work->sample;

    tprep=CMPC_tic(tprep);
    qp_generation(in, size, qp, out, qp_work);
    tprep=CMPC_toc(tprep);
    rti_work->cpt_prep = tprep.t;

    set_zeros((N+1)*nx, out->dx);
    for(jj=0;jj<nx;jj++)
        out->dx[jj] = x0[jj] - in->x[jj];
    set_zeros(N*nu, out->du);

    double *dx0 = rti_work->dx0;
    for(jj=0;jj<nx;jj++)
        dx0[jj] = x0[jj] - in->x[jj];

    /* call qp solver */
    switch (opt->qpsolver){
        case 0:
            rti_work->cpt_qp_sparse =0;
            tcond=CMPC_tic(tcond);
            full_condensing(dx0, size, qp, full_condensing_work);
            tcond=CMPC_toc(tcond);
            rti_work->cpt_cond = tcond.t;
        
            tqore=CMPC_tic(tqore);
            qpsolver_qore(size, qp, full_condensing_work,
                qore_work, out, rti_work->sample);
            tqore=CMPC_toc(tqore);
            rti_work->cpt_qp_dense = tqore.t;

            memcpy(out->dx, dx0, nx*sizeof(double));
            expand(size, qp, out);

            break;
        case 1:
            rti_work->cpt_cond = 0;
            rti_work->cpt_qp_dense =0;
            thpipm=CMPC_tic(thpipm);
            qpsolver_hpipm_ocp(size, qp, out, hpipm_ocp_work);
            thpipm=CMPC_toc(thpipm);
            rti_work->cpt_qp_sparse = thpipm.t;
            break;
        default:
            printf("Please choose a supported QP solver!");
    }

    /* update solution */
    for(jj=0;jj<nx*(N+1);jj++)
        in->x[jj] += out->dx[jj];

    for(jj=0;jj<nu*N;jj++)
        in->u[jj] += out->du[jj]; 

    /* store the first optimal input */
    for(jj=0;jj<nu;jj++)
        rti_work->u_opt[jj] = in->u[jj];

    /* apply shifting */
    if(opt->shifting){
        for(jj=0;jj<nx*N;jj++)
            in->x[jj] = in->x[jj+nx];
        for(jj=0;jj<nu*(N-1);jj++)
            in->u[jj] = in->u[jj+nu];
    }


    return 0;
}


void rti_step_workspace_free(rti_step_workspace *work)
{
    qp_in_free(work->in);
    qp_problem_free(work->qp);
    qp_generation_workspace_free(work->qp_work);
    qp_out_free(work->out);
    qpsolver_hpipm_ocp_workspace_free(work->hpipm_ocp_work);
    full_condensing_workspace_free(work->full_condensing_work);
    qpsolver_qore_workspace_free(work->qore_work);

    free(work->dx0);
    free(work->u_opt);
  
    free(work);
}