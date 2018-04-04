#include <stdlib.h>
#include <stdio.h>

#include "qpsolver_hpipm_pcond.h"
#include "qp_generation.h"
#include "common.h"

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_v_aux_ext_dep.h>
#include <blasfeo_d_aux_ext_dep.h>
#include <blasfeo_i_aux_ext_dep.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include "hpipm_d_ocp_qp_dim.h"
#include "hpipm_d_ocp_qp.h"
#include "hpipm_d_ocp_qp_sol.h"
#include "hpipm_d_ocp_qp_ipm.h"
#include "hpipm_d_part_cond.h"

qpsolver_hpipm_pcond_workspace* qpsolver_hpipm_pcond_workspace_create(model_size *size)
{
    int nx=size->nx;
    int nu=size->nu;
    int nbg = size->nbg;
    int N = size->N;
    int N2 = size->N2;

    qpsolver_hpipm_pcond_workspace *work = (qpsolver_hpipm_pcond_workspace*)malloc(sizeof(qpsolver_hpipm_pcond_workspace));

    int dim_size = d_memsize_ocp_qp_dim(N);
	work->dim_mem = calloc(dim_size,1);
    work->ocp_dim = malloc(sizeof(struct d_ocp_qp_dim));
    d_create_ocp_qp_dim(N, work->ocp_dim, work->dim_mem);

    int dim_size2 = d_memsize_ocp_qp_dim(N2);
	work->dim_mem2 = calloc(dim_size2,1);
    work->ocp_dim2 = malloc(sizeof(struct d_ocp_qp_dim));
    d_create_ocp_qp_dim(N2, work->ocp_dim2, work->dim_mem2);

    work->nx_v = (int*)malloc((N+1)*sizeof(int));
    work->nu_v = (int*)malloc((N+1)*sizeof(int));
    work->nbu_v = (int*)malloc((N+1)*sizeof(int));
    work->nbx_v = (int*)malloc((N+1)*sizeof(int));
    work->ng_v = (int*)malloc((N+1)*sizeof(int));
    work->ns_v = (int*)malloc((N+1)*sizeof(int));

    work->nx_v2 = (int*)malloc((N2+1)*sizeof(int));
    work->nu_v2 = (int*)malloc((N2+1)*sizeof(int));
    work->nbu_v2 = (int*)malloc((N2+1)*sizeof(int));
    work->nbx_v2 = (int*)malloc((N2+1)*sizeof(int));
    work->ng_v2 = (int*)malloc((N2+1)*sizeof(int));
    work->ns_v2 = (int*)malloc((N2+1)*sizeof(int));
    work->block_size = (int*)malloc((N2+1)*sizeof(int));

    work->b0 = (double*)malloc(nx*sizeof(double));
    work->r0 = (double*)malloc(nu*sizeof(double));
    work->lg0 = (double*)malloc(nbg*sizeof(double));
    work->ug0 = (double*)malloc(nbg*sizeof(double));

    work->hA = (double**)malloc(N*sizeof(double*));
    work->hB = (double**)malloc(N*sizeof(double*));
    work->hb = (double**)malloc(N*sizeof(double*));
    work->hQ = (double**)malloc((N+1)*sizeof(double*));
    work->hS = (double**)malloc((N+1)*sizeof(double*));
    work->hR = (double**)malloc((N+1)*sizeof(double*));
    work->hq = (double**)malloc((N+1)*sizeof(double*));
    work->hr = (double**)malloc((N+1)*sizeof(double*));
    work->hlb = (double**)malloc((N+1)*sizeof(double*));
    work->hub = (double**)malloc((N+1)*sizeof(double*));
    work->hC = (double**)malloc((N+1)*sizeof(double*));
    work->hD = (double**)malloc((N+1)*sizeof(double*));
    work->hlg = (double**)malloc((N+1)*sizeof(double*));
    work->hug = (double**)malloc((N+1)*sizeof(double*));
    work->hidxb = (int**)malloc((N+1)*sizeof(int*));

    work->hx = (double**)malloc((N+1)*sizeof(double*));
    work->hu = (double**)malloc((N+1)*sizeof(double*));
    work->hpi = (double**)malloc(N*sizeof(double*));
    work->hlam_lb = (double**)malloc((N+1)*sizeof(double*));
    work->hlam_ub = (double**)malloc((N+1)*sizeof(double*));
    work->hlam_lg = (double**)malloc((N+1)*sizeof(double*));
    work->hlam_ug = (double**)malloc((N+1)*sizeof(double*));

    work->idxc = (int*)malloc(N*sizeof(int));

    return work;
}

void qpsolver_hpipm_pcond_workspace_init(model_size *size,
     qp_problem *qp, qp_out *out, qpsolver_hpipm_pcond_workspace *work)
{
    int ii,jj;
    
    int nx = size->nx;
    int nu = size->nu;
    int nbx = size->nbx;
    int nbu = size->nbu;
    int nbg = size->nbg;
    int nbgN = size->nbgN;
    int N = size->N;
    int N2 = size->N2;
    int *nbx_idx = size->nbx_idx;
    int *nbu_idx = size->nbu_idx;

    int *nx_v = work->nx_v;
    int *nu_v = work->nu_v;
    int *nbx_v =work->nbx_v;
    int *nbu_v =work->nbu_v;
    int *ng_v = work->ng_v;
    int *ns_v = work->ns_v;
    int **hidxb = work->hidxb;
    double **hlam_lb = work->hlam_lb;
    double **hlam_ub = work->hlam_ub;
    double **hlam_lg = work->hlam_lg;
    double **hlam_ug = work->hlam_ug;

    int *nx_v2 = work->nx_v2;
    int *nu_v2 = work->nu_v2;
    int *nbx_v2 =work->nbx_v2;
    int *nbu_v2 =work->nbu_v2;
    int *ng_v2 = work->ng_v2;
    int *ns_v2 = work->ns_v2;
    int *block_size = work->block_size;
    
    nx_v[0] = 0;
	for(ii=1; ii<=N; ii++)
		nx_v[ii] = nx;

    for(ii=0; ii<N; ii++)
		nu_v[ii] = nu;
	nu_v[N] = 0;

    for(ii=0; ii<N; ii++)
        nbu_v[ii] = nbu;
	nbu_v[N] = 0;

    nbx_v[0] = 0;
	for(ii=1; ii<=N; ii++)
        nbx_v[ii] = nbx;

    for(ii=0; ii<N; ii++)
		ng_v[ii] = nbg;
	ng_v[N] = nbgN;

    for(ii=0; ii<=N; ii++)
		ns_v[ii] = 0;

    for(ii=0;ii<=N;ii++){
        int_zeros(hidxb+ii, nbu_v[ii]+nbx_v[ii], 1);
        d_zeros(hlam_lb+ii, nbu_v[ii]+nbx_v[ii], 1);
        d_zeros(hlam_ub+ii, nbu_v[ii]+nbx_v[ii], 1);
        d_zeros(hlam_lg+ii, ng_v[ii], 1);
        d_zeros(hlam_ug+ii, ng_v[ii], 1);
    }


    d_cvt_int_to_ocp_qp_dim(N, nx_v, nu_v, nbx_v, nbu_v, ng_v, ns_v, work->ocp_dim);
    d_cvt_int_to_ocp_qp_dim(N2, nx_v2, nu_v2, nbx_v2, nbu_v2, ng_v2, ns_v2, work->ocp_dim2);

    d_compute_block_size_cond_qp_ocp2ocp(N, N2, block_size);

    d_compute_qp_dim_ocp2ocp(work->ocp_dim, block_size, work->ocp_dim2);

    // for(ii=0; ii<=N2; ii++)
		// printf("\n%d %d %d %d %d\n", work->ocp_dim2->nx[ii], work->ocp_dim2->nu[ii], work->ocp_dim2->nbx[ii], work->ocp_dim2->nbu[ii], work->ocp_dim2->ng[ii]);

    int qp_size = d_memsize_ocp_qp(work->ocp_dim);
	work->qp_mem = calloc(qp_size,1);
    work->ocp_qp = malloc(sizeof(struct d_ocp_qp));
    d_create_ocp_qp(work->ocp_dim, work->ocp_qp, work->qp_mem);

    int qp_size2 = d_memsize_ocp_qp(work->ocp_dim2);
	work->qp_mem2 = calloc(qp_size2,1);
    work->ocp_qp2 = malloc(sizeof(struct d_ocp_qp));
    d_create_ocp_qp(work->ocp_dim2, work->ocp_qp2, work->qp_mem2);

    int qp_sol_size = d_memsize_ocp_qp_sol(work->ocp_dim);
	work->sol_mem = calloc(qp_sol_size,1);
    work->ocp_sol = malloc(sizeof(struct d_ocp_qp_sol));
	d_create_ocp_qp_sol(work->ocp_dim, work->ocp_sol, work->sol_mem);

    int qp_sol_size2 = d_memsize_ocp_qp_sol(work->ocp_dim2);
	work->sol_mem2 = calloc(qp_sol_size2,1);
    work->ocp_sol2 = malloc(sizeof(struct d_ocp_qp_sol));
	d_create_ocp_qp_sol(work->ocp_dim2, work->ocp_sol2, work->sol_mem2);

    int ipm_arg_size = d_memsize_ocp_qp_ipm_arg(work->ocp_dim2);
	work->arg_mem = calloc(ipm_arg_size,1);
    work->ocp_arg = malloc(sizeof(struct d_ocp_qp_ipm_arg));
	d_create_ocp_qp_ipm_arg(work->ocp_dim2, work->ocp_arg, work->arg_mem);

    d_set_default_ocp_qp_ipm_arg(work->ocp_arg);
    work->ocp_arg->res_g_max = 1e-4;
	work->ocp_arg->res_b_max = 1e-6;
	work->ocp_arg->res_d_max = 1e-6;
	work->ocp_arg->res_m_max = 1e-6;
	work->ocp_arg->mu0 = 10.0;
	work->ocp_arg->iter_max = 20;

    int ipm_size = d_memsize_ocp_qp_ipm(work->ocp_dim2, work->ocp_arg);
	work->ipm_mem = calloc(ipm_size,1);
    work->ocp_workspace = malloc(sizeof(struct d_ocp_qp_ipm_workspace));
	d_create_ocp_qp_ipm(work->ocp_dim2, work->ocp_arg, work->ocp_workspace, work->ipm_mem);	

    int part_cond_arg_size = d_memsize_cond_qp_ocp2ocp_arg(N2);
    work->part_cond_arg_mem = calloc(part_cond_arg_size,1);
    work->part_cond_arg = malloc(sizeof(struct d_cond_qp_ocp2ocp_arg));
    d_create_cond_qp_ocp2ocp_arg(N2, work->part_cond_arg, work->part_cond_arg_mem);
    d_set_default_cond_qp_ocp2ocp_arg(N2, work->part_cond_arg);

    int part_cond_size = d_memsize_cond_qp_ocp2ocp(work->ocp_dim, block_size, work->ocp_dim2, work->part_cond_arg);
	work->part_cond_mem = calloc(part_cond_size,1);
	work->part_cond_ws = malloc(sizeof(struct d_cond_qp_ocp2ocp_workspace));
	d_create_cond_qp_ocp2ocp(work->ocp_dim, block_size, work->ocp_dim2, work->part_cond_arg, work->part_cond_ws, work->part_cond_mem);

    for(ii=0; ii<N; ii++)
		work->hA[ii] = qp->A+ii*nx*nx;

	for(ii=0; ii<N; ii++)
		work->hB[ii] = qp->B+ii*nx*nu;

	work->hb[0] = work->b0;
	for(ii=1; ii<N; ii++)
		work->hb[ii] = qp->b+ii*nx;

	for(ii=0; ii<=N; ii++)
		work->hQ[ii] = qp->Q+ii*nx*nx;
        
	for(ii=0; ii<=N; ii++)
		work->hS[ii] = qp->S+ii*nx*nu;

    for(ii=0; ii<=N; ii++)
		work->hR[ii] = qp->R+ii*nu*nu;

	for(ii=0; ii<=N; ii++)
		work->hq[ii] = qp->gx+ii*nx;
		
	work->hr[0] = work->r0;
	for(ii=1; ii<=N; ii++)
		work->hr[ii] = qp->gu+ii*nu;

	for(ii=0; ii<=N; ii++)
		work->hlb[ii] = qp->lb+ii*(nbu+nbx);

	for(ii=0; ii<=N; ii++)
		work->hub[ii] = qp->ub+ii*(nbu+nbx);

	for(ii=0; ii<N; ii++)
		work->hC[ii] = qp->C+ii*nbg*nx;
	work->hC[N] = qp->CN;

	for(ii=0; ii<N; ii++)
		work->hD[ii] = qp->D+ii*nbg*nu;
    work->hD[N] = qp->DN;

	work->hlg[0] = work->lg0;
	for(ii=1; ii<=N; ii++)
		work->hlg[ii] = qp->lb_g+ii*nbg;

	work->hug[0] = work->ug0;
	for(ii=1; ii<=N; ii++)
		work->hug[ii] = qp->ub_g+ii*nbg;

    for(ii=0; ii<=N; ii++)
	{
		for(jj=0; jj<nbu_v[ii]; jj++)
			work->hidxb[ii][jj] = nbu_idx[jj];
        for(jj=0; jj<nbx_v[ii]; jj++)
			work->hidxb[ii][nbu_v[ii]+jj] = nbu_v[ii]+nbx_idx[jj];
	}
	
	for(ii=0; ii<=N; ii++)
		work->hx[ii] = out->dx+ii*nx;

	for(ii=0; ii<N; ii++)
		work->hu[ii] = out->du+ii*nu;
	
	for(ii=0; ii<N; ii++)
		work->hpi[ii] = out->lam+(ii+1)*nx;
    
}

int qpsolver_hpipm_pcond(model_size *size, qp_problem *qp, qp_out *out, qpsolver_hpipm_pcond_workspace *hpipm_pcond_work)
{
    Timer tcond, tqp;
    
    int nx = size->nx;
    int nu = size->nu;   
    int nbg = size->nbg;
    // int N2 = size->N2;   
    
    double **hA = hpipm_pcond_work->hA;
    double **hB = hpipm_pcond_work->hB;
    double **hb = hpipm_pcond_work->hb;
    double **hQ = hpipm_pcond_work->hQ;
    double **hS = hpipm_pcond_work->hS;
    double **hR = hpipm_pcond_work->hR;
    double **hq = hpipm_pcond_work->hq;
    double **hr = hpipm_pcond_work->hr;
    double **hlb = hpipm_pcond_work->hlb;
    double **hub = hpipm_pcond_work->hub;
    double **hC = hpipm_pcond_work->hC;
    double **hD = hpipm_pcond_work->hD;
    double **hlg = hpipm_pcond_work->hlg;
    double **hug = hpipm_pcond_work->hug;
    int **hidxb = hpipm_pcond_work->hidxb;
 
    double **hx = hpipm_pcond_work->hx;
    double **hu = hpipm_pcond_work->hu;
    double **hpi = hpipm_pcond_work->hpi;
    double **hlam_lb = hpipm_pcond_work->hlam_lb;
    double **hlam_ub = hpipm_pcond_work->hlam_ub;
    double **hlam_lg = hpipm_pcond_work->hlam_lg;
    double **hlam_ug = hpipm_pcond_work->hlam_ug;

    dgemv_n_3l(nx, nx, 1.0, qp->A, nx, out->dx, 1.0, qp->b, hpipm_pcond_work->b0); // b0 = 1.0*b+1.0*A*x	
	dgemv_n_3l(nu, nx, 1.0, qp->S, nu, out->dx, 1.0, qp->gu, hpipm_pcond_work->r0); // r0 = 1.0*r+1.0*S'*x
	dgemv_n_3l(nbg, nx, -1.0, qp->C, nbg, out->dx, 1.0, qp->lb_g, hpipm_pcond_work->lg0); // lg0 = 1.0*lg-1.0*C*x
	dgemv_n_3l(nbg, nx, -1.0, qp->C, nbg, out->dx, 1.0, qp->ub_g, hpipm_pcond_work->ug0); // ug0 = 1.0*ug-1.0*C*x

    d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, NULL, NULL, NULL, NULL, NULL, hpipm_pcond_work->ocp_qp);

    // printf("before condensing\n");

    // int *nx_v2 = hpipm_pcond_work->nx_v2;
    // int *nu_v2 = hpipm_pcond_work->nu_v2;
    // int *nbx_v2 =hpipm_pcond_work->nbx_v2;
    // int *nbu_v2 =hpipm_pcond_work->nbu_v2;
    // int *ng_v2 = hpipm_pcond_work->ng_v2;
    // int *ns_v2 = hpipm_pcond_work->ns_v2;
    // int ii;
    // for (ii=0; ii<=N2; ii++)
	// 	nx_v2[ii] =hpipm_pcond_work->ocp_dim2->nx[ii];
	// for (ii=0; ii<=N2; ii++)
	// 	nu_v2[ii] = hpipm_pcond_work->ocp_dim2->nu[ii];
	// for (ii=0; ii<=N2; ii++)
	// 	nbx_v2[ii] = hpipm_pcond_work->ocp_dim2->nbx[ii];
	// for (ii=0; ii<=N2; ii++)
	// 	nbu_v2[ii] = hpipm_pcond_work->ocp_dim2->nbu[ii];
	// for (ii=0; ii<=N2; ii++)
	// 	ng_v2[ii] = hpipm_pcond_work->ocp_dim2->ng[ii];
	// for (ii=0; ii<=N2; ii++)
	// 	ns_v2[ii] = hpipm_pcond_work->ocp_dim2->ns[ii];
	// for(ii=0; ii<=N2; ii++)
	// 	printf("\n%d %d %d %d %d\n", nx_v2[ii], nu_v2[ii], nbx_v2[ii], nbu_v2[ii], ng_v2[ii]);


    Tic(&tcond);
    d_cond_qp_ocp2ocp(hpipm_pcond_work->ocp_qp, hpipm_pcond_work->ocp_qp2, hpipm_pcond_work->part_cond_arg, hpipm_pcond_work->part_cond_ws);
    d_cond_rhs_qp_ocp2ocp(hpipm_pcond_work->ocp_qp, hpipm_pcond_work->ocp_qp2, hpipm_pcond_work->part_cond_arg, hpipm_pcond_work->part_cond_ws);
    hpipm_pcond_work->tcond = Toc(&tcond);

    // printf("after condensing\n");

    Tic(&tqp);
    d_solve_ocp_qp_ipm(hpipm_pcond_work->ocp_qp2, hpipm_pcond_work->ocp_sol2,
        hpipm_pcond_work->ocp_arg, hpipm_pcond_work->ocp_workspace);
    hpipm_pcond_work->tqp = Toc(&tqp);

    d_expand_sol_ocp2ocp(hpipm_pcond_work->ocp_qp, hpipm_pcond_work->ocp_qp2, hpipm_pcond_work->ocp_sol2, hpipm_pcond_work->ocp_sol, hpipm_pcond_work->part_cond_arg, hpipm_pcond_work->part_cond_ws);

    d_cvt_ocp_qp_sol_to_colmaj(hpipm_pcond_work->ocp_sol, hu, hx, NULL, NULL, hpi, hlam_lb, hlam_ub, hlam_lg, hlam_ug, NULL, NULL);

    return 0;
}

void qpsolver_hpipm_pcond_workspace_free(model_size *size, qpsolver_hpipm_pcond_workspace* work)
{ 
    int N = size->N;

    free(work->nx_v);
    free(work->nu_v);
    free(work->nbx_v);
    free(work->nbu_v);
    free(work->ng_v);
    free(work->ns_v);

    free(work->nx_v2);
    free(work->nu_v2);
    free(work->nbx_v2);
    free(work->nbu_v2);
    free(work->ng_v2);
    free(work->ns_v2);
    free(work->block_size);

    free(work->ocp_dim);
    free(work->ocp_qp);
    free(work->ocp_sol);
    free(work->ocp_arg);
    free(work->ocp_workspace);

    free(work->ocp_dim2);
    free(work->ocp_qp2);
    free(work->ocp_sol2);
    free(work->part_cond_arg);
    free(work->part_cond_ws);

    free(work->dim_mem);
    free(work->qp_mem);
    free(work->sol_mem);
    free(work->arg_mem);
    free(work->ipm_mem);

    free(work->dim_mem2);
    free(work->qp_mem2);
    free(work->sol_mem2);
    free(work->part_cond_arg_mem);
    free(work->part_cond_mem);

    free(work->hA);
    free(work->hB);
    free(work->hb);
    free(work->hQ);
    free(work->hS);
    free(work->hR);
    free(work->hq);
    free(work->hr);
    free(work->hlb);
    free(work->hub);
    free(work->hC);
    free(work->hD);
    free(work->hlg);
    free(work->hug);

    int i;
    for(i=0;i<=N;i++){
        int_free(work->hidxb[i]);
        d_free(work->hlam_lb[i]);
        d_free(work->hlam_ub[i]);
        d_free(work->hlam_lg[i]);
        d_free(work->hlam_ug[i]);
    }
    free(work->hidxb);
    free(work->hlam_lb);
    free(work->hlam_ub);
    free(work->hlam_lg);
    free(work->hlam_ug);

    free(work->b0);
    free(work->r0);
    free(work->lg0);
    free(work->ug0);

    free(work->hu);
    free(work->hx);
    free(work->hpi);

    free(work->idxc);
    
    free(work);
}