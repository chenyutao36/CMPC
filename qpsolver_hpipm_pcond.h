#ifndef QPSOLVER_HPIPM_PCOND_H_
#define QPSOLVER_HPIPM_PCOND_H_

#include "qp_generation.h"
#include "common.h"

typedef struct{
    int *nx_v;
    int *nu_v;
    int *nbx_v;
    int *nbu_v;
    int *ng_v;
    int *ns_v;

    int *nx_v2;
    int *nu_v2;
    int *nbx_v2;
    int *nbu_v2;
    int *ng_v2;
    int *ns_v2;
    int *block_size;

    struct d_ocp_qp_dim *ocp_dim;
    struct d_ocp_qp_dim *ocp_dim2;
    struct d_ocp_qp *ocp_qp;
    struct d_ocp_qp *ocp_qp2;
    struct d_ocp_qp_sol *ocp_sol;
    struct d_ocp_qp_sol *ocp_sol2;
    struct d_cond_qp_ocp2ocp_arg *part_cond_arg;
    struct d_cond_qp_ocp2ocp_workspace *part_cond_ws;
    struct d_ocp_qp_ipm_arg *ocp_arg;
    struct d_ocp_qp_ipm_workspace *ocp_workspace;
    void *dim_mem;
    void *dim_mem2;
    void *qp_mem;
    void *qp_mem2;
    void *sol_mem;
    void *sol_mem2;
    void *part_cond_arg_mem;
    void *part_cond_mem;
    void *arg_mem;
    void *ipm_mem;

    double *b0;
    double *r0;
    double *lg0;
    double *ug0;

    double **hA;
    double **hB;
    double **hb;
    double **hQ;
    double **hS;
    double **hR;
    double **hq;
    double **hr;
    double **hlb;
    double **hub;
    double **hC;
    double **hD;
    double **hlg;
    double **hug;
    int **hidxb;

    double **hx;
    double **hu;
    double **hpi;
    double **hlam_lb;
    double **hlam_ub;
    double **hlam_lg;
    double **hlam_ug;

    int *idxc;

    double tcond;
    double tqp;

}qpsolver_hpipm_pcond_workspace;

qpsolver_hpipm_pcond_workspace* qpsolver_hpipm_pcond_workspace_create(model_size *size);

void qpsolver_hpipm_pcond_workspace_init(model_size *size,
    qp_problem *qp, qp_out *out, qpsolver_hpipm_pcond_workspace *hpipm_pcond_work);

int qpsolver_hpipm_pcond(model_size *size, qp_problem *qp, qp_out *out, qpsolver_hpipm_pcond_workspace *hpipm_pcond_work);

void qpsolver_hpipm_pcond_workspace_free(model_size *size, qpsolver_hpipm_pcond_workspace* work);

#endif