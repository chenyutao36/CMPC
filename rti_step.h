#ifndef RTI_STEP_H_
#define RTI_STEP_H_

#include "qpsolver_hpipm_ocp.h"
#include "qpsolver_hpipm_pcond.h"
#include "qp_generation.h"
#include "common.h"
#include "full_condensing.h"
#include "qpsolver_qore.h"

typedef struct{
    qp_in *in;
    qp_problem *qp;
    qp_generation_workspace *qp_work;  
    qp_out *out;
    qpsolver_hpipm_ocp_workspace *hpipm_ocp_work;
    qpsolver_hpipm_pcond_workspace *hpipm_pcond_work;      

    double *dx0;
    double *u_opt;
    int sample;

    double cpt_prep;
    double cpt_cond;
    double cpt_ocp_qp;
}rti_step_workspace;

rti_step_workspace* rti_step_workspace_create(model_size *size);

void rti_step_init(model_size *size, rti_step_workspace *rti_work);

int rti_step(double *x0, model_size *size, rti_opt *opt, rti_step_workspace *rti_work);

void rti_step_workspace_free(model_size *size, rti_step_workspace *work);

#endif