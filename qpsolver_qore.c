#include <stdlib.h>
#include <stdio.h>
#include "string.h"

#include "qpsolver_qore.h"
#include "full_condensing.h"
#include "common.h"
#include "qpsolver_dense.h"
#include "math.h"

qpsolver_qore_workspace* qpsolver_qore_workspace_create(model_size *size)
{
    // int nx=size->nx;
    int nu=size->nu;
    int nbg = size->nbg;
    int nbgN = size->nbgN;
    int N = size->N;

    qpsolver_qore_workspace *work = (qpsolver_qore_workspace*)malloc(sizeof(qpsolver_qore_workspace));

    work->lb_qore=(double*)malloc((N*nu+N*nbg+nbgN)*sizeof(double));
    work->ub_qore=(double*)malloc((N*nu+N*nbg+nbgN)*sizeof(double));
    work->sol_qore=(double*)malloc((N*nu+N*nbg+nbgN)*sizeof(double));
    
    work->problem= malloc(sizeof(struct QoreProblemDense*));
    QPDenseNew(&work->problem[0], N*nu, N*nbg+nbgN);
    QPDenseSetInt(work->problem[0], "prtfreq", -1);

    return work;
}

int qpsolver_qore(model_size *size, qp_problem *qp, full_condensing_workspace *full_condensing_work,
 qpsolver_qore_workspace *qore_work, qp_out *out, int iter)
{
    // int nx=size->nx;
    int nu=size->nu;
    int nbg = size->nbg;
    int nbgN = size->nbgN;
    int N = size->N;

    double *Hc = full_condensing_work->Hc;
    double *gc = full_condensing_work->gc;
    double *Cc = full_condensing_work->Cc;
    double *lcc = full_condensing_work->lcc;
    double *ucc = full_condensing_work->ucc;

    double *lb_qore = qore_work->lb_qore;
    double *ub_qore = qore_work->ub_qore;
    double *sol_qore = qore_work->sol_qore;
    QoreProblemDense **problem = qore_work->problem;

    memcpy(lb_qore, qp->lb_u, N*nu*sizeof(double));
    memcpy(lb_qore+N*nu, lcc, (N*nbg+nbgN)*sizeof(double));
    memcpy(ub_qore, qp->ub_u, N*nu*sizeof(double));
    memcpy(ub_qore+N*nu, ucc, (N*nbg+nbgN)*sizeof(double));

    if (iter==0){
        QPDenseSetData(problem[0], N*nu, N*nbg+nbgN, Cc, Hc);
        QPDenseOptimize(problem[0], lb_qore, ub_qore, gc, NULL, NULL);
    }else{
        //  err = QPDenseSetInt(problem[0], "warmstrategy", 1);
        QPDenseUpdateMatrices(problem[0], N*nu, N*nbg+nbgN, Cc, Hc);
        QPDenseOptimize(problem[0], lb_qore, ub_qore, gc, NULL, NULL);
    }   
    QPDenseGetDblVector(problem[0], "primalsol", sol_qore);

    memcpy(out->du, sol_qore, N*nu*sizeof(double));

    return 0;
}

void qpsolver_qore_workspace_free(qpsolver_qore_workspace *work)
{

    free(work->lb_qore);
    free(work->ub_qore);
    free(work->sol_qore);
    QPDenseFree(&work->problem[0]);
    free(work->problem);

    free(work);

}