#ifndef QPSOLVER_QORE_H_
#define QPSOLVER_QORE_H_

#include "common.h"
#include "qp_generation.h"
#include "full_condensing.h"

typedef struct{
    struct QoreProblemDense **problem;
    double *lb_qore;
    double *ub_qore;
    double *sol_qore;
}qpsolver_qore_workspace;

qpsolver_qore_workspace* qpsolver_qore_workspace_create(model_size *size);

int qpsolver_qore(model_size *size, qp_problem *qp, full_condensing_workspace *full_condensing_work,
 qpsolver_qore_workspace *qore_work, int iter);

void qpsolver_qore_workspace_free(qpsolver_qore_workspace *work);

#endif