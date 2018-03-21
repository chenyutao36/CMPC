#ifndef FULL_CONDENSING_H_
#define FULL_CONDENSING_H_

#include "common.h"
#include "qp_generation.h"

typedef struct{
    double *G;
    double *W;
    double *w;
    double *L;
    double *Hi;
    double *Cci;
    double *CcN;

    double *Hc;
    double *gc;
    double *Cc;
    double *lcc;
    double *ucc;

}full_condensing_workspace;

full_condensing_workspace* full_condensing_workspace_create(model_size *size);

int full_condensing(double *dx0, model_size *size, qp_problem *qp, full_condensing_workspace *work);

void full_condensing_workspace_free(full_condensing_workspace* work);

#endif