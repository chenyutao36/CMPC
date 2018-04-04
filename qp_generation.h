#ifndef QP_GENERATION_H_
#define QP_GENERATION_H_

#include "common.h"

typedef struct{
    double *x;
    double *u;
    double *y;
    double *yN;
    double *W;
    double *WN;
    double *p;
    double *lb;
    double *ub;
    double *lbg;
    double *ubg;
    double *lbgN;
    double *ubgN;
}qp_in;

typedef struct{
    double *Q;
    double *S;
    double *R;
    double *A;
    double *B;
    double *b;
    double *gx;
    double *gu;
    double *C;
    double *CN;
    double *D;
    double *DN;
    double *lb;
    double *ub;
    double *lb_g;
    double *ub_g;
}qp_problem;

typedef struct{
    double **Jac;
    double *Jac_N;
    int sample;
}qp_generation_workspace;

typedef struct{
    double *dx;
    double *du;
    double *lam;
    double *mu;
}qp_out;

qp_in* qp_in_create(model_size *size);
qp_problem* qp_problem_create(model_size *size);
qp_generation_workspace* qp_generation_workspace_create(model_size *size);
qp_out* qp_out_create(model_size *size);

int qp_generation(qp_in *in, model_size *size, 
    qp_problem *qp, qp_out *out, qp_generation_workspace *work);

int expand(model_size *size, qp_problem *qp, qp_out *out);

void qp_in_free(qp_in *in);
void qp_problem_free(qp_problem *qp);
void qp_generation_workspace_free(qp_generation_workspace *work);
void qp_out_free(qp_out *out);
#endif