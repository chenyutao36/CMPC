#ifndef FULL_CONDENSING_H_
#define FULL_CONDENSING_H_

#include "common.h"
#include "qp_generation.h"

typedef struct{
    struct blasfeo_dmat *sStQ;
    struct blasfeo_dmat *sRt;
    struct blasfeo_dmat *sBtAt;
    struct blasfeo_dmat *sCt;
    
    struct blasfeo_dvec *srq;

    struct blasfeo_dmat *sGt;
    struct blasfeo_dmat *sHtWt;
    struct blasfeo_dmat *sWttmp;
    struct blasfeo_dvec *shw;
    struct blasfeo_dvec *sL;

    struct blasfeo_dmat *sHc;
    struct blasfeo_dmat *sCct;
    struct blasfeo_dvec *sgc;

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