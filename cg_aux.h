 #ifndef __CG_AUX_H__
#define __CG_AUX_H__

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sys/types.h>
#include "papi.h"

#include "cg_basics.h"
#include "cg_patterns.h"

int cg_config(int argc, char *argv[]);  
int cg_setup(char *matname, mat_t *mat);
int readrhs(double *b, char* matname, int size);

void printmat(mat_t *A);
void printarr(double *arr, int dim);

void papiinit(int retval, int EventSet);
void papiend(int EventSet, long_long *values);
void papicount(int dim, mat_t *G, mat_t *GtranspCOO, double *tmp, long long dcmnorm, long long xdcmnorm, long long dcmtransp, long long xdcmtransp, long_long *values, double *dumm, int *limits, int nthreads, int EventSet, double *r, double *s);

void precond(mat_t *mat, mat_t *G, double *xfinal);
void precond_vf(mat_t *mat, mat_t *G, double *xfinal);
mat_t *inverseD(mat_t *mat);

double *cg_construct(mat_t *A, double *x, double *b, int maxiter, double error);

void transpose(mat_t *G, mat_t *Gtransp);
void transposeCSC(mat_t *G, mat_t *Gtransp);
void CSCtoCOO(mat_t *Gtransp, mat_t *GtranspCOO, int *limits, int nthreads);

void filterG(mat_t *G, mat_t *A, int dim, double *r);

void smoothmat(mat_t *A);
void smootharr(double *b, int dim);
int sumupdown(mat_t *A, int *col, int ract);
void smoothmatLLt(mat_t *A);
void smoothTOmatLLt(mat_t *A, mat_t *Af);
int sumLLt(mat_t *A, int *col, int ract);
void smoothmatL(mat_t *A);
void smoothTOmatL(mat_t *A, mat_t *Af);
int sumL(mat_t *A, int *col, int ract);
void smootharrL(double *b, int dim);
void backsmootharrLt(double *b, int dim);
void backsmootharrLtV2(double *b, int dim, int *arrayindex);

void solveLU(mat_t *A, int dim, double *b);
void impLU(mat_t *A, int dim, double *b);

void printfigdatax(char *matname, int patternpower, int rhs, mat_t *G, mat_t *A, int percentpattern, int bniter, double btottime, double btimeiter, double btimemult, double btimetransp, long long bdcmnorm, long long bxdcmnorm, long long bdifmult, long long bdcmtransp, long long bxdcmtransp, long long bdiftransp, int m, double *x, double *initx);

void checkbesttime(double end_usec, double start_usec, double tmin, int bniter, int i, double btottime, double btimeiter, double btimemult, double btimetransp, long long bdcmnorm, long long bxdcmnorm, long long bdifmult, long long bdcmtransp, long long bxdcmtransp, long long bdiftransp, int brep, int m, double tmult, double ttransp, long long dcmnorm, long long xdcmnorm, long long dcmtransp, long long xdcmtransp);

#endif // __CG_AUX_H__
