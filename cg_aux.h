 
#ifndef __CG_AUX_H__
#define __CG_AUX_H__

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sys/types.h>

#include "cg_basics.h"
#include "cg_patterns.h"


int cg_config(int argc, char *argv[]);  
int cg_setup(char *matname, mat_t *mat);
int readrhs(double *b, char* matname, unsigned int size);
void report(mat_t *A, double *r, double *b, unsigned int i, double start_usec, double end_usec);
void printmat(mat_t *A);
void printarr(double *arr, unsigned int dim);
void reperror(mat_t *A, double *r, double *b, unsigned int i);
void precond(mat_t *mat, mat_t *G, double *xfinal);
double *cg_construct(mat_t *A, double *x, double *b, int maxiter, double error);
void transpose(mat_t *G, mat_t *Gtransp);
void transposeCSC(mat_t *G, mat_t *Gtransp);
mat_t *inverseD(mat_t *mat);
void CSCtoCOO(mat_t *Gtransp, mat_t *GtranspCOO, unsigned int *limits, unsigned int nthreads);
void filterG(mat_t *G, mat_t *A, unsigned int dim, double *r);

void smoothmat(mat_t *A);
void smootharr(double *b, unsigned int dim);
int sumupdown(mat_t *A, int *col, int ract);

void smoothmatLLt(mat_t *A);
int sumLLt(mat_t *A, int *col, int ract);
void smoothmatL(mat_t *A);
int sumL(mat_t *A, int *col, int ract);
void smootharrL(double *b, unsigned int dim);
void solveLU(mat_t *A, unsigned int dim, double *b);
void impLU(mat_t *A, unsigned int dim, double *b);
void backsmootharrLt(double *b, unsigned int dim);

#endif // __CG_AUX_H__
