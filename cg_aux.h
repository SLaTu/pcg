 
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
void report(mat_t *A, double *r, double *b, unsigned int i, double start_usec, double end_usec);
void reperror(mat_t *A, double *r, double *b, unsigned int i);
void precond(mat_t *mat, mat_t *G, double *xfinal);
double *cg_construct(mat_t *A, double *x, double *b, int maxiter, double error);
void transpose(mat_t *G, mat_t *Gtransp);
void transposeCSC(mat_t *G, mat_t *Gtransp);
mat_t *inverseD(mat_t *mat);

#endif // __CG_AUX_H__
