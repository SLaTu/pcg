 
#ifndef __CG_PATTERNS_H__
#define __CG_PATTERNS_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h> 

#include "cg_basics.h"

int compd (const void *a, const void *b);
int comp (const void *a, const void *b);
void ltp(mat_t *mat, pat_t *pattern, pat_t *expanded_patt);
void flt(mat_t *mat, pat_t *pattern, pat_t *expanded_patt);
void perclt(mat_t *mat, pat_t *pattern, pat_t *expanded_patt, double *xfinal);
void powerA(mat_t *mat, pat_t *pattern, pat_t *expanded_patt, double *xfinal);
void powerAf(mat_t *mat, pat_t *pattern, pat_t *expanded_patt, double *xfinal);
void OptpowerA(mat_t *mat, pat_t *pattern, pat_t *expanded_patt, double *xfinal);
void poweredA(mat_t *mat, unsigned int dim);


#endif // __CG_PATTERNS_H__
