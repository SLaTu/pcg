 
#ifndef __CG_PATTERNS_H__
#define __CG_PATTERNS_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h> 

#include "cg_basics.h"

void ltp(mat_t *mat, pat_t *pattern, pat_t *expanded_patt);
void flt(mat_t *mat, pat_t *pattern, pat_t *expanded_patt);
void powerA(mat_t *mat, pat_t *pattern, pat_t *expanded_patt, double *xfinal);
void perclt(mat_t *mat, pat_t *pattern, pat_t *expanded_patt, double *xfinal);


#endif // __CG_PATTERNS_H__
