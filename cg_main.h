#ifndef __CG_MAIN_H__
#define __CG_MAIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "papi.h"
#include "cg_basics.h"
#include "cg_aux.h"

static inline __attribute__((always_inline)) void dump_info(char *name, int k, double *residuals, double *elapse)
{
	FILE *log = fopen(name, "w");
	for ( int i = 0; i <= k; i++ ) {
		fprintf(log, "%d\t%.8e\t%.8e\n", i, residuals[i], elapse[i]);
	}
	fclose(log);
}

static inline __attribute__((always_inline)) void print_pattern(char *name, mat_t *G, int lim)
{
	FILE *log = fopen(name, "w");
	for ( int i = 0; i < lim; i++ ) {
		for(int j = G->rows[i]; j < G->rows[i + 1]; j++){
			fprintf(log, "%i %i %lf\n", i, G->cols[j], G->values[j]);
		}
	}
	fclose(log);
}

double *cg_base(mat_t *A, double *x, double *b, char *argv);  
double *cg_precond_diag(mat_t *A, double *x, double *b, char *argv);
double *cg_precond(mat_t *A, double *x, double *b, char *argv, double *initx);
void lu(mat_t *A, double *x, double *b, double *initx, double *initb);

#endif //__CG_MAIN_H__ 
