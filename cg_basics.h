 
#ifndef __CG_BASICS_H__
#define __CG_BASICS_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#ifdef INTEL_MKL

#include "mkl.h"
#define BLAS_cp(n, dx, incx, dy, incy) 						cblas_dcopy(n, dx, incx, dy, incy)
#define BLAS_dot(n, dx, incx, dy, incy) 					cblas_ddot(n, dx, incx, dy, incy) 
#define BLAS_axpy(n, da, dx, incx, dy, incy) 				cblas_daxpy(n, da, dx, incx, dy, incy)
#define SBLAS_csrmv(trans, m, n, alpha, matdescra, avval, avpos, avptr, avptr1, Bptr, beta, Cptr) \
	mkl_dcsrmv(trans, &m, &n, &alpha, matdescra, avval, avpos, avptr, avptr1, Bptr, &beta, Cptr)

#elif defined LAPACK

#include "cblas.h"
#define BLAS_cp(n, dx, incx, dy, incy) 						cblas_dcopy(n, dx, incx, dy, incy)
#define BLAS_dot(n, dx, incx, dy, incy) 					cblas_ddot(n, dx, incx, dy, incy) 
#define BLAS_axpy(n, da, dx, incx, dy, incy) 				cblas_daxpy(n, da, dx, incx, dy, incy)
#define SBLAS_csrmv(trans, m, n, alpha, matdescra, avval, avpos, avptr, avptr1, Bptr, beta, Cptr) \
	manual_csrmv(trans, m, n, alpha, avval, avpos, avptr, Bptr, beta, Cptr)

#endif

typedef struct sparse_matrix {

	int size, nnz;
	int *rows, *cols;
	double *values; 

} mat_t;

typedef struct pattern_matrix {

	int *rows;
	int *cols;
	int nnz;
	
} pat_t;


static inline __attribute__((always_inline)) double rnum(){
	return (double)rand()/(double)RAND_MAX;
}

double getvalue_mat(int row, int col, mat_t *mat);
double getvalue_matCSC(int row, int col, mat_t *mat);

void pmultMatVect(double *tmp, double *vect, mat_t *mat);
void pmultMatVectCSC(double *tmp, double *vect, mat_t *mat);
void pmultMatVectCOO(double *tmp, double *vect, mat_t *mat, int *limits, int nthreads);
void pmultMatVect_DUMM(double *tmp, mat_t *mat);
void pmultMatVectCSC_DUMM(double *tmp, mat_t *mat);
void pmultMatVectCOO_DUMM(double *tmp, mat_t *mat, int *limits, int nthreads);
void pcopyvect(double *IN, double *OUT, int max);
void pzerovect(double *vect, int max);
double pmultVectVect(double *V1, double *V2, int max);
void pscaleVect(double scalar, double *V1, double *tmp, int max);
void paddToVect(double *V1, double *V2, int max);
void psubToVect(double *V1, double *V2, int max);
void paddVects(double *V1, double *V2, double *Vres, int max);
void psubtVects(double *V1, double *V2, double *result, int max);
void pequalvects(double *V1, double *V2, int max);
void cleararrays(double *r, double *rold, double *d, double *q, double *s, double *tmp, double *x, int dim);
void clearlogs(double *elapses, double *residuals, int imax);

#endif // __CG_BASICS_H__ 
