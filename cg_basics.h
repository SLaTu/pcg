 
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

	unsigned int size, nnz;
	unsigned int *rows, *cols;
	double *values; 

} mat_t;

typedef struct pattern_matrix {

	unsigned int *rows;
	unsigned int *cols;
	unsigned int nnz;
	
} pat_t;

void multMatVect(double *tmp, double *vect, mat_t *mat);
void multMatVectCSC(double *tmp, double *vect, mat_t *mat);
void zerovect(double *vect, unsigned int max);
double multVectVect(double *V1, double *V2, unsigned int max);
void scaleVect(double scalar, double *V1, double *tmp, unsigned int max);
void addToVect(double *V1, double *V2, unsigned int max);
void subToVect(double *V1, double *V2, unsigned int max);
void addVects(double *V1, double *V2, double *Vres, unsigned int max);
void subtVects(double *V1, double *V2, double *result, unsigned int max);
double getvalue_mat(unsigned int row, unsigned int col, mat_t *mat);



void pmultMatVect(double *tmp, double *vect, mat_t *mat);
void pmultMatVectCSC(double *tmp, double *vect, mat_t *mat);
void pmultMatVect_DUMM(double *tmp, mat_t *mat);
void pmultMatVectCSC_DUMM(double *tmp, mat_t *mat);
void pzerovect(double *vect, unsigned int max);
double pmultVectVect(double *V1, double *V2, unsigned int max);
void pscaleVect(double scalar, double *V1, double *tmp, unsigned int max);
void paddToVect(double *V1, double *V2, unsigned int max);
void psubToVect(double *V1, double *V2, unsigned int max);
void paddVects(double *V1, double *V2, double *Vres, unsigned int max);
void psubtVects(double *V1, double *V2, double *result, unsigned int max);
void pequalvects(double *V1, double *V2, unsigned int max);

#endif // __CG_BASICS_H__ 
