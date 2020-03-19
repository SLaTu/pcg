
#include "cg_main.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define PI 3.14159265
#define meannums 1

char matname[256];
int imax;
int mode;
double err;
FILE *fp;
FILE *outmn;
int patternmode;
int percentpattern;
int patternpower;
int bsize;
int reps;
int rhs;
mat_t *powmat;

/*
 * 
 * 
 * 				MAIN
 * 
 * 
 */

int main(int argc, char *argv[]){
	
	mat_t *A;
	A = malloc(sizeof(mat_t));
	
	if (cg_config(argc, argv)){
		return 1;
	}
	if (cg_setup(matname, A)){
		return 2;
	}
	
	int dim = A->size;
	double *b = calloc(dim, sizeof(double));
	double *x = calloc(dim, sizeof(double));
	double *initx = malloc(dim*sizeof(double));
	double *initb = malloc(dim*sizeof(double));
	int j;
	
	if (rhs == 1){
		readrhs(b, argv[1], dim);
		for (j = 0; j < dim; j++) {
			initb[j] = b[j];
		}
	}
	else {
		for (j = 0; j < dim; j++) {
			x[j] = rnum();
			initx[j] = x[j];
		}
		pmultMatVect(b, x, A);
		for (j = 0; j < dim; j++) {
			initb[j] = b[j];
			x[j] = 0.0; 
		}
	}
	
	printf("Init X\n");
	printarr(initx, dim);
	
	printf("B\n");
	printarr(b, dim);
	
	printf("NNZ:\t%u\n", A->nnz);
	printf("Size:\t%u\n\n", dim);
	
// 	char buf[256];
// 	snprintf(buf, sizeof buf, "../Outputs/cg/DataCG/DATA_%s_%.1e", matname, err);
// 	
// 	
// 	outmn = fopen(buf, "a");
// 	if (outmn == NULL){
// 		printf("Could not open writing file.");
// 		return 0;
// 	}
	
	switch ( mode ) {
		case 1: 
			printf("Mode:\tBase CG\n");
			x = cg_base(A, x, b, argv[1]);
			break;
		case 2:
			printf("Mode:\tD⁻1 Preconditioned CG\n");
			fprintf(stderr, "Mode:\tD⁻1 Preconditioned CG\n");
			x = cg_precond_diag(A, x, b, argv[1]);
			break;
		case 3: 
			printf("Mode:\tPCG\n");
			fprintf(stderr, "Mode:\tPCG");
			x = cg_precond(A, x, b, argv[1], initx);
			break;
		case 4: 
			printf("Mode:\tLU\n");
			lu(A, x, b, initx, initb);
			break;
		default:
			printf("No algorithm selected.\n");
			break;
	}
	
	return 0;
}
	
/*
 * 
 * 
 * 				BASE CG
 * 
 * 
 */
	
double *cg_base(mat_t *A, double *x, double *b, char *argv){
	
	int i = 0;
	double d_new = 0.0;
	double alpha = 0.0;
	double d_old = 0.0;
	double beta = 0.0;
	double start_usec = 0.0, end_usec = 0.0;
	double start = 0.0, end = 0.0;
	double norm_r = 0.0;
	int dim = A->size;
	
	double *r = malloc(dim*sizeof(double));
	double *d = malloc(dim*sizeof(double));
	double *q = malloc(dim*sizeof(double));
	double *tmp = malloc(dim*sizeof(double));
	double *residuals = malloc(imax * sizeof(double));
	double *elapses = malloc(imax * sizeof(double));
	
	pmultMatVect(tmp, x, A);
	psubtVects(b, tmp, r, dim);
	pequalvects(d, r, dim);
	d_new = pmultVectVect(r, r, dim);
	double norm_b = sqrt(pmultVectVect(b, b, dim));
	
	start_usec = omp_get_wtime();
	for (i = 0; i < imax; i++){
		start = omp_get_wtime();
		pmultMatVect(q, d, A);
		alpha = d_new / pmultVectVect(d, q, dim);
		pscaleVect(alpha, d, tmp, dim);
		paddToVect(x, tmp, dim);
		if (i%50 == 0){
			pmultMatVect(tmp, x, A);
			psubtVects(b, tmp, r, dim);
		}
		else {
			pscaleVect(alpha, q, tmp, dim);
			psubToVect(r, tmp, dim);
		}
		d_old = d_new;
		d_new = pmultVectVect(r, r, dim);
		beta = d_new/d_old;
		pscaleVect(beta, d, tmp, dim);
		paddVects(r, tmp, d, dim);
		end = omp_get_wtime();
		elapses[i] = end - start;
		norm_r = sqrt(pmultVectVect(r, r, dim));
		residuals[i] = norm_r/norm_b;
		if (isless(residuals[i], err)){
			end_usec = omp_get_wtime();
			fprintf(stdout, "Precision reached; Iter: %i; Error: %.2e; Total time: %lf\n\n", i, residuals[i], end_usec - start_usec);
			fprintf(stderr, "\nPrecision reached; Iter: %i; Error: %.2e; Total time: %lf\n\n", i, residuals[i], end_usec - start_usec);
			break;}
	}
	if (i == imax) printf("No convergence. Residual = %.4e\n\n", residuals[i - 1]);
	
	char buf_t[256];
	snprintf(buf_t, sizeof buf_t, "../Outputs/cg/DataCG/logs/basecg_%s.log", argv);
	dump_info(buf_t, i, residuals, elapses);
	
	free(r); free(d); free(q); free(tmp); free(residuals); free(elapses);
	return x;
}

/*
 * 
 * 
 * 				DPCG - Preconditioned CG with D⁻1
 * 
 * 
 */

double *cg_precond_diag(mat_t *A, double *x, double *b, char *argv){

	int i = 0;
	double d_new = 0.0; 
	double d_old = 0.0;
	double alpha = 0.0;
	double beta = 0.0;
	double start_usec = 0.0, end_usec = 0.0;
	double start = 0.0, end = 0.0;
	int dim = A->size;
	
	double *r = calloc(dim, sizeof(double));
	double *d = calloc(dim, sizeof(double));
	double *q = calloc(dim, sizeof(double));
	double *s = calloc(dim, sizeof(double));
	double *tmp = calloc(dim, sizeof(double));
	double *residuals = malloc(imax * sizeof(double));
	double *elapses = malloc(imax * sizeof(double));
	
	mat_t *G;
	G = malloc(sizeof(mat_t));
	G->values = malloc(dim*sizeof(double));
	G->cols = malloc(dim*sizeof(int));
	G->rows = malloc((dim + 1)*sizeof(int));
	
	G = inverseD(A);									/* Precond Inverse Diagonal */
	
	pmultMatVect(tmp, x, A);
	psubtVects(b, tmp, r, dim);
	pmultMatVect(d, r, G);
	d_new = pmultVectVect(d, r, dim);
	double norm_b = sqrt(pmultVectVect(b, b, dim));
	double norm_r = sqrt(pmultVectVect(r, r, dim));
	start_usec = omp_get_wtime();
	for (i = 0; i < imax; i++){
		start = omp_get_wtime();
		pmultMatVect(q, d, A);
		alpha = d_new / pmultVectVect(d, q, dim);
		pscaleVect(alpha, d, tmp, dim);
		paddToVect(x, tmp, dim);
		if (i%50 == 0){
			pmultMatVect(tmp, x, A);
			psubtVects(b, tmp, r, dim);
		}
		else {
			pscaleVect(alpha, q, tmp, dim);
			psubToVect(r, tmp, dim);
		}
		pmultMatVect(s, r, G);
		d_old = d_new;
		d_new = pmultVectVect(r, s, dim);
		beta = d_new/d_old;
		pscaleVect(beta, d, tmp, dim);
		paddVects(s, tmp, d, dim);
		end = omp_get_wtime();
		elapses[i] = end - start;
		norm_r = sqrt(pmultVectVect(r, r, dim));
		residuals[i] = norm_r/norm_b;
			
		if ( isless(residuals[i], err) ) {
			end_usec = omp_get_wtime();
			fprintf(stdout, "Precision reached; Iter: %i; Error: %.2e; Total time: %lf\n\n", i, residuals[i], end_usec - start_usec);
			fprintf(stderr, "\nPrecision reached; Iter: %i; Error: %.2e; Total time: %lf\n\n", i, residuals[i], end_usec - start_usec);
			break;
		}
	}
	if (i == (imax)) printf("No convergence. Residual = %.4e\n\n", residuals[i-1]);
	
	char buf_t[256];
	snprintf(buf_t, sizeof buf_t, "../Outputs/cg/DataCG/logs/diagpcg_%s.log", argv);
	dump_info(buf_t, i, residuals, elapses);
	
	free(d); free(r); free(q); free(s); free(G); free(tmp);
	
	return x;
}

/*
 * 
 * 
 * 				PCG - Includes powering and pattern extending
 * 
 * 
 */

double *cg_precond(mat_t *A, double *x, double *b, char *argv, double *initx){

	int i = 0, m = 0;
	int dim = A->size;
	char buf_t[256];
	double d_new = 0.0, alpha = 0.0, d_old = 0.0; 
	double start_usec = 0.0, end_usec = 0.0, start = 0.0, end = 0.0, ttmp = 0.0, tmult = 0.0, ttransp = 0.0;
	double tmin = 1000.0, btottime = 0.0, btimeiter = 0.0, btimemult = 0.0, btimetransp = 0.0;
	int bniter = 0, brep = 0;;
	long long bdcmnorm = 0, bxdcmnorm = 0, bdifmult = 0, bdcmtransp = 0, bxdcmtransp = 0, bdiftransp = 0;
	long long dcmnorm = 0, dcmtransp = 0, xdcmnorm = 0, xdcmtransp = 0;
	int retval = 0, EventSet = PAPI_NULL;
	long_long values[1];
	
	omp_sched_t kind;
	int chunk;
	omp_get_schedule(&kind, &chunk);
	
	mat_t *G, *Gtransp;
	
	G = malloc(sizeof(mat_t));
	Gtransp = malloc(sizeof(mat_t));
	G->size = dim;
	Gtransp->size = dim;
	
	
	double *s = calloc(dim, sizeof(double));
	double *r = calloc(dim, sizeof(double));
	double *rold = calloc(dim, sizeof(double));
	double *d = calloc(dim, sizeof(double));
	double *q = calloc(dim, sizeof(double));
	double *tmp = calloc(dim, sizeof(double));
	double *dumm = calloc(dim, sizeof(double));
	double *residuals = malloc(imax * sizeof(double));
	double *elapses = malloc(imax * sizeof(double));
	
	
// 	poweredA(A, A->size);
	start_usec = omp_get_wtime();											/* PRECONDITIONER + TRANSPOSED */
	precond(A, G, r);
// 	precond_vf(A, G, r);
	end_usec = omp_get_wtime();
	
	printf("\nPreconditioning Time: %lf\n", end_usec - start_usec);
	
	transposeCSC(G, Gtransp);
	start_usec = omp_get_wtime();
	printf("Transposing Time: %lf\n", start_usec - end_usec);
	
// 	char buf_p[256];
// 	snprintf(buf_p, sizeof buf_p, "../Outputs/cg/patterns/pattern_%s_%i_%i.mtx", argv, percentpattern, patternpower);
// 	print_pattern(buf_p, G, G->size/10);
	
	mat_t *GtranspCOO = malloc(sizeof(mat_t));
	GtranspCOO->size = G->size;
	GtranspCOO->nnz = G->nnz;
	GtranspCOO->values = malloc(G->nnz*sizeof(double));
	GtranspCOO->cols = malloc(G->nnz*sizeof(int));
	GtranspCOO->rows = malloc(G->nnz*sizeof(int));
	
	int nthreads = omp_get_max_threads(); 
	int *limits = calloc(nthreads + 1, sizeof(int));
	
	CSCtoCOO(Gtransp, GtranspCOO, limits, nthreads);
	
	free(Gtransp);
	
	for (m = 0; m < reps; m++){														/* REPEAT LOOP reps TIMES */
		
		cleararrays(r, rold, d, q, s, tmp, x, dim);
		clearlogs(elapses, residuals, imax);
		tmult = 0.0; ttransp = 0.0;
		
		pmultMatVect(tmp, x, A);														/* r = b - Ax */
		psubtVects(b, tmp, r, A->size);
		pmultMatVect(tmp, r, G);														/* d = M⁻1*r */
		pmultMatVectCOO(d, tmp, GtranspCOO, limits, nthreads);
		
		papiinit(retval, EventSet);
		papicount(dim, G, GtranspCOO, tmp, dcmnorm, xdcmnorm, dcmtransp, xdcmtransp, values, dumm, limits, nthreads, EventSet, r, s);
		papiend(EventSet, values);
		
		d_new = pmultVectVect(d, r, A->size);						/* d_new = d * r */
		double norm_b = sqrt(pmultVectVect(b, b, A->size));
		double norm_r = sqrt(pmultVectVect(r, r, A->size));
		start_usec = omp_get_wtime();
		for (i = 0; i < imax; i++){									/* BEGIN CG ITERATION */
			start = omp_get_wtime();
			pmultMatVect(q, d, A);									/* q = A * d */
			alpha = d_new / pmultVectVect(d, q, A->size);			/* alpha = d_new / (d * q) */
			pscaleVect(alpha, d, tmp, A->size);						/* x += alpha * d */
			paddToVect(x, tmp, A->size);
			if (i%50 == 0){											/* Residual correction */
				pmultMatVect(tmp, x, A);							/* r = b - Ax */
				psubtVects(b, tmp, r, A->size);
			}
			else {
				pscaleVect(alpha, q, tmp, A->size);					/* r -= alpha * q */
				psubToVect(r, tmp, A->size);
			}
			pequalvects(rold, r, A->size);
			ttmp = omp_get_wtime();
			pmultMatVect(tmp, r, G);								/* s = M⁻1 * r    --->   s = G * G^t * r */
			tmult += omp_get_wtime() - ttmp;
			pzerovect(r, A->size);
			ttmp = omp_get_wtime();
			pmultMatVectCOO(r, tmp, GtranspCOO, limits, nthreads);
			ttransp += omp_get_wtime() - ttmp;
			pequalvects(s, r, A->size);
			pequalvects(r, rold, A->size);
			d_old = d_new;											/* d_old = d_new */
			d_new = pmultVectVect(s, r, A->size);					/* d_new = r * s */
			pscaleVect(d_new/d_old, d, tmp, A->size);				/* d = s + (d_new/d_old) * d */
			paddVects(s, tmp, d, A->size);
																	/* END CG ITERATION */
			end = omp_get_wtime();
			elapses[i] = end - start;
			norm_r = sqrt(pmultVectVect(r, r, dim));
			residuals[i] = norm_r/norm_b;
			if (isless(residuals[i], err)) {
				end_usec = omp_get_wtime();
				fprintf(stdout, "Precision reached; Iter: %i; Error: %.2e; Total time: %lf\n\n", i, residuals[i], end_usec - start_usec);
				fprintf(stderr, "Precision reached; Iter: %i; Error: %.2e; Total time: %lf\n\n", i, residuals[i], end_usec - start_usec);
				printf("RESULTING X\n\n");
				printarr(x, A->size);
				break;
			}
		}
		if (i == (imax)) printf("No convergence. Residual = %.4e\n\n", residuals[i - 1]);
		
		if ((end_usec - start_usec) < tmin){
			tmin = end_usec - start_usec;
			bniter = i;
			btottime = end_usec - start_usec;
			btimeiter = tmin / (double) bniter;
			btimemult = tmult / (double) bniter;
			btimetransp = ttransp / (double) bniter;
			bdcmnorm = dcmnorm;
			bxdcmnorm = xdcmnorm;
			bdifmult = dcmnorm - xdcmnorm;
			bdcmtransp = dcmtransp;
			bxdcmtransp = xdcmtransp;
			bdiftransp = dcmtransp - xdcmtransp;
			brep = m;
		}
		
		snprintf(buf_t, sizeof buf_t, "../Outputs/cg/DataCG/logs/pcg_%s_%i_%i_%i_%i.log", argv, percentpattern, patternpower, m, rhs);
		dump_info(buf_t, i, residuals, elapses);
		
	}
	
	printfigdatax(matname, patternpower, rhs, G, A, percentpattern, bniter, btottime, btimeiter, btimemult, btimetransp, bdcmnorm, bxdcmnorm, bdifmult, bdcmtransp, bxdcmtransp, bdiftransp, brep, x, initx);
	
	free(dumm); free(r); free(d); free(q); free(s); free(tmp); free(G); free(residuals); free(elapses);
	
	return x;
}



void lu(mat_t *A, double *x, double *b, double *initx, double *initb){
	
	int dim = A->size;
	
	
	printmat(A);
	smoothmatLLt(A);
	
	
	printf("\n\nREAL X\n");
	printarr(initx, dim);
	printf("B\n");
	printarr(b, dim);
	
	smootharrL(b, dim);
	smootharrL(initb, dim);
	
	impLU(A, dim, b);
	
	printf("SMOOTH B\n");
	printarr(b, dim);
	
	solveLU(A, dim, b);
	
	printf("LU X (Z)\n");
	printarr(b, dim);
	
	
	double *bax = calloc(dim, sizeof(double));
	double *tmp = calloc(dim, sizeof(double));
	pmultMatVect(tmp, b, A);
	backsmootharrLt(b, dim);
	
	printf("X\n");
	printarr(b, dim);
	
	
	double norm_bax = 0.0;
	for (int i = 0; i < dim; i++){
		bax[i] = initb[i] - tmp[i];
		norm_bax += (initb[i] - tmp[i])*(initb[i] - tmp[i]);
	}
	
	printf("BAX\n");
	printarr(bax, dim);
	
	printf("NORM BAX:\t%.2e\n\n", sqrt(norm_bax));
	
	char buf_x[256];
	snprintf(buf_x, sizeof buf_x, "../Outputs/cg/DataCG/x_%s.log", matname);
	FILE *logt = fopen(buf_x, "w");
	for ( int i = 0; i < dim; i++ ) {
		fprintf(logt, "%i\t%.8e\t%.8e\n", i, b[i], initx[i]);
	}
	fclose(logt);
	
	
	
}




















