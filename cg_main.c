
#include "cg_main.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define PI 3.14159265
#define meannums 1

char matname[256];
unsigned int imax;
int mode;
double err;
FILE *fp;
FILE *outmn;
unsigned int patternmode;
unsigned int percentpattern;
unsigned int patternpower;
unsigned int bsize;
unsigned int reps;
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
	
	unsigned int j;
	mat_t *A;
	A = malloc(sizeof(mat_t));
	
	
	if (cg_config(argc, argv)){
		return 1;
	}
	if (cg_setup(matname, A)){
		return 2;
	}
	double *b = calloc(A->size, sizeof(double));
	double *x = calloc(A->size, sizeof(double));
	double *initx = malloc(A->size*sizeof(double));
	double *initb = malloc(A->size*sizeof(double));
	unsigned int dim = A->size;
	
	fprintf(stdout, "RHS:\t%i\n", rhs);
	
	
	
	
	
// 	if (rhs == 1){
// 		readrhs(b, argv[1], A->size);
// 	}
// 	else {
		for (j = 0; j < A->size; j++) {
			x[j] = sin(((double) j) * PI/180000) + sin(((double) j) * PI/18000) + sin(((double) j) * PI/1800) + sin(((double) j) * PI/180);
			initx[j] = x[j];
		}
// 		x[A->size - 1] = 1.0;
// 		initx[A->size - 1] = 1.0;
		pmultMatVect(b, x, A);
		for (j = 0; j < A->size; j++) {
			initb[j] = b[j];
			x[j] = 0.0; 
		}
// 	}
	
	
	printf("Init X\n");
	printarr(initx, dim);
	
	printf("B\n");
	printarr(b, dim);
	
	
// if (rhs==1){
// 	printmat(A);
// 	smoothmatLLt(A);
// 	fprintf(stderr, "FILTERED -> NNZ: %i\tLT: %i\n\n", A->nnz, (A->nnz - A->size)/2 + A->size);
// 	
// 	smootharrL(b, dim);
// 	
// 	printf("Smooth B\n");
// 	printarr(b, dim);
// }
// else {
// 	fprintf(stderr, "NOT FILTERED\n\n");
// }
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	printf("NNZ:\t%u\n", A->nnz);
	printf("Size:\t%u\n\n", A->size);
	
	
	char buf[256];
	snprintf(buf, sizeof buf, "../Outputs/cg/DataCG/DATA_%s_%.1e", matname, err);
	
	
	outmn = fopen(buf, "a");
	if (outmn == NULL){
		printf("Could not open writing file.");
		return 0;
	}
	
	
	
	switch ( mode ) {
		case 1: 
			printf("Mode:\tBase CG\n");
			x = cg_base(A, x, b);
			break;
		case 2:
			printf("Mode:\tD⁻1 Preconditioned CG\n");
			fprintf(stderr, "Mode:\tD⁻1 Preconditioned CG");
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
	
	
// if (rhs==1){
// 	backsmootharrLt(x, dim);
// 	
// 	printf("Final X\n");
// 	printarr(x, dim);
// }
	
	
// 	fclose(outmn);
// 	free(A);
// 	free(b);
// 	free(x);
	return 0;
}

/*
 * 
 * 
 * 				BASE CG
 * 
 * 
 */

double *cg_base(mat_t *A, double *x, double *b){

	int i = 0;
	double *d;
	double d_new = 0.0;
	double *q;
	double *r;
	double alpha = 0.0;
	double d_old = 0.0;
	double beta = 0.0;
	double start_usec = 0.0, end_usec = 0.0;
	double *tmp;
	
	r = malloc(A->size*sizeof(double));
	d = malloc(A->size*sizeof(double));
	q = malloc(A->size*sizeof(double));
	tmp = malloc(A->size*sizeof(double));
	
	// 2
	multMatVect(tmp, x, A);
	subtVects(b, tmp, r, A->size);
	
	// 3
	for (i = 0; i < A->size; i++) d[i] = r[i];
	//d = r;
	// 4
	d_new = multVectVect(r, r, A->size);
	double norm_b = multVectVect(b, b, A->size);
	
	// 6
	i = 0;
	start_usec = omp_get_wtime();
	while ((i < imax)&&(sqrt(d_new)/norm_b > err)){
		
		// 7
		multMatVect(q, d, A);
		
		// 8
		alpha = d_new / multVectVect(d, q, A->size);
		
		// 9
		scaleVect(alpha, d, tmp, A->size);
		addToVect(x, tmp, A->size);
		
		// 10
		if (i%50 == 0){
			multMatVect(tmp, x, A);
			subtVects(b, tmp, r, A->size);
		}
		else {
			scaleVect(alpha, q, tmp, A->size);
			subToVect(r, tmp, A->size);
		}
		
		// 14
		d_old = d_new;
		
		// 15
		d_new = multVectVect(r, r, A->size);
		
		// 16
		beta = d_new/d_old;
		
		// 17
		scaleVect(beta, d, tmp, A->size);
		addVects(r, tmp, d, A->size);
		
		// 18
		
		if (i%200==0){
			printf("Iter: %i\tError: %.2e\n", i, sqrt(d_new)/norm_b);
		}
		
		
		i++;
	}
	end_usec = omp_get_wtime();
	
	report(A, r, b, i, start_usec, end_usec);

        fprintf(outmn, "%u\t", i);
        fprintf(outmn, "%lf\t", end_usec - start_usec);
        fprintf(outmn, "%lf\n", (end_usec - start_usec)/((double) i));
	
	free(r); free(d); free(q); free(tmp);
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

	int i = 0, j = 0;
	double *d;
	double *r;
	double d_new = 0.0; 
	double d_old = 0.0;
	double *q;
	double *s;
	double alpha = 0.0;
	double beta = 0.0;
	double start_usec = 0.0, end_usec = 0.0;
	double *tmp;
	double sumelapses = 0.0;
	
	double *residuals = malloc(imax * sizeof(double));
	double *elapses = malloc(imax * sizeof(double));
	double *dnew = malloc(imax * sizeof(double));
	
// 	printmat(A);
	
	mat_t *D_1;
	
	D_1 = malloc(sizeof(mat_t));
	
	D_1->values = malloc(A->size*sizeof(double));
	D_1->cols = malloc(A->size*sizeof(unsigned int));
	D_1->rows = malloc((A->size + 1)*sizeof(unsigned int));
	
	r = calloc(A->size, sizeof(double));
	d = calloc(A->size, sizeof(double));
	q = calloc(A->size, sizeof(double));
	s = calloc(A->size, sizeof(double));
	tmp = calloc(A->size, sizeof(double));
	
	// Inverse diagonal A
	D_1 = inverseD(A);

	// r = b - A * x
	pmultMatVect(tmp, x, A);

	psubtVects(b, tmp, r, A->size);

	// d = M^-1 * r      // MODIFY
	
	pmultMatVect(d, r, D_1); // instead of A D^-1

	// dnew = r^T * d
	d_new = multVectVect(d, r, A->size);
	double norm_b = sqrt(multVectVect(b, b, A->size));
	double norm_r = sqrt(multVectVect(r, r, A->size));
	for (i = 0; i < imax; i++){
		start_usec = omp_get_wtime();											/* BEGIN CG ITERATION */
		// q = A * d
		pmultMatVect(q, d, A);
		
		// alpha = dnew / (d^T * q)
		alpha = d_new / pmultVectVect(d, q, A->size);
		
		// x = x + alpha * d
		pscaleVect(alpha, d, tmp, A->size);
		paddToVect(x, tmp, A->size);
		
		// 
		if (i%50 == 0){
			// r = b - A * x
			pmultMatVect(tmp, x, A);
			psubtVects(b, tmp, r, A->size);
		}
		else {
			// r = r - alpha * q
			pscaleVect(alpha, q, tmp, A->size);
			psubToVect(r, tmp, A->size);
		}
		// PRECONDITIONING STEP -- s = M-1 * r -- s = G^t * G * r
		// G = lower triangular part of matrix A^L 
		
		pmultMatVect(s, r, D_1); // instead of A M^-1
		// dold = dnew
		d_old = d_new;
		
		// dnew = r^T * s
		d_new = pmultVectVect(r, s, A->size);
		dnew[i] = d_new;
		// beta = dnew / dold
		beta = d_new/d_old;
		
		// d = s + beta * d
		pscaleVect(beta, d, tmp, A->size);
		paddVects(s, tmp, d, A->size);
		
		end_usec = omp_get_wtime();
		elapses[i] = end_usec - start_usec;
		norm_r = sqrt(cblas_ddot(A->size, r, 1, r, 1));
		residuals[i] = norm_r/norm_b;
			
		if ( isless(residuals[i], err) ) {
			end_usec = 0.0;
			for (int j = 0; j < i; j++) end_usec += elapses[j];
			fprintf(stdout, "Precision reached; Iter: %i; Error: %.2e; Total time: %lf\n\n", i, residuals[i], end_usec);
			fprintf(stderr, "Precision reached; Iter: %i; Error: %.2e; Total time: %lf\n\n", i, residuals[i], end_usec);
			break;
		}
	}
	if (i == (imax)) printf("No convergence. Residual = %.4e\n\n", residuals[i-1]);
	
	for (j = 0; j <= i; j++) sumelapses += elapses[j]; 
	
	
	fprintf(outmn, "%u\t", i);
	fprintf(outmn, "%lf\t", sumelapses);
	fprintf(outmn, "%lf\n", sumelapses/((double) i));
	
	char buf_t[256];
	snprintf(buf_t, sizeof buf_t, "../Outputs/cg/DataCG/logs/diagpcg_%s.log", argv);
	dump_info(buf_t, i, residuals, elapses, dnew);
	
	
		
	free(d); free(r); free(q); free(s); free(D_1); free(tmp);
	
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

	unsigned int i = 0;
	int m = 0;
	double d_new = 0.0; 
	double alpha = 0.0;
	double d_old = 0.0;
	double start_usec = 0.0, end_usec = 0.0;
	
	int *sumiter = calloc(reps, sizeof(int));
	double *sumtottime = calloc(reps, sizeof(double));
	double *sumtimeiter = calloc(reps, sizeof(double));
	double *sumtmult = calloc(reps, sizeof(double));
	double *sumttransp = calloc(reps, sizeof(double));
	long long *sumdcmnorm = calloc(reps, sizeof(long long));
	long long *sumxdcmnorm = calloc(reps, sizeof(long long));
	long long *sumdifmult = calloc(reps, sizeof(long long));
	long long *sumdcmtransp = calloc(reps, sizeof(long long));
	long long *sumxdcmtransp = calloc(reps, sizeof(long long));
	long long *sumdiftransp = calloc(reps, sizeof(long long));
	
	double tmult = 0.0, ttransp = 0.0;
// 	int retval, EventSet = PAPI_NULL;
// 	long_long values[1];
	long long dcmnorm = 0, dcmtransp = 0;
	long long xdcmnorm = 0, xdcmtransp = 0;

// 	retval = PAPI_library_init(PAPI_VER_CURRENT);
// 	if (retval != PAPI_VER_CURRENT){
// 		fprintf(stderr, "PAPI library init error!\n");
// 		exit(1);
// 		
// 	}
// 	if (PAPI_thread_init((long unsigned int (*)(void)) omp_get_thread_num) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
// 	if (PAPI_create_eventset(&EventSet) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
// 	if (PAPI_add_event(EventSet, PAPI_L1_DCM) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
// 	if (PAPI_start(EventSet) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	
	omp_sched_t kind;
	int chunk;
	omp_get_schedule(&kind, &chunk);
	
	char buf_t[256];
	
	mat_t *G, *Gtransp;
	
	G = malloc(sizeof(mat_t));
	Gtransp = malloc(sizeof(mat_t));
	G->size = A->size;
	Gtransp->size = A->size;
	
	double *xrow = calloc(A->size, sizeof(double));
	
	double *s = calloc(A->size, sizeof(double));
	double *r = calloc(A->size, sizeof(double));
	double *rold = calloc(A->size, sizeof(double));
	double *d = calloc(A->size, sizeof(double));
	double *q = calloc(A->size, sizeof(double));
	double *tmp = calloc(A->size, sizeof(double));
	double *dumm = calloc(A->size, sizeof(double));
	double *residuals = malloc(imax * sizeof(double));
	double *elapses = malloc(imax * sizeof(double));
	double *dnew = malloc(imax * sizeof(double));
	
	
// 	poweredA(A, A->size);
	
	
	
	start_usec = omp_get_wtime();											/* PRECONDITIONER + TRANSPOSED */
	
	
	
	
	
	precond(A, G, r);
// 	precond_vf(A, G, r);
	
	
	
	
	end_usec = omp_get_wtime();
	printf("\nPreconditioning Time: %lf\n", end_usec - start_usec);
	
	
	transposeCSC(G, Gtransp);
	start_usec = omp_get_wtime();
	printf("Transposing Time: %lf\n", start_usec - end_usec);
	
	
	char buf_p[256];
	snprintf(buf_p, sizeof buf_p, "../Outputs/cg/DataCG/patterns/pattern_%s_%i_%i.mtx", argv, percentpattern, patternpower);
	print_pattern(buf_p, G, G->size/10);
	
	
	mat_t *GtranspCOO = malloc(sizeof(mat_t));
	GtranspCOO->size = G->size;
	GtranspCOO->nnz = G->nnz;
	GtranspCOO->values = malloc(G->nnz*sizeof(double));
	GtranspCOO->cols = malloc(G->nnz*sizeof(unsigned int));
	GtranspCOO->rows = malloc(G->nnz*sizeof(unsigned int));
	
	
// 	printf("%u\n\n", G->nnz);
	
	unsigned int nthreads = omp_get_max_threads(); 
	unsigned int *limits = calloc(nthreads + 1, sizeof(unsigned int));
	
	CSCtoCOO(Gtransp, GtranspCOO, limits, nthreads);
	
	free(Gtransp);
	
	
	for (m = 0; m < reps; m++){														/* REPEAT LOOP reps TIMES */
		for (i = 0; i < A->size; i++){
			r[i] = 0.0;
			rold[i] = 0.0;
			d[i] = 0.0;
			q[i] = 0.0;
			s[i] = 0.0;
			tmp[i] = 0.0;
			x[i] = 0.0;
		}
		for (i = 0; i < imax; i++){
			elapses[i] = 0.0;
			residuals[i] = 0.0;
			dnew[i] = 0.0;
		}
		tmult = 0.0;
		ttransp = 0.0;
		
		
		pmultMatVect(tmp, r, G);														/* GENERATE INITIAL X */
		pmultMatVectCOO(x, tmp, GtranspCOO, limits, nthreads);
		
		pmultMatVect(tmp, x, A);														/* r = b - Ax */
		psubtVects(b, tmp, r, A->size);
		
		pmultMatVect(tmp, r, G);														/* d = M⁻1*r */
		pmultMatVectCOO(d, tmp, GtranspCOO, limits, nthreads);
		
		/* GET CACHE MISSES */
		/* DUMMY OPS TO GET xdcmnorm */
		
// 		zerovect(tmp, A->size);
// 		PAPI_reset(EventSet);
// 		pmultMatVect_DUMM(tmp, G);
// 		if (PAPI_read(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
// 		xdcmnorm = values[0];
// 		zerovect(dumm, A->size);
// 		PAPI_reset(EventSet);
// 		pmultMatVectCOO_DUMM(dumm, GtranspCOO, limits, nthreads);
// 		if (PAPI_read(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
// 		xdcmtransp = values[0];
		
		
		/* END DUMM */
		/* GET TOTAL DCM */
			

// 		pzerovect(tmp, A->size);
// 		PAPI_reset(EventSet);
// 		pmultMatVect(tmp, r, G);
// 		if (PAPI_read(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
// 		dcmnorm = values[0];
// 		pzerovect(s, A->size);
// 		PAPI_reset(EventSet);
// 		pmultMatVectCOO(s, tmp, GtranspCOO, limits, nthreads);
// 		if (PAPI_read(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
// 		dcmtransp = values[0];
		
		/* END GET CACHE MISSES */
		
		
		
		
		d_new = multVectVect(d, r, A->size);										/* d_new = d * r */
		double norm_b = sqrt(multVectVect(b, b, A->size));
		double norm_r = sqrt(multVectVect(r, r, A->size));
		start_usec = omp_get_wtime();
		for (i = 0; i < imax; i++){
// 			start_usec = omp_get_wtime();											/* BEGIN CG ITERATION */
			
			pmultMatVect(q, d, A);													/* q = A * d */
			alpha = d_new / pmultVectVect(d, q, A->size);							/* alpha = d_new / (d * q) */
			pscaleVect(alpha, d, tmp, A->size);										/* x += alpha * d */
			paddToVect(x, tmp, A->size);
			
			if (i%50 == 0){															/* Residual correction */
				pmultMatVect(tmp, x, A);											/* r = b - Ax */
				psubtVects(b, tmp, r, A->size);
			}
			else {
				pscaleVect(alpha, q, tmp, A->size);									/* r -= alpha * q */
				psubToVect(r, tmp, A->size);
			}
		
		
		
// // 			ttmp = omp_get_wtime();													/* Preconditioning residual */
// 			pmultMatVect(tmp, r, G);												/* s = M⁻1 * r    --->   s = G * G^t * r */
// // 			tmult += omp_get_wtime() - ttmp;
// 			pzerovect(s, A->size);
// // 			ttmp = omp_get_wtime();
// 			pmultMatVectCOO(s, tmp, GtranspCOO, limits, nthreads);
// // 			ttransp += omp_get_wtime() - ttmp;
			
			
			pequalvects(rold, r, A->size);
			pmultMatVect(tmp, r, G);												/* s = M⁻1 * r    --->   s = G * G^t * r */
			pzerovect(r, A->size);
			pmultMatVectCOO(r, tmp, GtranspCOO, limits, nthreads);
			pequalvects(s, r, A->size);
			pequalvects(r, rold, A->size);
			
			
			
			d_old = d_new;															/* d_old = d_new */
			d_new = pmultVectVect(s, r, A->size);									/* d_new = r * s */
			dnew[i] = d_new;
			pscaleVect(d_new/d_old, d, tmp, A->size);								/* d = s + (d_new/d_old) * d */
			paddVects(s, tmp, d, A->size);
																					/* END CG ITERATION */
// 			end_usec = omp_get_wtime();
// 			elapses[i] = end_usec - start_usec;
			norm_r = sqrt(cblas_ddot(A->size, r, 1, r, 1));
			residuals[i] = norm_r/norm_b;
			
			
			if ((i%10 == 0)){
				
				pmultMatVect(tmp, b, G);
				pzerovect(xrow, A->size);
				pmultMatVectCOO(xrow, tmp, GtranspCOO, limits, nthreads);
				
				char bufres[256];
				snprintf(bufres, sizeof bufres, "../Outputs/cg/DataCG/residual/resid_%s_%i_%i_%i", matname, patternpower, percentpattern, i);
				FILE *resfig = fopen(bufres, "a");
				double sum = 0.0, sumxrow = 0.0;
				int pos = 0;
				double avg = 0.0, avgxrows = 0.0;
				for(int j = 0; j < A->size; j++){
					sum += r[j]*r[j];
					sumxrow += (x[j] - xrow[j]) * (x[j] - xrow[j]);
					if ((j%1000 == 0) && (j > 0)){
						avg = sum / 1000.0;
						avgxrows = sumxrow / 1000.0;
						fprintf(resfig, "%i %.20lf %.20lf\n", pos, avg, avgxrows);
						pos++;
						sum = 0.0;
						sumxrow = 0.0;
					}
				}
				fclose(resfig);
			}
			
			if ( isless(residuals[i], err) ) {
// 				end_usec = 0.0;
// 				for (int j = 0; j < i; j++) end_usec += elapses[j];
				fprintf(stdout, "Precision reached; Iter: %i; Error: %.2e\n\n", i, residuals[i]);
				fprintf(stderr, "Precision reached; Iter: %i; Error: %.2e\n\n", i, residuals[i]);
				printf("RESULTING X\n\n");
				printarr(x, A->size);
				break;
			}
		}
		end_usec = omp_get_wtime();
		
// 		printf("tmult: %lf\t\tdcmmult: %.4e\txdcmmult: %.4e\nttransp: %lf\tdcmtransp: %.4e\txdcmtransp: %.4e\n\n", tmult/((double) i), (double) dcmnorm/((double) i), (double) xdcmnorm, ttransp/((double) i), (double) dcmtransp/((double) i), (double) xdcmtransp);
		
		if (i == (imax)) printf("No convergence. Residual = %.4e\n\n", residuals[i-1]);
		
		
		
		
		fprintf(outmn, "%.2lf\t", 100.0 * ((double) G->nnz) / (double) (((A->nnz - A->size)/2) + A->size) - 100.0);
		
		fprintf(outmn, "%u\t", percentpattern);
		fprintf(outmn, "%u\t", G->nnz);
		
		fprintf(outmn, "%u\t", i);
		fprintf(outmn, "%lf\t", end_usec - start_usec);
		fprintf(outmn, "%lf\t", (end_usec - start_usec)/((double) i));
		
		fprintf(outmn, "%lf\t", tmult/((double) i));
		fprintf(outmn, "%lf\t", ttransp/((double) i));
		
		fprintf(outmn, "%lli\t", dcmnorm);
		fprintf(outmn, "%lli\t", xdcmnorm);
		fprintf(outmn, "%lli\t", dcmnorm - xdcmnorm);
		
		fprintf(outmn, "%lli\t", dcmtransp);
		fprintf(outmn, "%lli\t", xdcmtransp);
		fprintf(outmn, "%lli\n", dcmtransp - xdcmtransp);
		
		
		sumiter[m] = i;
		sumtottime[m] = end_usec - start_usec;
		sumtimeiter[m] = (end_usec - start_usec)/((double) i);
		
		sumtmult[m] = tmult/((double) i);
		sumttransp[m] = ttransp/((double) i);
		
		sumdcmnorm[m] = dcmnorm;
		sumxdcmnorm[m] = xdcmnorm;
		sumdifmult[m] = dcmnorm - xdcmnorm;
		sumdcmtransp[m] = dcmtransp;
		sumxdcmtransp[m] = xdcmtransp;
		sumdiftransp[m] = dcmtransp - xdcmtransp;
		
		
		
		if (m >= 0) {
			snprintf(buf_t, sizeof buf_t, "../Outputs/cg/DataCG/logs/pcg_%s_%i_%i_%i_%i.log", argv, percentpattern, patternpower, m, rhs);
			dump_info(buf_t, i, residuals, elapses, dnew);
		}
	}
	
	
	// IF GETTING MEDIAN
	
// 	double *sumtottimeOrdered = calloc(reps, sizeof(double));
// 	
// 	m = floor(((double) reps)/2.0);
// 	
// 	for (i = 0; i < reps; i++) sumtottimeOrdered[i] = sumtottime[i];
// 	
// 	qsort(sumtottimeOrdered, reps, sizeof(double), compd);
// 	
// 	for (i = 0; i < reps; i++){
// // 		printf("%lf, ", sumtottimeOrdered[i]);
// // 		printf("%lf, ", sumtottime[i]);
// // 		printf("\n");
// 		if (sumtottime[i] == sumtottimeOrdered[m]) break;
// 	}
	
	// IF GETTING LOWEST
	
	double tmin = 1000.0;
	for (unsigned int j = 0; j < reps; j++){
		if (sumtottime[j] < tmin) {
			i = j;
			tmin = sumtottime[j];
		}
	}
	
	
	char buffig[256];
	snprintf(buffig, sizeof buffig, "../Outputs/cg/DataCG/Fig_DATA_%s_%i_%i", matname, patternpower, rhs);
	
	
	FILE *outmnfig = fopen(buffig, "a");
	if (outmn == NULL){
		printf("Could not open writing file.");
		return 0;
	}
	
	
	
	
// 	fprintf(outmn, "%i\t", omp_get_max_threads());

	fprintf(outmnfig, "%.2lf\t", 100.0 * ((double) G->nnz) / (double) (((A->nnz - A->size)/2) + A->size) - 100.0);
	fprintf(outmnfig, "%u\t", percentpattern);
	fprintf(outmnfig, "%u\t", G->nnz);
// 	fprintf(outmnfig, "%u\t", (int) floor(((double) sumiter[m]) / ((double) meannums)));
	fprintf(outmnfig, "%u\t", (int) floor(((double) sumiter[i]) / ((double) meannums)));
	
	
	fprintf(outmnfig, "%lf\t", sumtottime[i] / ((double) meannums));
	fprintf(outmnfig, "%lf\t", sumtimeiter[i] / ((double) meannums));
	fprintf(outmnfig, "%lf\t", sumtmult[i] / ((double) meannums));
	fprintf(outmnfig, "%lf\t", sumttransp[i] / ((double) meannums));
	
	fprintf(outmnfig, "%lli\t", (long long) floor(((double) sumdcmnorm[i]) / ((double) meannums)));
	fprintf(outmnfig, "%lli\t", (long long) floor(((double) sumxdcmnorm[i]) / ((double) meannums)));
	fprintf(outmnfig, "%lli\t", (long long) floor(((double) sumdifmult[i]) / ((double) meannums)));
	
	fprintf(outmnfig, "%lli\t", (long long) floor(((double) sumdcmtransp[i]) / ((double) meannums)));
	fprintf(outmnfig, "%lli\t", (long long) floor(((double) sumxdcmtransp[i]) / ((double) meannums)));
	fprintf(outmnfig, "%lli\t", (long long) floor(((double) sumdiftransp[i]) / ((double) meannums)));
	fprintf(outmnfig, "%u\n", i);
	
	fclose(outmnfig);
	
// 	if (PAPI_stop(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	
	
	
	
	
/*	
	
	
	unsigned int dim = A->size;
	int testcol = 0;
	int cachesize = 8;
	
	unsigned int *positions = calloc(dim, sizeof(unsigned int));
	
	int cachecol = 0;
	long long cacheline = 0;
	
	for (i = A->rows[dim - 1]; i < A->rows[dim]; i++){
		cachecol = (int) A->cols[i];
		cacheline = (((long long) &r[A->cols[i]]%64)/8);
		
		for (int k = 0; k < cachesize; k++){
			testcol = cachecol - (int) cacheline + k;
			if ((testcol >= 0) && (testcol <= i)){
				positions[testcol] = 1;
			}
			
		}
		
	}
	
	
	for (i = 0; i < dim; i++){
		if (positions[i] == 1){
			printf("%i\t", i);
		}
	}
	printf("\n");
	
	for (i = 0; i < dim; i++){
		if (positions[i] == 1){
// 			printf("%.2lf\t", x[i]);
			printf("%.2lf\t", getvalue_mat(dim - 1, i, G));
			
		}
	}
	printf("\n");
	
	
	
	
	
	*/
	
	
	
	
	char buf_x[256];
	snprintf(buf_x, sizeof buf_x, "../Outputs/cg/DataCG/x_%s.log", matname);
	FILE *logt = fopen(buf_x, "w");
	for ( int i = 0; i < A->size; i++ ) {
		fprintf(logt, "%i\t%.8e\t%.8e\n", i, x[i], initx[i]);
	}
	fclose(logt);
	
	
	
	
	
	
	
	
	free(dumm); free(r); free(d); free(q); free(s); free(tmp); free(G); free(residuals); free(elapses); free(dnew);
	
	free(sumiter); free(sumtottime); free(sumtimeiter); free(sumtmult); free(sumttransp); free(sumdcmnorm); free(sumxdcmnorm); free(sumdifmult); free(sumdcmtransp); free(sumxdcmtransp); free(sumdiftransp); 
// 	free(sumtottimeOrdered);
	
	return x;
}





void lu(mat_t *A, double *x, double *b, double *initx, double *initb){
	
	unsigned int dim = A->size;
	
	
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




















