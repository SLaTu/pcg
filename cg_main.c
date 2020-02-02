
#include "cg_main.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define PI 3.14159265
#define meannums 10

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
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	double *b = calloc(A->size, sizeof(double));
	double *x = calloc(A->size, sizeof(double));
	
	
	printf("NNZ:\t%u\n", A->nnz);
	printf("Size:\t%u\n\n", A->size);
	
	char buf[256];
	snprintf(buf, sizeof buf, "../Outputs/cg/DataCG/DATA_%s.txt", argv[1]);
	
	outmn = fopen(buf, "a");
	if (outmn == NULL){
		printf("Could not open writing file.");
		return 0;
	}
	
	for (j = 0; j < A->size; j++) x[j] = sin(((double) j) * PI/18000);
	pmultMatVect(b, x, A);
	for (j = 0; j < A->size; j++) x[j] = 0.0;

	switch ( mode ) {
		case 1: 
			printf("Mode:\tBase CG\n");
			x = cg_base(A, x, b);
			break;
			printf("Mode:\tD⁻1 Preconditioned CG\n");
			x = cg_precond_diag(A, x, b);
			break;
		case 3: 
			printf("Mode:\tPCG\n");
			x = cg_precond(A, x, b, argv[1]);
			break;
		default:
			printf("No algorithm selected.\n");
			break;
	}

	
	fclose (outmn);
	free(A);
	free(b);
	free(x);
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

double *cg_precond_diag(mat_t *A, double *x, double *b){

	int i = 0;
	double *d;
	double *r;
	double d_new = 0.0; 
	double d_0 = 0.0;
	double *q;
	double *s;
	double alpha = 0.0;
	double d_old = 0.0;
	double beta = 0.0;
	double start_usec = 0.0, end_usec = 0.0;
	double *tmp;
	
	mat_t *D_1;
	
	D_1 = malloc(sizeof(mat_t));
	
	D_1->values = malloc(A->size*sizeof(double));
	D_1->cols = malloc(A->size*sizeof(unsigned int));
	D_1->rows = malloc((A->size + 1)*sizeof(unsigned int));
	
	r = malloc(A->size*sizeof(double));
	d = malloc(A->size*sizeof(double));
	q = malloc(A->size*sizeof(double));
	s = malloc(A->size*sizeof(double));
	tmp = malloc(A->size*sizeof(double));
	
	// Inverse diagonal A
	D_1 = inverseD(A);
	
	// r = b - A * x
	multMatVect(tmp, x, A);
	subtVects(b, tmp, r, A->size);
	
	// d = M^-1 * r      // MODIFY
	
	multMatVect(d, r, D_1); // instead of A D^-1
	
	// dnew = r^T * d
	d_new = multVectVect(r, d, A->size);
	d_0 = d_new;

	start_usec = omp_get_wtime();
	while ((i < imax)&&(d_new > err * err * d_0)){
		
		// q = A * d
		multMatVect(q, d, A);
		
		// alpha = dnew / (d^T * q)
		alpha = d_new / multVectVect(d, q, A->size);
		
		// x = x + alpha * d
		scaleVect(alpha, d, tmp, A->size);
		addToVect(x, tmp, A->size);
		
		// 
		if (i%50 == 0){
			// r = b - A * x
			multMatVect(tmp, x, A);
			subtVects(b, tmp, r, A->size);
		}
		else {
			// r = r - alpha * q
			scaleVect(alpha, q, tmp, A->size);
			subToVect(r, tmp, A->size);
		}
		// PRECONDITIONING STEP -- s = M-1 * r -- s = G^t * G * r
		// G = lower triangular part of matrix A^L 
		
		multMatVect(s, r, D_1); // instead of A M^-1
		// dold = dnew
		d_old = d_new;
		
		// dnew = r^T * s
		d_new = multVectVect(r, s, A->size);
		
		// beta = dnew / dold
		beta = d_new/d_old;
		
		// d = s + beta * d
		scaleVect(beta, d, tmp, A->size);
		addVects(s, tmp, d, A->size);
		
		//
		if (i%200==0){
			reperror(A, r, b, i);
		}
		
		i++;
	}
	end_usec = omp_get_wtime();
	
	report(A, r, b, i, start_usec, end_usec);
	
        fprintf(outmn, "%u\t", i);
        fprintf(outmn, "%lf\t", end_usec - start_usec);
        fprintf(outmn, "%lf\n", (end_usec - start_usec)/((double) i));

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

double *cg_precond(mat_t *A, double *x, double *b, char *argv){

	unsigned int i = 0;
	double d_new = 0.0; 
	double alpha = 0.0;
	double d_old = 0.0;
	double start_usec = 0.0, end_usec = 0.0;
	
	
	int sumiter = 0;
	double sumtottime = 0.0;
	double sumtimeiter = 0.0;
	double sumtmult = 0.0;
	double sumttransp = 0.0;
	long long sumdcmnorm = 0;
	long long sumxdcmnorm = 0;
	long long sumdifmult = 0;
	long long sumdcmtransp = 0;
	long long sumxdcmtransp = 0;
	long long sumdiftransp = 0;
	
	double tmult = 0.0, ttransp = 0.0, ttmp = 0.0;
	int retval, EventSet = PAPI_NULL;
	long_long values[1];
	long long dcmnorm = 0, dcmtransp = 0;
	long long xdcmnorm = 0, xdcmtransp = 0;

	retval = PAPI_library_init(PAPI_VER_CURRENT);
	
	if (retval != PAPI_VER_CURRENT){
		fprintf(stderr, "PAPI library init error!\n");
		exit(1);
		
	}

	if (PAPI_thread_init((long unsigned int (*)(void)) omp_get_thread_num) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	
	if (PAPI_create_eventset(&EventSet) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	if (PAPI_add_event(EventSet, PAPI_L1_DCM) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	
	if (PAPI_start(EventSet) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	
	mat_t *G, *Gtransp;
	
	G = malloc(sizeof(mat_t));
	Gtransp = malloc(sizeof(mat_t));
	G->size = A->size;
	Gtransp->size = A->size;
	

	double *s = (double*) aligned_alloc(64, A->size*sizeof(double));
	double *r = (double*) aligned_alloc(64, A->size*sizeof(double));
	double *d = calloc(A->size, sizeof(double));
	double *q = calloc(A->size, sizeof(double));
	double *tmp = calloc(A->size, sizeof(double));
	double *dumm = calloc(A->size, sizeof(double));
	double *residuals = malloc(imax * sizeof(double));
	double *elapses = malloc(imax * sizeof(double));
	double *dnew = malloc(imax * sizeof(double));
	
	start_usec = omp_get_wtime();														/* PRECONDITIONER + TRANSPOSED */
	precond(A, G, r);
	end_usec = omp_get_wtime();
	printf("Preconditioning Time: %lf", end_usec - start_usec);
	transposeCSC(G, Gtransp);
	start_usec = omp_get_wtime();
	printf("\nTransposing Time: %lf\n\n", start_usec - end_usec);
	
	
	for (int m = 0; m < reps; m++){														/* REPEAT LOOP reps TIMES */
		for (i = 0; i < A->size; i++){
			r[i] = 0.0;
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
		pmultMatVectCSC(x, tmp, Gtransp);
		
		pmultMatVect(tmp, x, A);														/* r = b - Ax */
		psubtVects(b, tmp, r, A->size);
		
		pmultMatVect(tmp, r, G);														/* d = M⁻1*r */
		pmultMatVectCSC(d, tmp, Gtransp);
		
		/* GET CACHE MISSES */
		/* DUMMY OPS TO GET xdcmnorm */
		
		zerovect(tmp, A->size);
		PAPI_reset(EventSet);
		pmultMatVect_DUMM(tmp, G);
		if (PAPI_read(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
		xdcmnorm = values[0];
		zerovect(dumm, A->size);
		PAPI_reset(EventSet);
		pmultMatVectCSC_DUMM(dumm, Gtransp);
		if (PAPI_read(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
		xdcmtransp = values[0];
		
		
		/* END DUMM */
		/* GET TOTAL DCM */
			

		pzerovect(tmp, A->size);
		PAPI_reset(EventSet);
		pmultMatVect(tmp, r, G);
		if (PAPI_read(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
		dcmnorm = values[0];
		pzerovect(s, A->size);
		PAPI_reset(EventSet);
		pmultMatVectCSC(s, tmp, Gtransp);
		if (PAPI_read(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
		dcmtransp = values[0];
		
		/* END GET CACHE MISSES */
		
		
		d_new = multVectVect(d, r, A->size);										/* d_new = d * r */
		double norm_b = multVectVect(b, b, A->size);
		double norm_r = sqrt(multVectVect(r, r, A->size));
		for (i = 0; i < imax; i++){
			start_usec = omp_get_wtime();											/* BEGIN CG ITERATION */
			
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
			
			ttmp = omp_get_wtime();													/* Preconditioning residual */
			pmultMatVect(tmp, r, G);												/* s = M⁻1 * r    --->   s = G * G^t * r */
			tmult += omp_get_wtime() - ttmp;
			pzerovect(s, A->size);
			ttmp = omp_get_wtime();
			pmultMatVectCSC(s, tmp, Gtransp);
			ttransp += omp_get_wtime() - ttmp;
			
			d_old = d_new;															/* d_old = d_new */
			d_new = pmultVectVect(s, r, A->size);									/* d_new = r * s */
			dnew[i] = d_new;
			pscaleVect(d_new/d_old, d, tmp, A->size);								/* d = s + (d_new/d_old) * d */
			paddVects(s, tmp, d, A->size);
																					/* END CG ITERATION */
			end_usec = omp_get_wtime();
			elapses[i] = end_usec - start_usec;
			norm_r = sqrt(cblas_ddot(A->size, r, 1, r, 1));
			residuals[i] = norm_r/norm_b;
			
			if ( isless(residuals[i], err) ) {
				end_usec = 0.0;
				for (int j = 0; j < i; j++) end_usec += elapses[j];
				fprintf(stdout, "Precision reached; Iter: %i; Error: %.2e; Total time: %lf\n\n", i, residuals[i], end_usec);
				break;
			}
		}
		
		
		printf("tmult: %lf\t\tdcmmult: %.4e\txdcmmult: %.4e\nttransp: %lf\tdcmtransp: %.4e\txdcmtransp: %.4e\n\n", tmult/((double) i), (double) dcmnorm/((double) i), (double) xdcmnorm, ttransp/((double) i), (double) dcmtransp/((double) i), (double) xdcmtransp);
		
		if (i == (imax)) printf("No convergence. Residual = %.4e\n\n", residuals[i]);
		
		
		
		
		fprintf(outmn, "%i\t", omp_get_max_threads());
		
		fprintf(outmn, "%u\t", percentpattern);
		fprintf(outmn, "%u\t", G->nnz);
		
		fprintf(outmn, "%u\t", i);
		fprintf(outmn, "%lf\t", end_usec);
		fprintf(outmn, "%lf\t", (end_usec)/((double) i));
		
		fprintf(outmn, "%lf\t", tmult/((double) i));
		fprintf(outmn, "%lf\t", ttransp/((double) i));
		
		fprintf(outmn, "%lli\t", dcmnorm);
		fprintf(outmn, "%lli\t", xdcmnorm);
		fprintf(outmn, "%lli\t", dcmnorm - xdcmnorm);
		
		fprintf(outmn, "%lli\t", dcmtransp);
		fprintf(outmn, "%lli\t", xdcmtransp);
		fprintf(outmn, "%lli\n", dcmtransp - xdcmtransp);
		
		
		
		
		if((m >= (reps - meannums)) & (m < reps)){
			sumiter += i;
			sumtottime += end_usec;
			sumtimeiter += (end_usec)/((double) i);
			sumtmult += tmult/((double) i);
			sumttransp += ttransp/((double) i);
			
			sumdcmnorm += dcmnorm;
			sumxdcmnorm += xdcmnorm;
			sumdifmult += dcmnorm - xdcmnorm;
			sumdcmtransp += dcmtransp;
			sumxdcmtransp += xdcmtransp;
			sumdiftransp += dcmtransp - xdcmtransp;
			
		}
	}
	
	
	
	
	fprintf(outmn, "\n%i\t", omp_get_max_threads());
	fprintf(outmn, "%u\t", percentpattern);
	fprintf(outmn, "%u\t", G->nnz);
	fprintf(outmn, "%u\t", (int) floor(((double) sumiter) / ((double) meannums)));
	
	
	
	fprintf(outmn, "%lf\t", sumtottime / ((double) meannums));
	fprintf(outmn, "%lf\t", sumtimeiter / ((double) meannums));
	fprintf(outmn, "%lf\t", sumtmult / ((double) meannums));
	fprintf(outmn, "%lf\t", sumttransp / ((double) meannums));
	
	fprintf(outmn, "%lli\t", (long long) floor(((double) sumdcmnorm) / ((double) meannums)));
	fprintf(outmn, "%lli\t", (long long) floor(((double) sumxdcmnorm) / ((double) meannums)));
	fprintf(outmn, "%lli\t", (long long) floor(((double) sumdifmult) / ((double) meannums)));
	
	fprintf(outmn, "%lli\t", (long long) floor(((double) sumdcmtransp) / ((double) meannums)));
	fprintf(outmn, "%lli\t", (long long) floor(((double) sumxdcmtransp) / ((double) meannums)));
	fprintf(outmn, "%lli\n\n", (long long) floor(((double) sumdiftransp) / ((double) meannums)));
	
	
	
	
	omp_sched_t kind;
	int chunk;
	omp_get_schedule(&kind, &chunk);
	
	char buf_t[256];
	snprintf(buf_t, sizeof buf_t, "../Outputs/cg/DataCG/pcg_%s_%d_%d.log", argv, kind, chunk);
	dump_info(buf_t, i, residuals, elapses, dnew);
	
	
	if (PAPI_stop(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	
	free(dumm); free(r); free(d); free(q); free(s); free(tmp); free(G); free(Gtransp); free(residuals); free(elapses); free(dnew);

	
	return x;
}


