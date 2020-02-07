
#include "cg_aux.h"

extern char matname[256];
extern int rhs;
extern int mode;
extern int imax;
extern double err;
extern unsigned int bsize; 
extern FILE *fp;
extern unsigned int patternmode;
extern unsigned int percentpattern;
extern unsigned int patternpower;
extern unsigned int reps;

int cg_config(int argc, char *argv[]) {

	if (argc < 7){
		printf("\nAdd as arguments a matrix file, block size, maximum iterations, tolerance (<1) and mode. For mode 3 (PCG) add PCG pattern. For PCG pattern 3 add power and percentage. For PCG pattern 4 add percentage.\n\n");
		return 1;
	}

	snprintf(matname, sizeof matname, "%s", argv[1]);
	rhs = atoi(argv[2]);
	reps = atoi(argv[3]);
	imax = atoi(argv[4]);
	err = atof(argv[5]);
	mode = atoi(argv[6]);
	

	if ((mode == 3) && (argc >= 8)) {
		patternmode = atoi(argv[6]);
		
		if ((patternmode == 4) && (argc >= 9)) {
			percentpattern = atoi(argv[8]);
		}
		else if ((patternmode == 4) && (argc < 9)) {
			fprintf(stderr, "\nError: Add percentage of elements to add with PCG 4 (e.g.: 100).\n\n");
			return 1;
		}
		
		
		if ((patternmode == 3) && (argc >= 10)) {
			percentpattern = atoi(argv[8]);
			patternpower = atoi(argv[9]);
		}
		else if ((patternmode == 3) && (argc < 10)) {
			fprintf(stderr, "\nError: Add percentage of elements to add with PCG 3 (e.g.: 100) and pattern power (e.g.: 2).\n\n");
			return 1;
		}
	}
	else if ((mode == 3) && (argc < 8)) {
		patternmode = 0;
	}
	
	if (mode != 3) {patternmode = 0; percentpattern = 0; patternpower = 0;}
	
	printf("\nMatrix:\t%s\n", matname);
	printf("Imax:\t%d\n", imax);
	printf("Error:\t%.2e\n", err);
	
	return 0;
}


int cg_setup(char* matname, mat_t *mat){

	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	unsigned int j = 0, i = 0;
	unsigned int rowvalue = 0, cmp = 0, acum = 0, position = 0, count = 0;
	double num = 0;
	
	char buff[256];
	
	snprintf(buff, sizeof buff, "../Inputs/Matrices/%s", matname);
	fp = fopen (buff, "r");
	if (!fp){
		fprintf(stderr, "Error: file open failed '%s'.\n\n", matname);
		return 1;
	}

	while ((read = getline(&line, &len, fp)) != -1){

		num = 0;

		while(fscanf(fp, "%lf", &num) == 1){

			if (count == 0){
				mat->size = (unsigned int) num;
			}
			else if (count == 1){}
			else if (count == 2){
				mat->nnz = (unsigned int) num;
				mat->values = malloc(mat->nnz*sizeof(double));
				mat->cols = malloc(mat->nnz*sizeof(unsigned int));
				mat->rows = malloc((mat->size + 1) * sizeof(unsigned int));
				mat->rows[0] = 0;
			}
			else if (count == 4){
				mat->cols[((i - 1)/3) - 1] = (unsigned int) num - 1;
			}
			else if (count == 3){
				rowvalue = num - 1;
				if ((rowvalue-cmp) > 1){
					for (j = cmp; j < rowvalue; j++){
						position++;
						mat->rows[position] = acum;
					}
					cmp = rowvalue;
					acum++;
				}
				else if ((rowvalue - cmp) == 1){
					position++;
					cmp = rowvalue;
					mat->rows[position] = acum;
					acum++;
				}
				else{acum++;}
			}
			else if (count == 5){
				mat->values[(i - 2)/3 - 1] = num;
			}
			else {}
			i++;
			count++;
			if (count == 6) count = 3;
		}
	}

	for (i = (position + 1); i < (mat->size + 1); i++){mat->rows[i] = acum;}

	if (fp != stdin) fclose (fp);   /* close file if not stdin */
	if (line) free(line);
	
	
	
	
	for (i = 0; i < mat->nnz; i++) printf("%lf, ", mat->values[i]);
	printf("\n");
	for (i = 0; i < mat->nnz; i++) printf("%i, ", mat->cols[i]);
	printf("\n");
	for (i = 0; i < mat->size + 1; i++) printf("%i, ", mat->rows[i]);
	printf("\n");
	
	

	return 0;
}


int readrhs(double *b, char* matname, unsigned int size){

	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	unsigned int i = 0;
	unsigned int index = 0;
	double num = 0;
	
	char buff[256];
	unsigned int count = 0;
	
	snprintf(buff, sizeof buff, "../Inputs/RHS/RHS_%s", matname);
	fp = fopen (buff, "r");
	if (!fp){
		fprintf(stderr, "Error: file open failed 'RHS_%s'.\n\n", matname);
		return 1;
	}
	
	while ((read = getline(&line, &len, fp)) != -1){
		num = 0;
		while(fscanf(fp, "%lf", &num) == 1){
			if (count == 0){
				if ((unsigned int) num == size) {printf("\nOK\n");}
				else {return 1;}
			}
			else if (count == 1){}
			else if (count == 2){
				index = (unsigned int) num;
			}
			else if (count == 3){
				b[index - 1] = (double) num;
			}
			else {}
			count++;
			if (count == 4) count = 2;
		}
	}

	if (fp != stdin) fclose (fp);   /* close file if not stdin */
	if (line) free(line);
	
	for (i = 0; i < size; i++) printf("%lf, ", b[i]);
	printf("\n");

	return 0;
}


void report(mat_t *A, double *r, double *b, unsigned int i, double start_usec, double end_usec){

	unsigned int j = 0;
	double norm_b = 0.0;
	double norm_r = 0.0;

	if (i == 0) {printf("Something bad happened. Probably related to patterns.\n");}
	else {
		for (j = 0; j < A->size; j++) norm_r += (r[j]*r[j]);
		norm_r = sqrt(norm_r);
		
		for (j = 0; j < A->size; j++) norm_b += (b[j]*b[j]);
		
		if (i == imax) {printf("\nSolution wouldn't converge.\n"); printf("Error:\t%.2e\n", norm_r / norm_b);}
		else {printf("\nSolution converged in %u iteration/s.\n" , i); printf("Error:\t%.2e\n", norm_r / norm_b);}
	}
	
	printf("Total Time:\t%.4lf s\tIteration Time:\t%.4lf s\n\n", end_usec - start_usec, (end_usec - start_usec) / ((double) i));

}



void reperror(mat_t *A, double *r, double *b, unsigned int i){
	
	double berror = 0.0;
	double error = 0.0;
	unsigned int j = 0;
	for (j = 0; j < A->size; j++) error += (r[j]*r[j]);
	error = sqrt(error);
	for (j = 0; j < A->size; j++) berror += (b[j]*b[j]);
	berror = sqrt(berror);
	printf("Iter: %d;\tError: %.2e\n", i, error / berror);
	
}



mat_t *inverseD(mat_t *mat){

	unsigned int i = 0, j = 0;
	mat_t *D;

	D = malloc(sizeof(mat_t));
	D->size = mat->size;
	D->nnz = D->size;

	D->values = malloc(D->nnz*sizeof(double));
	D->cols = malloc(D->nnz*sizeof(unsigned int));
	D->rows = malloc((D->size + 1)*sizeof(unsigned int));

	for (i = 0; i < D->size; i++){
		for (j = mat->rows[i]; j < mat->rows[i + 1]; j++)
		{
			if (mat->cols[j] == i){
				D->values[i] = 1.0/mat->values[j];
			}
		}
		D->cols[i] = i;
		D->rows[i + 1] = i + 1;
	}
	return D;
}



double *cg_construct(mat_t *A, double *x, double *b, int maxiter, double error){

	int i = 0, j = 0;
	double *d;
	double d_new = 0.0;
	double *q;
	double *r;
	double alpha = 0.0;
	double d_old = 0.0;
	double beta = 0.0;
	double *tmp;

	r = malloc(A->size*sizeof(double));
	d = malloc(A->size*sizeof(double));
	q = malloc(A->size*sizeof(double));
	tmp = malloc(A->size*sizeof(double));

	multMatVect(tmp, x, A);
	subtVects(b, tmp, r, A->size);
	for (j = 0; j < A->size; j++) d[j] = r[j];
	d_new = multVectVect(r, r, A->size);
	double norm_b = multVectVect(b, b, A->size);
	double norm_r = sqrt(multVectVect(r, r, A->size));
	while ((i < maxiter)&&(norm_r/norm_b > error)){
		multMatVect(q, d, A);
		alpha = d_new / multVectVect(d, q, A->size);
		scaleVect(alpha, d, tmp, A->size);
		addToVect(x, tmp, A->size);
		scaleVect(alpha, q, tmp, A->size);
		subToVect(r, tmp, A->size);
		d_old = d_new;
		d_new = multVectVect(r, r, A->size);
		beta = d_new/d_old;
		scaleVect(beta, d, tmp, A->size);
		addVects(r, tmp, d, A->size);
		norm_r = sqrt(multVectVect(r, r, A->size));
		i++;
	}
	free(r); free(d); free(q); free(tmp);
	return x;
}












int LUPDecompose(double **A, int N, double Tol, int *P) {

    int i, j, k, imax; 
    double maxA, *ptr, absA;

    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA) { 
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) return 0; //failure, matrix is degenerate

        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            //counting pivots starting from N (for determinant)
            P[N]++;
        }

        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }

    return 1;  //decomposition done 
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
void LUPSolve(double **A, int *P, double *b, int N, double *x) {

	int i, k;
    for (i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }

    for (i = N - 1; i >= 0; i--) {
        for (k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];

        x[i] = x[i] / A[i][i];
    }
}






void precond(mat_t *mat, mat_t *G, double *xfinal){

	int i = 0, j = 0, k = 0;

	pat_t *pattern = malloc(sizeof(pat_t));
	pattern->rows = calloc((G->size + 1), sizeof(unsigned int));

	pat_t *expanded_patt = malloc(sizeof(pat_t));
	expanded_patt->rows = calloc((G->size + 1), sizeof(unsigned int));


	// ---------------------------------------------------------------

	// PATTERNS


	switch ( patternmode ) {
		case 1:
			printf("PCG:\tLOWER TRIANGLE PATTERN\n\n");
			ltp(mat, pattern, expanded_patt);
			break;
		case 2:
			printf("PCG:\tFULL LOWER TRIANGLE\n\n");
			flt(mat, pattern, expanded_patt);
			break;
		case 3:
			printf("PCG:\tPOWER %u OF A\n\n", patternpower);
			powerA(mat, pattern, expanded_patt, xfinal);
			break;
		case 4:
			printf("PCG:\tADD PERCENTAGE -> %u\n\n", percentpattern);
			perclt(mat, pattern, expanded_patt, xfinal);
			break;
		default:
			printf("No pattern selected. Running default LOWER TRIANGLE PATTERN\n\n");
			ltp(mat, pattern, expanded_patt);
			break;
	}
	
	G->nnz = expanded_patt->nnz;
	G->values = calloc(G->nnz, sizeof(double));
	G->cols = malloc(G->nnz * sizeof(unsigned int));
	G->rows = malloc((G->size + 1) * sizeof(unsigned int));
	
	//	3.
	//	Compute nonzero entries in G
	// Init G
	for (i = 0; i < G->nnz; i++){
		G->cols[i] = expanded_patt->cols[i];
	}
	for (i = 0; i < G->size + 1; i++){
		G->rows[i] = expanded_patt->rows[i];
	}
	
	double totalresid = 0.0;							/* to reduce when parallelizing */
	
	
// 	#pragma omp parallel for private(j, k) reduction(+:totalresid)
	for (i = 0; i < mat->size; i++){					// For every row in the matrix
		
		unsigned int rowelems = expanded_patt->rows[i + 1] - expanded_patt->rows[i];
		unsigned int *arrayindex = calloc(rowelems, sizeof(unsigned int));
		int n = rowelems, nrhs = 1, lda = rowelems, ldb = rowelems, info;
		int ipiv[rowelems];
		double *a = malloc(rowelems*rowelems*sizeof(double));
		double *a_res = malloc(rowelems*rowelems*sizeof(double));
		double *b = calloc(rowelems, sizeof(double));
		double *b_res = calloc(rowelems, sizeof(double));
		double *r = malloc(rowelems*sizeof(double));
		double *tmp = calloc(rowelems, sizeof(double));
		double norm_r = 0.0;
		
		for (j = expanded_patt->rows[i]; j < expanded_patt->rows[i + 1]; j++){
			arrayindex[j - expanded_patt->rows[i]] = expanded_patt->cols[j];
		}
		
		unsigned int poscount = 0;
		
		for (j = 0; j < rowelems; j++){					// Rows construct_mat
			for (k = 0; k < rowelems; k++){				// Cols construct_mat
				a[poscount] = getvalue_mat(arrayindex[j], arrayindex[k], mat);
				a_res[poscount] = a[poscount];
				poscount++;
			}
			b[j] = 0.0;
			b_res[j] = 0.0;
// 			b[j] = a[poscount - 1];
// 			b_res[j] = a[poscount - 1];
		}
		b[rowelems - 1] = 1.0;
		b_res[rowelems - 1] = 1.0;
		
		
		/****************************************************************/
// 		for (j = 0; j < rowelems; j++) printf("%.16lf, ", b[j]);
// 		printf("\n");
		
		dgesv( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
		
// 		for (j = 0; j < rowelems; j++) printf("%.16lf, ", b[j]);
// 		printf("\n\n\n");
		
		double norm_x = cblas_ddot(rowelems, b, 1, b, 1);
		
		for (j = 0; j < rowelems; j++){
			for (k = 0; k < rowelems; k++){
				tmp[j] += a_res[j*rowelems + k] * b[k]; 
			}
			r[j] = b_res[j] - tmp[j];
			norm_r += r[j] * r[j];
		}
		
		totalresid += sqrt(norm_r)/sqrt(norm_x);
		
		
		for (j = 0; j < rowelems; j++){
			G->values[expanded_patt->rows[i] + j] = b[j];
		}
		
		free(b_res);
		free(tmp);
		free(a);
		free(a_res);
		free(b);
		free(r);
		free(arrayindex);
		
		
		/****************************************************************/	
/*		
		int counter = 0;
		double **Ap = malloc(rowelems*sizeof(double*));
		double tol = 1E-16;
		int *P = calloc(rowelems + 1, sizeof(int));
		double *x = calloc(rowelems, sizeof(double));
		
		
		for (j = 0; j < rowelems; j++){
			*(Ap + j) = malloc(rowelems * sizeof(double));
			for (k = 0; k < rowelems; k++){
				*(*(Ap + j) + k) = a[counter];
				counter++;
			}
		}
		
		LUPDecompose(Ap, rowelems, tol, P);
		LUPSolve(Ap, P, b, rowelems, x);
		
		double norm_x = multVectVect(x, x, rowelems);
		
		for (j = 0; j < rowelems; j++){
			for (k = 0; k < rowelems; k++){
				tmp[j] += a_res[j*rowelems + k] * x[k]; 
			}
			r[j] = b_res[j] - tmp[j];
			norm_r += r[j] * r[j];
		}
		
		totalresid += sqrt(norm_r)/sqrt(norm_x);
		
		for (j = 0; j < rowelems; j++){
			G->values[expanded_patt->rows[i] + j] = x[j];
		}
		
		free(Ap);
		free(P);
		free(x);
		free(b_res);
		free(tmp);
		free(a);
		free(a_res);
		free(b);
		free(r);
		free(arrayindex);
*/
		
	}
	
	
	printf("INVERSE: Accumulated Residual: %.2e\n", totalresid);
		
		
	double *diag;
	unsigned int diagcounter = 0, zerocounter = 0;
	diag = calloc(G->size, sizeof(double));
	double zero = 0.0000000000000001;
	unsigned int counter = 0;
	

	for (i = 0; i < G->size; i++){
		for (j = G->rows[i]; j < G->rows[i + 1]; j++){
			if (G->cols[j] == i){diag[i] = sqrt(fabs(G->values[j]));}
		}
		if (fabs(diag[i]) < zero) {
			diagcounter++;
		}
	}
	
	for (i = 0; i < G->size; i++){
		for (j = G->rows[i]; j < G->rows[i + 1]; j++){
			G->values[j] = G->values[j]/diag[i];
			if (fabs(G->values[j]) < zero) {
				zerocounter++;
			}
		}
	}
	
	printf("INVERSE: %u zeroes in diag\nINVERSE: %u zeroes in matrix\n\n", diagcounter, zerocounter);
	
	// REMOVE ZEROS
	
	mat_t *gnoz = malloc(sizeof(mat_t));
	
	gnoz->size = G->size;
	gnoz->nnz = G->nnz - zerocounter - diagcounter;
	gnoz->values = calloc(gnoz->nnz, sizeof(double));
	gnoz->cols = calloc(gnoz->nnz, sizeof(unsigned int));
	gnoz->rows = calloc((gnoz->size + 1), sizeof(unsigned int));

	
	for (i = 0; i < G->size; i++){
		for (j = G->rows[i]; j < G->rows[i + 1]; j++){
			if (fabs(G->values[j]) > zero) {
				gnoz->values[counter] = G->values[j];
				gnoz->cols[counter] = G->cols[j];
				counter++;
			}
		}
		gnoz->rows[i + 1] = counter;
	}

	
	printf("%u -> %u\n", pattern->nnz, expanded_patt->nnz);
	printf("%u -> ^%.2lf%%\n\n", gnoz->nnz, 100.0 * (((double) gnoz->nnz - (double) pattern->nnz)/ (double) pattern->nnz));
	
// 	G->nnz = gnoz->nnz;
// 
// 	G->values = realloc(G->values, G->nnz*sizeof(double));
// 	G->cols = realloc(G->cols, G->nnz*sizeof(unsigned int));
// 	
// 	for (i = 0; i < G->nnz; i++){
// 		G->values[i] = gnoz->values[i];
// 		G->cols[i] = gnoz->cols[i];
// 	}
// 	
// 	for (i = 0; i < G->size; i++){
// 		G->rows[i] = gnoz->rows[i];
// 	}
	
	// 	4.
	// 	Drop small entries in G and rescale
	//
	free(diag); free(pattern); free(expanded_patt); free(gnoz);
}




void transpose(mat_t *G, mat_t *Gtransp){

	unsigned int i = 0, j = 0, k = 0;
	unsigned int *contc = calloc(G->size, sizeof(unsigned int));
	
	Gtransp->nnz = G->nnz;
	Gtransp->values = calloc(Gtransp->nnz, sizeof(double));
	Gtransp->cols = calloc(Gtransp->nnz, sizeof(unsigned int));
	Gtransp->rows = calloc((Gtransp->size + 1), sizeof(unsigned int));
	
	for (i = 0; i < G->nnz; i++){contc[G->cols[i]] += 1;}
	for (i = 1; i < G->size + 1; i++){
		Gtransp->rows[i] = Gtransp->rows[i - 1] + contc[i - 1];
	}
	for (i = 0; i < G->size; i++) contc[i] = 0;
	for (i = 0; i < G->size; i++){
		for (j = G->rows[i]; j < G->rows[i + 1]; j++){
			Gtransp->cols[Gtransp->rows[G->cols[k]] + contc[G->cols[k]]] = i;
			Gtransp->values[Gtransp->rows[G->cols[k]] + contc[G->cols[k]]] = G->values[k];
			contc[G->cols[k]] += 1;
			k++;
		}
	}
	free(contc);
}



void transposeCSC(mat_t *G, mat_t *Gtransp){

	unsigned int i = 0;
	
	Gtransp->nnz = G->nnz;
	Gtransp->values = calloc(Gtransp->nnz, sizeof(double));
	Gtransp->rows = calloc(Gtransp->nnz, sizeof(unsigned int));
	Gtransp->cols = calloc((Gtransp->size + 1), sizeof(unsigned int));

	for (i = 0; i < Gtransp->nnz; i++){
		Gtransp->values = G->values;
		Gtransp->rows = G->cols;
	}
	for (i = 0; i < Gtransp->size + 1; i++){
		Gtransp->cols = G->rows;
	}
	
}






void CSCtoCOO(mat_t *Gtransp, mat_t *GtranspCOO, unsigned int *limits, unsigned int nthreads){
	
	unsigned int i, j, k;
	unsigned int counter = 0;
	unsigned int sum = 0;
	unsigned int *contr = calloc(Gtransp->size, sizeof(unsigned int));
	unsigned int *rowlim = calloc(nthreads + 1, sizeof(unsigned int));
	unsigned int *counterrows = calloc(nthreads, sizeof(unsigned int));
	unsigned int elems = (unsigned int) ceil(((double) Gtransp->nnz) / ((double) nthreads));
	
	for (i = 0; i < Gtransp->nnz; i++){contr[Gtransp->rows[i]] += 1;}
	
	for (i = 0; i < Gtransp->size; i++){
		sum += contr[i];
		if ((sum >= elems) || i == (Gtransp->size - 1)){
			counter++;
			rowlim[counter] = i;
			limits[counter] = limits[counter - 1] + sum;
			sum = 0;
		}
	}
	
// 	for (i = 0; i < (nthreads + 1); i++) printf("%u, ", limits[i]);
// 	printf("\n");
// 	
// 	for (i = 0; i < (nthreads + 1); i++) printf("%u, ", rowlim[i]);
// 	printf("\n");
	
	counter = 0;
	
	for (i = 0; i < Gtransp->size; i++){
		for (j = Gtransp->cols[i]; j < Gtransp->cols[i + 1]; j++){
			for (k = 0; k < nthreads; k++){
				if (Gtransp->rows[j] <= rowlim[k + 1]){
					GtranspCOO->values[limits[k] + counterrows[k]] = Gtransp->values[j];
					GtranspCOO->cols[limits[k] + counterrows[k]] = i;
					GtranspCOO->rows[limits[k] + counterrows[k]] = Gtransp->rows[j];
					counterrows[k] += 1;
					break;
				}
			}
		}
	}
	
// 	for (i = 0; i < (nthreads); i++) printf("%u, ", counterrows[i]);
// 	printf("\n");
	
	
	
// 	for (i = 0; i < Gtransp->nnz; i++) printf("%u\t%u\t%lf\t\t\t%lf\n", GtranspCOO->rows[i], GtranspCOO->cols[i], GtranspCOO->values[i], getvalue_matCSC(GtranspCOO->rows[i], GtranspCOO->cols[i], Gtransp));
	
	
	free(counterrows);
	free(rowlim);
	free(contr);
}
	
	
