
#include "cg_aux.h"
#define cachesize 8
#define widthprint 8

double mainval = 16.0;
double sideval = 1.0;

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
extern mat_t *powmat;

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
		patternmode = atoi(argv[7]);
		
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
		
		if ((patternmode == 5) && (argc >= 10)) {
			percentpattern = atoi(argv[8]);
			patternpower = atoi(argv[9]);
		}
		else if ((patternmode == 3) && (argc < 10)) {
			fprintf(stderr, "\nError: Add number of adding sweeps with PCG 5 (e.g.: 1, 2...) and pattern power (e.g.: 2).\n\n");
			return 1;
		}
		
		if ((patternmode == 6) && (argc >= 10)) {
			percentpattern = atoi(argv[8]);
			patternpower = atoi(argv[9]);
		}
		else if ((patternmode == 3) && (argc < 10)) {
			fprintf(stderr, "\nError: Add percentage of elements to add with PCG 6 (e.g.: 100) and pattern power (e.g.: 2).\n\n");
			return 1;
		}
	}
	else if ((mode == 3) && (argc < 8)) {
		patternmode = 0;
	}
	
	if (mode != 3) {patternmode = 0; percentpattern = 0; patternpower = 0;}
	
	
	printf("\n---------------------------------------------------------\n");
	fprintf(stderr, "---------------------------------------------------------\n\n");
	printf("\nMatrix:\t%s\n", matname);
	fprintf(stderr, "Matrix:\t%s\n\n", matname);
	printf("Imax:\t%d\n", imax);
	printf("Error:\t%.2e\n", err);
	printf("\n---------------------------------------------------------\n\n");
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
		
// 		while(fscanf(fp, "%lf", &num) == 1){
// 			
// 			if (count == 0){
// 				mat->size = (unsigned int) num;
// 			}
// 			else if (count == 1){}
// 			else if (count == 2){
// 				mat->nnz = (unsigned int) num;
// 				mat->values = malloc(mat->nnz*sizeof(double));
// 				mat->cols = malloc(mat->nnz*sizeof(unsigned int));
// 				mat->rows = malloc((mat->size + 1) * sizeof(unsigned int));
// 				mat->rows[0] = 0;
// 			}
// 			else if (count == 4){
// 				mat->cols[((i - 1)/3) - 1] = (unsigned int) num - 1;
// 			}
// 			else if (count == 3){
// 				rowvalue = num - 1;
// 				if ((rowvalue-cmp) > 1){
// 					for (j = cmp; j < rowvalue; j++){
// 						position++;
// 						mat->rows[position] = acum;
// 					}
// 					cmp = rowvalue;
// 					acum++;
// 				}
// 				else if ((rowvalue - cmp) == 1){
// 					position++;
// 					cmp = rowvalue;
// 					mat->rows[position] = acum;
// 					acum++;
// 				}
// 				else{acum++;}
// 			}
// 			else if (count == 5){
// 				mat->values[(i - 2)/3 - 1] = num;
// 			}
// 			else {}
// 			i++;
// 			count++;
// 			if (count == 6) count = 3;
// 		}
		
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
			else if (count == 3){
				mat->cols[i/3 - 1] = (unsigned int) num - 1;
			}
			else if (count == 4){
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
	
	fprintf(stderr, "INITIAL NNZ: %i\tLT: %i\n", mat->nnz, (mat->nnz-mat->size)/2 + mat->size);
	return 0;
}


int readrhs(double *b, char* matname, unsigned int size){
	
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
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


void printmat(mat_t *A){

	unsigned int lim = A->size;
	int place = 0;
	
	if (lim < 10){
		printf("\t");
		for (int i = 0; i < lim; i++){
			printf("%*i", widthprint, i);
		}
		printf("\n");
		printf("\n");
		
		for (int i = 0; i < lim; i++){
			printf("%*i", widthprint, i);
			place = 0;
			for (int j = 0; j < lim; j++){
				if ((j == A->cols[A->rows[i] + place]) && (place < (A->rows[i + 1] - A->rows[i]))){
					printf("%*.1e", widthprint, A->values[A->rows[i] + place]);
					place++;
				}
				else {
					printf("\t");
				}
			}
			printf("\n");
		}
		
		printf("\n");
		printf("\n");
	}
	else {
		lim = 10;
		printf("\t");
		for (int i = 0; i < lim; i++){
			printf("%*i", widthprint, i);
		}
		printf("\n");
		printf("\n");
		
		for (int i = 0; i < lim; i++){
			printf("%*i", widthprint, i);
			place = 0;
			for (int j = 0; j < lim; j++){
				if ((j == A->cols[A->rows[i] + place]) && (place < (A->rows[i + 1] - A->rows[i]))){
					printf("%*.1e", widthprint, A->values[A->rows[i] + place]);
					place++;
				}
				else {
					printf("\t");
				}
			}
			printf("\n");
		}
		
		printf("\n");
		printf("\n");
	}
}


void printarr(double *arr, unsigned int dim){
	
	
	if (dim > 10){
		printf("\t");
		for (int j = 0; j < 5; j++){
			printf("%*.2lf", widthprint, arr[j]);
		}
		printf("\n\t");
		for (int j = dim/2; j < dim/2 + 5; j++){
			printf("%*.2lf", widthprint, arr[j]);
		}
		printf("\n\t");
		for (int j = dim - 5; j < dim; j++){
			printf("%*.2lf", widthprint, arr[j]);
		}
		printf("\n\n");
	}
	else {
		printf("\t");
		for (int j = 0; j < dim; j++){
			printf("%*.2lf", widthprint, arr[j]);
		}
		printf("\n\n");
	}
	
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
			fprintf(stderr, "\tPOWER %u OF A\n", patternpower);
			powerA(mat, pattern, expanded_patt, xfinal);
			break;
		case 4:
			printf("PCG:\tParallel ADD PERCENTAGE -> %u\n\n", percentpattern);
			perclt(mat, pattern, expanded_patt, xfinal);
			break;
		case 5:
			printf("PCG:\tOPTIMIZED POWER %u OF A; SWEEPS: %i\n\n", patternpower, percentpattern);
			fprintf(stderr, "\tOPTIMIZED POWER %u OF A; SWEEPS: %i\n", patternpower, percentpattern);
			OptpowerA(mat, pattern, expanded_patt, xfinal);
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
	double *errorrow = calloc(mat->size, sizeof(double));
	double *residrow = calloc(mat->size, sizeof(double));
	double *errrow = calloc(mat->size, sizeof(double));
	
	#pragma omp parallel for private(j, k) reduction(+:totalresid)
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
		
// 		double *atransp = calloc(rowelems*rowelems, sizeof(double));
		
		for (j = expanded_patt->rows[i]; j < expanded_patt->rows[i + 1]; j++){
			arrayindex[j - expanded_patt->rows[i]] = expanded_patt->cols[j];
		}
		
		unsigned int poscount = 0;
		
		/****************************************************************/
		
		
		
		for (j = 0; j < rowelems; j++){					// Rows construct_mat
			for (k = 0; k < rowelems; k++){				// Cols construct_mat
				a[poscount] = getvalue_mat(arrayindex[j], arrayindex[k], mat);
				a_res[poscount] = a[poscount];
				poscount++;
			}
			b[j] = 0.0;
			b_res[j] = 0.0;
		}
		b[rowelems - 1] = 1.0;
		b_res[rowelems - 1] = 1.0;
		

		
// 		for (j = 0; j < rowelems; j++){					// Rows construct_mat
// 			for (k = 0; k < rowelems; k++){				// Cols construct_mat
// 				a[poscount] = getvalue_mat(arrayindex[j], arrayindex[k], powmat);
// 				a_res[poscount] = a[poscount];
// 				poscount++;
// 			}
// 			b[j] = getvalue_mat(arrayindex[j], arrayindex[rowelems - 1], mat);
// 			b_res[j] = getvalue_mat(arrayindex[j], arrayindex[rowelems - 1], mat);
// 		}
		
		
		/****************************************************************/
		
		
// 		if ((i >= 322) && (i < 323)){
// 			printf("\n\n\t");
// 			for (j = 0; j < rowelems; j++){printf("%i\t", arrayindex[j]);}
// 			printf("\n\n");
// 			
// 			for (j = 0; j < rowelems; j++){					// Rows construct_mat
// 				
// 				printf("%i\t", arrayindex[j]);
// 				for (k = 0; k < rowelems; k++){				// Cols construct_mat
// 					printf("%.4lf\t", a[rowelems*j + k]);
// 				}
// 				printf("\n");
// 			}
// 		}
// 		
// 		
// 		if ((i >= 322) && (i < 323)){
// 			printf("\n\n\t");
// 			for (j = 0; j < rowelems; j++) printf("%.4lf\t", b[j]);
// 			printf("\n");
// 		}
		
		dgesv( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
		
// 		if ((i >= 322) && (i < 323)){
// 			printf("%i\t", i);
// 			for (j = 0; j < rowelems; j++) printf("%.4lf\t", b[j]);
// 			printf("\n\n");
// 		}
// 		
// 		
// 		if ((i >= 322) && (i < 323)){
// 			for (j = mat->rows[i]; j < mat->rows[i + 1]; j++){
// 				printf("%i\t", mat->cols[j]);
// 				for (k = mat->rows[mat->cols[j]]; k < mat->rows[mat->cols[j] + 1]; k++){
// 					printf("%i\t", mat->cols[k]);
// 				}
// 				printf("\n");
// 			}
// 			printf("\n\n");
// 		}
		
		
// 		if ((i >= (mat->size - 1)) && (i < mat->size)){
// 			char buf_x[256];
// 			snprintf(buf_x, sizeof buf_x, "../Outputs/cg/DataCG/x.log");
// 			FILE *log = fopen(buf_x, "w");
// 			for ( int j = 0; j < rowelems; j++ ) {
// 				fprintf(log, "%i\t%.8e\n", arrayindex[j], b[j]);
// 			}
// 			fclose(log);
// 		}
		
		
		
		
		
		double norm_x = cblas_ddot(rowelems, b, 1, b, 1);
		
		for (j = 0; j < rowelems; j++){
			for (k = 0; k < rowelems; k++){
				tmp[j] += a_res[j*rowelems + k] * b[k]; 
			}
			r[j] = b_res[j] - tmp[j];
			norm_r += r[j] * r[j];
		}
		
		residrow[i] = sqrt(norm_r);
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
	
	
	printf("\nINVERSE: Accumulated Residual: %.2e\n", totalresid);
		
	double *diag;
	unsigned int diagcounter = 0, zerocounter = 0;
	diag = calloc(G->size, sizeof(double));
	double zero = 1E-10;
	

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
	
	printf("INVERSE: %u zeroes in diag\nINVERSE: %u zeroes in matrix\n", diagcounter, zerocounter);
	
	fprintf(stderr, "INVERSE: %u zeroes in diag\tINVERSE: %u zeroes in matrix\n", diagcounter, zerocounter);
	
	
	// REMOVE ZEROS
	
	int counter = 0;
	
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

	printf("%u -> %u\n", pattern->nnz, counter);
	printf("%u -> %u\n", pattern->nnz, expanded_patt->nnz);
	printf("%u -> ^%.2lf%%\n\n", gnoz->nnz, 100.0 * (((double) gnoz->nnz - (double) pattern->nnz)/ (double) pattern->nnz));
	
	
	fprintf(stderr, "%u -> %u\t%u -> ^%.2lf%%\n\n", pattern->nnz, counter, gnoz->nnz, 100.0 * (((double) gnoz->nnz - (double) pattern->nnz)/ (double) pattern->nnz));
	
	
// 		4.
// 		Drop small entries in G and rescale
	
	
	
	G->nnz = gnoz->nnz;
	G->values = gnoz->values;
	G->cols = gnoz->cols;
	G->rows = gnoz->rows;
	
	
	
	
	
	
	// OTHER STUFF
	
	char buf_res[256];
	snprintf(buf_res, sizeof buf_res, "../Outputs/cg/DataCG/addedelems/res_%s_%i_%i.mtx", matname, percentpattern, patternpower);
	FILE *testres= fopen(buf_res, "w");
	double sum = 0.0;
	int pos = 0;
	double avg = 0.0;
	for(i = 0; i < mat->size; i++){
		sum += residrow[i];
		if ((i%1000 == 0) && (i > 0)){
			avg = sum;
			fprintf(testres, "%i %.20lf\n", pos, avg);
			pos++;
			sum = 0.0;
		}
	}
	fclose(testres);
	
	
	free(diag); free(pattern); free(expanded_patt); 
	free(errorrow); free(residrow); free(errrow);
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
	
	free(counterrows);
	free(rowlim);
	free(contr);
}


// Limits G to cache lines of A
void filterG(mat_t *G, mat_t *A, unsigned int dim, double *r){
	
	int i = 0, j = 0, k = 0, l = 0;
	int cachecol = 0;
	long long cacheline = 0;
	unsigned int testcol = 0;
	int counter = 0;
	
	
	mat_t *Gfinal = malloc(sizeof(mat_t));
	Gfinal->size = dim;
	Gfinal->nnz = G->nnz;
	Gfinal->rows = calloc(dim + 1, sizeof(unsigned int));
	Gfinal->cols = calloc(Gfinal->nnz, sizeof(unsigned int));
	
	printf("\n");
	
	for (i = dim - 10; i < dim; i++){
		printf("%i\t", G->rows[i + 1]);
		for (j = G->rows[i]; j < G->rows[i + 1]; j++){
			printf("%i\t", G->cols[j]);
		}
		printf("\n");
	}
	
	printf("\n");
	
	for (i = 0; i < dim; i++){
		for(j = A->rows[i]; j < A->rows[i + 1]; j++){
			
			cachecol = (int) A->cols[j];
			cacheline = (((long long) &r[A->cols[j]]%64)/8);
			
			for (k = 0; k < cachesize; k++){
				
				testcol = cachecol - (int) cacheline + k;
				
				for (l = G->rows[i]; l < G->rows[i + 1]; l++){
					if (G->cols[l] == testcol){
						Gfinal->cols[counter] = testcol;
						Gfinal->rows[i + 1]++;
						counter++;
						break;
					}
				}
				
				for (l = A->rows[i]; l < A->rows[i + 1]; l++){
					if (((int) A->cols[l] == testcol) && (testcol != cachecol)){
						j++;
						break;
					}
				}
			}
		}
	}
	
	
	for (i = 1; i < dim + 1; i++){
		Gfinal->rows[i] += Gfinal->rows[i - 1];
	}
	
	Gfinal->nnz = counter;
	Gfinal->cols = realloc(Gfinal->cols, Gfinal->nnz*sizeof(unsigned int));
	
	printf("Gvalues: %i\nGfinalvalues: %i\n\n", G->nnz, counter);
	
	
	
	printf("\n");
	for (i = dim - 10; i < dim; i++){
		printf("%i\t", Gfinal->rows[i + 1]);
		for (j = Gfinal->rows[i]; j < Gfinal->rows[i + 1]; j++){
			printf("%i\t", Gfinal->cols[j]);
		}
		printf("\n");
	}
	printf("\n");
	
	G->nnz = Gfinal->nnz;
	G->cols = Gfinal->cols;
	G->rows = Gfinal->rows;
	
	
	
	for (i = dim - 10; i < dim; i++){
		printf("%i\t", G->rows[i + 1]);
		for (j = G->rows[i]; j < G->rows[i + 1]; j++){
			printf("%i\t", G->cols[j]);
		}
		printf("\n");
	}
	
	printf("\n");
	
// 	free(Gfinal);
}



int sumupdown(mat_t *A, int *col, int ract){
	
	int i = 0; 
	unsigned int dim = A->size;
	int row = 0;
	int elements = 1;
	int counter = 0;
	
	row = ract;
	for (i = A->rows[row]; i < A->rows[row + 1]; i++){
		col[counter] = A->cols[i];
		counter++;
	}
	
	if (ract > 0){
		row = ract - 1;
		for (i = A->rows[row]; i < A->rows[row + 1]; i++){
			col[counter] = A->cols[i];
			counter++;
		}
	}
	
	if (ract < (dim - 1)){
		row = ract + 1;
		for (i = A->rows[row]; i < A->rows[row + 1]; i++){
			col[counter] = A->cols[i];
			counter++;
		}
	}
	
	qsort(col, counter, sizeof(int), comp);
	
	unsigned int *tmpcol = calloc(counter, sizeof(int));
	tmpcol[0] = col[0];
	
	for (i = 1; i < counter; i++){
		if (col[i] != col[i - 1]){
			tmpcol[elements] = col[i];
			elements++;
		}
	}
	
	col = realloc(col, elements*sizeof(int));
	for (i = 0; i < elements; i++){
		col[i] = tmpcol[i];
	}
	
	
	free(tmpcol);
	
	return elements;
	
}






/* ------------------------- LLt Filter -------------------------- */
/* ------------------------- LALtx = Lb -------------------------- */


int sumLLt(mat_t *A, int *col, int ract){
	
	int i = 0; 
	unsigned int dim = A->size;
	int row = 0;
	int elements = 1;
	int counter = 0;
	
	row = ract;
	for (i = A->rows[row]; i < A->rows[row + 1]; i++){
		col[counter] = A->cols[i];
		counter++;
		if (A->cols[i] > 0){
			col[counter] = A->cols[i] - 1;
			counter++;
		}
	}
	
	if (ract < (dim - 1)){
		row = ract + 1;
		for (i = A->rows[row]; i < A->rows[row + 1]; i++){
			col[counter] = A->cols[i];
			counter++;
			if (A->cols[i] > 0){
				col[counter] = A->cols[i] - 1;
				counter++;
			}
		}
	}
	
	qsort(col, counter, sizeof(int), comp);
	
	unsigned int *tmpcol = calloc(counter, sizeof(int));
	tmpcol[0] = col[0];
	
	for (i = 1; i < counter; i++){
		if (col[i] != col[i - 1]){
			tmpcol[elements] = col[i];
			elements++;
		}
	}
	
	col = realloc(col, elements*sizeof(int));
	for (i = 0; i < elements; i++){
		col[i] = tmpcol[i];
	}
	
	free(tmpcol);
	
	return elements;
	
}


void smoothmatLLt(mat_t *A){
	
	unsigned int dim = A->size;
	
	unsigned int maxAsizemult = (A->nnz/A->size + 1) * 10;
	double *storevals = calloc(maxAsizemult*A->nnz, sizeof(double));
	int *storecols = calloc(maxAsizemult*A->nnz, sizeof(int));
	int *storerows = calloc(dim + 1, sizeof(int));
	unsigned int pos = 0;
	int i = 0, j = 0;
	int counter = 0;
	
	for (i = 0; i < maxAsizemult*A->nnz; i++) storecols[i] = -10;
	
	#pragma omp parallel for private(j, counter)
	for (i = 0; i < dim; i++){
		
		counter = 0;
		
		int *col = calloc(maxAsizemult, sizeof(int));
		
		counter = sumLLt(A, col, i);
		
		if (i == (dim - 1)){
			
			for (j = 0; j < counter; j++){
				
				if (col[j] < (dim - 1)){
					storevals[maxAsizemult*i + j] = (mainval*mainval*getvalue_mat(i, col[j], A) + sideval*mainval*getvalue_mat(i, col[j] + 1, A));
					storecols[maxAsizemult*i + j] = col[j];
					storerows[i + 1]++;
				}
				else {
					storevals[maxAsizemult*i + j] = mainval*mainval*getvalue_mat(i, col[j], A);
					storecols[maxAsizemult*i + j] = col[j];
					storerows[i + 1]++;
				}
				
			}
			
		}
		else{
			
			for (j = 0; j < counter; j++){
				
				if (col[j] < (dim - 1)){
					storevals[maxAsizemult*i + j] = (mainval*mainval*getvalue_mat(i, col[j], A) + sideval*mainval*getvalue_mat(i, col[j] + 1, A) + sideval*mainval*getvalue_mat(i + 1, col[j], A) + sideval*sideval*getvalue_mat(i + 1, col[j] + 1, A));
					storecols[maxAsizemult*i + j] = col[j];
					storerows[i + 1]++;
				}
				else {
					storevals[maxAsizemult*i + j] = (mainval*mainval*getvalue_mat(i, col[j], A) + sideval*mainval*getvalue_mat(i + 1, col[j], A));
					storecols[maxAsizemult*i + j] = col[j];
					storerows[i + 1]++;
				}
				
			}
			
		}
			
		free(col);
	}
	
	
	
	for (i = 0; i < maxAsizemult*A->nnz; i++){
		if (storecols[i] > -1){
			storecols[pos] = storecols[i];
			storevals[pos] = storevals[i];
			pos++;
		}
	}
	for (i = 1; i < (dim + 1); i++) storerows[i] += storerows[i - 1];
	
	
	
	storecols = realloc(storecols, pos*sizeof(int));
	storevals = realloc(storevals, pos*sizeof(double));
	
	
	A->values = realloc(A->values, pos*sizeof(double));
	A->cols = realloc(A->cols, pos*sizeof(unsigned int));
	
	
	for (i = 0; i < pos; i++){
		A->values[i] = storevals[i];
		A->cols[i] = (unsigned int) storecols[i];
	}
	
	for (i = 0; i < (dim + 1); i++) A->rows[i] = (unsigned int) storerows[i];
	
	A->nnz = pos;
	
	
	for (i = 0; i < dim; i++){
		for (j = A->rows[i]; j < A->rows[i + 1]; j++){
			if (A->cols[j] == i) A->values[j] = A->values[j];
		}
	}	
	
	
	printmat(A);
	
	
	printf("Checking for non-symmetric points\n\n");
	for (i = 0; i < dim; i++){
		for (j = A->rows[i]; j < A->rows[i + 1]; j++){
			if (i >= A->cols[j]){
				if (fabs((getvalue_mat(i, A->cols[j], A) - getvalue_mat(A->cols[j], i, A))) > 1E-6){
					printf("%i, %i\t", i, A->cols[j]);
					printf("%*.4lf ", (int) 8, (getvalue_mat(i, A->cols[j], A)));
					printf("%*.4lf\n", (int) 8, (getvalue_mat(A->cols[j], i, A)));
				}
			}
		}
	}
	printf("\nEnd\n\n");
	
	
	
	free(storecols); free(storerows); free(storevals);
}



void smoothTOmatLLt(mat_t *A, mat_t *Af){
	
	unsigned int dim = A->size;
	
	unsigned int maxAsizemult = (A->nnz/A->size + 1) * 10;
	double *storevals = calloc(maxAsizemult*A->nnz, sizeof(double));
	int *storecols = calloc(maxAsizemult*A->nnz, sizeof(int));
	int *storerows = calloc(dim + 1, sizeof(int));
	unsigned int pos = 0;
	int i = 0, j = 0;
	int counter = 0;
	
	for (i = 0; i < maxAsizemult*A->nnz; i++) storecols[i] = -10;
	
	#pragma omp parallel for private(j, counter)
	for (i = 0; i < dim; i++){
		
		counter = 0;
		
		int *col = calloc(maxAsizemult, sizeof(int));
		
		counter = sumLLt(A, col, i);
		
		if (i == (dim - 1)){
			
			for (j = 0; j < counter; j++){
				
				if (col[j] < (dim - 1)){
					storevals[maxAsizemult*i + j] = (mainval*mainval*getvalue_mat(i, col[j], A) + sideval*mainval*getvalue_mat(i, col[j] + 1, A));
					storecols[maxAsizemult*i + j] = col[j];
					storerows[i + 1]++;
				}
				else {
					storevals[maxAsizemult*i + j] = mainval*mainval*getvalue_mat(i, col[j], A);
					storecols[maxAsizemult*i + j] = col[j];
					storerows[i + 1]++;
				}
				
			}
			
		}
		else{
			
			for (j = 0; j < counter; j++){
				
				if (col[j] < (dim - 1)){
					storevals[maxAsizemult*i + j] = (mainval*mainval*getvalue_mat(i, col[j], A) + sideval*mainval*getvalue_mat(i, col[j] + 1, A) + sideval*mainval*getvalue_mat(i + 1, col[j], A) + sideval*sideval*getvalue_mat(i + 1, col[j] + 1, A));
					storecols[maxAsizemult*i + j] = col[j];
					storerows[i + 1]++;
				}
				else {
					storevals[maxAsizemult*i + j] = (mainval*mainval*getvalue_mat(i, col[j], A) + sideval*mainval*getvalue_mat(i + 1, col[j], A));
					storecols[maxAsizemult*i + j] = col[j];
					storerows[i + 1]++;
				}
				
			}
			
		}
			
		free(col);
	}
	
	
	
	for (i = 0; i < maxAsizemult*A->nnz; i++){
		if (storecols[i] > -1){
			storecols[pos] = storecols[i];
			storevals[pos] = storevals[i];
			pos++;
		}
	}
	for (i = 1; i < (dim + 1); i++) storerows[i] += storerows[i - 1];
	
	
	
	storecols = realloc(storecols, pos*sizeof(int));
	storevals = realloc(storevals, pos*sizeof(double));
	
	
	Af->values = calloc(pos, sizeof(double));
	Af->cols = calloc(pos, sizeof(unsigned int));
	Af->rows = calloc(dim + 1, sizeof(unsigned int));
	Af->size = dim;
	
	for (i = 0; i < pos; i++){
		Af->values[i] = storevals[i];
		Af->cols[i] = (unsigned int) storecols[i];
	}
	
	for (i = 0; i < (dim + 1); i++) Af->rows[i] = (unsigned int) storerows[i];
	
	Af->nnz = pos;
	
	
	printf("Checking for non-symmetric points\n\n");
	for (i = 0; i < dim; i++){
		for (j = Af->rows[i]; j < Af->rows[i + 1]; j++){
			if (i >= Af->cols[j]){
				if (fabs((getvalue_mat(i, Af->cols[j], Af) - getvalue_mat(Af->cols[j], i, Af))) > 1E-6){
					printf("%i, %i\t", i, Af->cols[j]);
					printf("%*.4lf ", (int) 8, (getvalue_mat(i, Af->cols[j], Af)));
					printf("%*.4lf\n", (int) 8, (getvalue_mat(Af->cols[j], i, Af)));
				}
			}
		}
	}
	printf("\nEnd\n\n");
	
	
	
	free(storecols); free(storerows); free(storevals);
}




/* ------------------------- L Filter -------------------------- */
/* ------------------------- LAx = Lb -------------------------- */


int sumL(mat_t *A, int *col, int ract){
	
	int i = 0; 
	unsigned int dim = A->size;
	int row = 0;
	int elements = 1;
	int counter = 0;
	
	row = ract;
	for (i = A->rows[row]; i < A->rows[row + 1]; i++){
		col[counter] = A->cols[i];
		counter++;
	}
	
	if (ract < (dim - 1)){
		row = ract + 1;
		for (i = A->rows[row]; i < A->rows[row + 1]; i++){
			col[counter] = A->cols[i];
			counter++;
		}
	}
	
	qsort(col, counter, sizeof(int), comp);
	
	unsigned int *tmpcol = calloc(counter, sizeof(int));
	tmpcol[0] = col[0];
	
	for (i = 1; i < counter; i++){
		if (col[i] != col[i - 1]){
			tmpcol[elements] = col[i];
			elements++;
		}
	}
	
	col = realloc(col, elements*sizeof(int));
	for (i = 0; i < elements; i++){
		col[i] = tmpcol[i];
	}
	
	free(tmpcol);
	
	return elements;
	
}



void smoothmatL(mat_t *A){
	
	unsigned int dim = A->size;
	
	unsigned int maxAsizemult = (A->nnz/A->size + 1) * 10;
	double *storevals = calloc(maxAsizemult*A->nnz, sizeof(double));
	int *storecols = calloc(maxAsizemult*A->nnz, sizeof(int));
	int *storerows = calloc(dim + 1, sizeof(int));
	unsigned int pos = 0;
	int i = 0, j = 0;
	int counter = 0;
	double zero = 1E-6;
	
	for (i = 0; i < maxAsizemult*A->nnz; i++) storecols[i] = -10;
	
	#pragma omp parallel for private(j, counter)
	for (i = 0; i < dim; i++){
		
		counter = 0;
		
		int *col = calloc(maxAsizemult, sizeof(int));
		
		counter = sumL(A, col, i);
		
		if (i == (dim - 1)){
			
			for (j = 0; j < counter; j++){
				
				storevals[maxAsizemult*i + j] = mainval*getvalue_mat(i, col[j], A);
				storecols[maxAsizemult*i + j] = col[j];
				storerows[i + 1]++;

			}
			
		}
		else{
			
			for (j = 0; j < counter; j++){
				
				storevals[maxAsizemult*i + j] = (mainval*getvalue_mat(i, col[j], A) + sideval*getvalue_mat(i + 1, col[j], A));
				storecols[maxAsizemult*i + j] = col[j];
				storerows[i + 1]++;
				
			}
			
		}
			
		free(col);
	}
	
	
	
	for (i = 0; i < maxAsizemult*A->nnz; i++){
		if (storecols[i] > -1){
			storecols[pos] = storecols[i];
			storevals[pos] = storevals[i];
			pos++;
		}
	}
	for (i = 1; i < (dim + 1); i++) storerows[i] += storerows[i - 1];
	
	
	
	storecols = realloc(storecols, pos*sizeof(int));
	storevals = realloc(storevals, pos*sizeof(double));
	
	
	A->values = realloc(A->values, pos*sizeof(double));
	A->cols = realloc(A->cols, pos*sizeof(unsigned int));
	
	
	for (i = 0; i < pos; i++){
		A->values[i] = storevals[i];
		A->cols[i] = (unsigned int) storecols[i];
	}
	
	for (i = 0; i < (dim + 1); i++) A->rows[i] = (unsigned int) storerows[i];
	
	A->nnz = pos;
	
	
	pos = 0;
	for (i = 0; i < A->nnz; i++){
		if (fabs(A->values[i]) < zero)
			pos++;
	}
	
	fprintf(stderr, "\n\nFiltered zeroes: %i\n\n", pos);
	
	
	for (i = 0; i < dim; i++){
		for (j = A->rows[i]; j < A->rows[i + 1]; j++){
			if (A->cols[j] == i) A->values[j] = A->values[j];
		}
	}
	
	
	printmat(A);
	
	
	printf("Checking for non-symmetric points\n\n");
	for (i = 0; i < dim; i++){
		for (j = A->rows[i]; j < A->rows[i + 1]; j++){
			if (i >= A->cols[j]){
				if (fabs((getvalue_mat(i, A->cols[j], A) - getvalue_mat(A->cols[j], i, A))) > 1E-6){
					printf("%i, %i\t", i, A->cols[j]);
					printf("%*.4lf ", (int) 8, (getvalue_mat(i, A->cols[j], A)));
					printf("%*.4lf\n", (int) 8, (getvalue_mat(A->cols[j], i, A)));
				}
			}
		}
	}
	printf("\nEnd\n\n");
	
	
	
	free(storecols); free(storerows); free(storevals);
}






void smoothTOmatL(mat_t *A, mat_t *Af){
	
	unsigned int dim = A->size;
	
	unsigned int maxAsizemult = (A->nnz/A->size + 1) * 10;
	double *storevals = calloc(maxAsizemult*A->nnz, sizeof(double));
	int *storecols = calloc(maxAsizemult*A->nnz, sizeof(int));
	int *storerows = calloc(dim + 1, sizeof(int));
	unsigned int pos = 0;
	int i = 0, j = 0;
	int counter = 0;
	
	for (i = 0; i < maxAsizemult*A->nnz; i++) storecols[i] = -10;
	
	#pragma omp parallel for private(j, counter)
	for (i = 0; i < dim; i++){
		
		counter = 0;
		
		int *col = calloc(maxAsizemult, sizeof(int));
		
		counter = sumL(A, col, i);
		
		if (i == (dim - 1)){
			
			for (j = 0; j < counter; j++){
				
				storevals[maxAsizemult*i + j] = mainval*getvalue_mat(i, col[j], A);
				storecols[maxAsizemult*i + j] = col[j];
				storerows[i + 1]++;

			}
			
		}
		else{
			
			for (j = 0; j < counter; j++){
				
				storevals[maxAsizemult*i + j] = (mainval*getvalue_mat(i, col[j], A) + sideval*getvalue_mat(i + 1, col[j], A));
				storecols[maxAsizemult*i + j] = col[j];
				storerows[i + 1]++;
				
			}
			
		}
			
		free(col);
	}
	
	
	
	for (i = 0; i < maxAsizemult*A->nnz; i++){
		if (storecols[i] > -1){
			storecols[pos] = storecols[i];
			storevals[pos] = storevals[i];
			pos++;
		}
	}
	for (i = 1; i < (dim + 1); i++) storerows[i] += storerows[i - 1];
	
	
	
	storecols = realloc(storecols, pos*sizeof(int));
	storevals = realloc(storevals, pos*sizeof(double));
	
	
	Af->values = calloc(pos, sizeof(double));
	Af->cols = calloc(pos, sizeof(unsigned int));
	Af->rows = calloc(dim + 1, sizeof(unsigned int));
	Af->size = dim;
	
	for (i = 0; i < pos; i++){
		Af->values[i] = storevals[i];
		Af->cols[i] = (unsigned int) storecols[i];
	}
	
	for (i = 0; i < (dim + 1); i++) Af->rows[i] = (unsigned int) storerows[i];
	
	Af->nnz = pos;
	
	
	printf("Checking for non-symmetric points\n\n");
	for (i = 0; i < dim; i++){
		for (j = Af->rows[i]; j < Af->rows[i + 1]; j++){
			if (i >= Af->cols[j]){
				if (fabs((getvalue_mat(i, Af->cols[j], Af) - getvalue_mat(Af->cols[j], i, Af))) > 1E-6){
					printf("%i, %i\t", i, Af->cols[j]);
					printf("%*.4lf ", (int) 8, (getvalue_mat(i, Af->cols[j], Af)));
					printf("%*.4lf\n", (int) 8, (getvalue_mat(Af->cols[j], i, Af)));
				}
			}
		}
	}
	printf("\nEnd\n\n");
	
	
	
	free(storecols); free(storerows); free(storevals);
}








/* ------------------------- Array Filter -------------------------- */


void smootharrL(double *b, unsigned int dim){
	
	double *storevals = calloc(dim, sizeof(double));
	
	int i = 0;
	
	#pragma omp parallel for
	for (i = 0; i < dim; i++){
		
		if (i == (dim - 1)){
			storevals[i] = mainval*b[i];
		}
		else{
			storevals[i] = (mainval*b[i] + sideval*b[i + 1]);
		}
			
	}
	
	for (i = 0; i < dim; i++) b[i] = storevals[i];
	
	
	free(storevals);
}


void backsmootharrLt(double *b, unsigned int dim){
	
	double *storevals = calloc(dim, sizeof(double));
	
	int i = 0;
	
	#pragma omp parallel for
	for (i = 0; i < dim; i++){
		
		if (i == 0){
			storevals[i] = mainval*b[i];
		}
		else{
			storevals[i] = (mainval*b[i] + sideval*b[i - 1]);
		}
			
	}
	
	for (i = 0; i < dim; i++) b[i] = storevals[i];
	
	
	free(storevals);
}




void backsmootharrLtV2(double *b, unsigned int dim, unsigned int *arrayindex){
	
	double *storevals = calloc(dim, sizeof(double));
	
	int i = 0;
	
	#pragma omp parallel for
	for (i = 0; i < dim; i++){
		
		if (i == 0){
			storevals[i] = mainval*b[i];
		}
		else{
			
			if ((arrayindex[i] - arrayindex[i - 1]) == 1){
				storevals[i] = (mainval*b[i] + sideval*b[i - 1]);
			}
			else {
				storevals[i] = mainval*b[i];
			}
		}
			
	}
	
	for (i = 0; i < dim; i++) b[i] = storevals[i];
	
	
	free(storevals);
}





void solveLU(mat_t *A, unsigned int dim, double *b){
	
	int i, j;
	int rowelems = (int) dim;
	int n = rowelems, nrhs = 1, lda = rowelems, ldb = rowelems, info;
	int ipiv[rowelems];
	double *a = malloc(rowelems*rowelems*sizeof(double));
	int posact = 0;
	int poscheck = 0;
	
	/****************************************************************/
	
	for (i = 0; i < dim; i++){
		for (j = 0; j < dim; j++){
			if (j == A->cols[poscheck]){
				a[posact] = A->values[poscheck];
				posact++;
				poscheck++;
			}
			else {
				a[posact] = 0.0;
				posact++;
			}
		}
	}
	
	
// 	for (i = 0; i < dim; i++){
// 		for (j = 0; j < dim; j++){
// 			printf("%.4lf\t", a[dim*i + j]);
// 		}
// 		printf("\n");
// 	}

	printf("V1\n");
	printarr(b, dim);
	
	dgesv( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
	
	printf("V2\n");
	printarr(b, dim);
	
// 	printf("INFO: %i\n", info);
}













void impLU(mat_t *A, unsigned int dim, double *b){
	
	int i, j;
	int rowelems = (int) dim;
	double *a = malloc(rowelems*rowelems*sizeof(double));
	int posact = 0;
	int poscheck = 0;
	
	/****************************************************************/
	
	for (i = 0; i < dim; i++){
		for (j = 0; j < dim; j++){
			if (j == A->cols[poscheck]){
				a[posact] = A->values[poscheck];
				posact++;
				poscheck++;
			}
			else {
				a[posact] = 0.0;
				posact++;
			}
		}
	}
	


	int counter = 0;
	double **Ap = malloc(rowelems*sizeof(double*));
	double tol = 1E-16;
	int *P = calloc(rowelems + 1, sizeof(int));
	double *x = calloc(rowelems, sizeof(double));
	
	
	for (j = 0; j < rowelems; j++){
		*(Ap + j) = malloc(rowelems * sizeof(double));
		for (int k = 0; k < rowelems; k++){
			*(*(Ap + j) + k) = a[counter];
			counter++;
		}
	}
	
	LUPDecompose(Ap, rowelems, tol, P);
	LUPSolve(Ap, P, b, rowelems, x);
	
	printf("XV2\n");
	printarr(x, dim);
	
}

































































































void precond_vf(mat_t *mat, mat_t *G, double *xfinal){
	
	int i = 0, j = 0, k = 0;
	
	pat_t *pattern = malloc(sizeof(pat_t));
	pattern->rows = calloc((G->size + 1), sizeof(unsigned int));
	
	pat_t *expanded_patt = malloc(sizeof(pat_t));
	expanded_patt->rows = calloc((G->size + 1), sizeof(unsigned int));
	
	mat_t *Af = malloc(sizeof(mat_t));
	
	if (rhs==1){
		smoothTOmatL(mat, Af);
	}
	else if (rhs == 2){
		smoothTOmatLLt(mat, Af);
	}
	else {
		Af = mat;
	}
	
	
	printmat(mat);
	printmat(Af);
	
	
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
			fprintf(stderr, "\tPOWER %u OF A\n", patternpower);
			powerA(mat, pattern, expanded_patt, xfinal);
			break;
		case 4:
			printf("PCG:\tADD PERCENTAGE -> %u\n\n", percentpattern);
			perclt(mat, pattern, expanded_patt, xfinal);
			break;
		case 5:
			printf("PCG:\tOPTIMIZED POWER %u OF A; SWEEPS: %i\n\n", patternpower, percentpattern);
			fprintf(stderr, "\tOPTIMIZED POWER %u OF A; SWEEPS: %i\n", patternpower, percentpattern);
			OptpowerA(mat, pattern, expanded_patt, xfinal);
			break;
		case 6:
			printf("PCG:\tPATTERN OF AF\n\n");
			fprintf(stderr, "\tPATTERN OF AF\n");
			powerAf(Af, pattern, expanded_patt, xfinal);
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
	double *errorrow = calloc(mat->size, sizeof(double));
	double *residrow = calloc(mat->size, sizeof(double));
	double *errrow = calloc(mat->size, sizeof(double));
	
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
		
		/****************************************************************/
		
		
		
		for (j = 0; j < rowelems; j++){					// Rows construct_mat
			for (k = 0; k < rowelems; k++){				// Cols construct_mat
				a[poscount] = getvalue_mat(arrayindex[j], arrayindex[k], Af);
				a_res[poscount] = a[poscount];
				poscount++;
			}
			b[j] = 0.0;
			b_res[j] = 0.0;
		}
		b[rowelems - 1] = 1.0;
		b_res[rowelems - 1] = 1.0;
		
		if (rhs != 0){
			smootharrL(b, rowelems);
			smootharrL(b_res, rowelems);
		}
		
		if ((i >= 32200) && (i < 32201)){
			printf("\n\n\t");
			for (j = 0; j < rowelems; j++){printf("%i\t", arrayindex[j]);}
			printf("\n\n");
			
			for (j = 0; j < rowelems; j++){					// Rows construct_mat
				
				printf("%i\t", arrayindex[j]);
				for (k = 0; k < rowelems; k++){				// Cols construct_mat
					printf("%.4lf\t", a[rowelems*j + k]);
				}
				printf("\n");
			}
		}
		
		
		if ((i >= 32200) && (i < 32201)){
			printf("\n\n\t");
			for (j = 0; j < rowelems; j++) printf("%.4lf\t", b[j]);
			printf("\n");
		}
		
		dgesv( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
		
		if ((i >= 32200) && (i < 32201)){
			printf("%i\t", i);
			for (j = 0; j < rowelems; j++) printf("%.4lf\t", b[j]);
			printf("\n\n");
		}
		
		
		if ((i >= 32200) && (i < 32201)){
			for (j = mat->rows[i]; j < mat->rows[i + 1]; j++){
				printf("%i\t", mat->cols[j]);
				for (k = mat->rows[mat->cols[j]]; k < mat->rows[mat->cols[j] + 1]; k++){
					printf("%i\t", mat->cols[k]);
				}
				printf("\n");
			}
			printf("\n\n");
		}
		
		
		
		double norm_x = cblas_ddot(rowelems, b, 1, b, 1);
		
		for (j = 0; j < rowelems; j++){
			for (k = 0; k < rowelems; k++){
				tmp[j] += a_res[j*rowelems + k] * b[k]; 
			}
			r[j] = b_res[j] - tmp[j];
			norm_r += r[j] * r[j];
		}
		
		residrow[i] = sqrt(norm_r);
		totalresid += sqrt(norm_r)/sqrt(norm_x);
		
		
		if (rhs == 2){
			backsmootharrLtV2(b, rowelems, arrayindex);
		}
		else {
			
		}
		
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
		
		
	}
	
	
	printf("\nINVERSE: Accumulated Residual: %.2e\n", totalresid);
		
		
	double *diag;
	unsigned int diagcounter = 0, zerocounter = 0;
	diag = calloc(G->size, sizeof(double));
	double zero = 1E-10;
	

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
	
	printf("INVERSE: %u zeroes in diag\nINVERSE: %u zeroes in matrix\n", diagcounter, zerocounter);
	
	fprintf(stderr, "INVERSE: %u zeroes in diag\tINVERSE: %u zeroes in matrix\n", diagcounter, zerocounter);
	
	
	// REMOVE ZEROS
	
	int counter = 0;
	
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

	printf("%u -> %u\n", pattern->nnz, counter);
	printf("%u -> %u\n", pattern->nnz, expanded_patt->nnz);
	printf("%u -> ^%.2lf%%\n\n", gnoz->nnz, 100.0 * (((double) gnoz->nnz - (double) pattern->nnz)/ (double) pattern->nnz));
	
	
	fprintf(stderr, "%u -> %u\t%u -> ^%.2lf%%\n\n", pattern->nnz, counter, gnoz->nnz, 100.0 * (((double) gnoz->nnz - (double) pattern->nnz)/ (double) pattern->nnz));
	
	
// 		4.
// 		Drop small entries in G and rescale
	
	
	
	G->nnz = gnoz->nnz;
	G->values = gnoz->values;
	G->cols = gnoz->cols;
	G->rows = gnoz->rows;
	
	
	
	
	
	
	// OTHER STUFF
	
	char buf_res[256];
	snprintf(buf_res, sizeof buf_res, "../Outputs/cg/DataCG/addedelems/res_%s_%i_%i.mtx", matname, percentpattern, patternpower);
	FILE *testres= fopen(buf_res, "w");
	double sum = 0.0;
	int pos = 0;
	double avg = 0.0;
	for(i = 0; i < mat->size; i++){
		sum += residrow[i];
		if ((i%1000 == 0) && (i > 0)){
			avg = sum;
			fprintf(testres, "%i %.20lf\n", pos, avg);
			pos++;
			sum = 0.0;
		}
	}
	fclose(testres);
	
	
	
	free(diag); free(pattern); free(expanded_patt); 
	free(errorrow); free(residrow); free(errrow);
}
