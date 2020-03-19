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
extern int bsize; 
extern FILE *fp;
extern int patternmode;
extern int percentpattern;
extern int patternpower;
extern int reps;
extern mat_t *powmat;


/* 		Function name: cg_config
 * 		Purpose: Read input arguments and set adequate variables.
 * 
 * 		Inputs: Execution arguments
 * 
 * 		Outputs:
 * 			- 0 on success
 * 			- 1 on failure
 * 
 * 		TODO:
 *
 * 		- Rearrange inputs
 * 		- Maybe rethink how to run code
 * 		- Adapt to filtered matrices
 * 
 * 		Others: 
 */

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



/* 		Function name: cg_setup
 * 		Purpose: Read matrix file and store data to a matrix type
 * 
 * 		Inputs: Matrix name, Pointer to matrix (mat_t)
 * 
 * 		Outputs:
 * 			- 0 on success
 * 			- 1 on failure
 * 
 * 		TODO:
 *
 * 		- Adapt to matrices with swapped rows and cols.
 * 		- Optimize code
 * 
 * 		Others: 
 */

int cg_setup(char* matname, mat_t *mat){
	
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int j = 0, i = 0;
	int rowvalue = 0, cmp = 0, acum = 0, position = 0, count = 0;
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
// 				mat->size = (int) num;
// 			}
// 			else if (count == 1){}
// 			else if (count == 2){
// 				mat->nnz = (int) num;
// 				mat->values = malloc(mat->nnz*sizeof(double));
// 				mat->cols = malloc(mat->nnz*sizeof(int));
// 				mat->rows = malloc((mat->size + 1) * sizeof(int));
// 				mat->rows[0] = 0;
// 			}
// 			else if (count == 4){
// 				mat->cols[((i - 1)/3) - 1] = (int) num - 1;
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
				mat->size = (int) num;
			}
			else if (count == 1){}
			else if (count == 2){
				mat->nnz = (int) num;
				mat->values = malloc(mat->nnz*sizeof(double));
				mat->cols = malloc(mat->nnz*sizeof(int));
				mat->rows = malloc((mat->size + 1) * sizeof(int));
				mat->rows[0] = 0;
			}
			else if (count == 3){
				mat->cols[i/3 - 1] = (int) num - 1;
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



/* 		Function name: readrhs
 * 		Purpose: Reads a RHS if provided
 * 
 * 		Inputs: Pointer to array, Matrix name, Matrix size
 * 
 * 		Outputs:
 * 			- 0 on success
 * 			- 1 on failure
 * 
 * 		TODO:
 *
 * 		- Adapt to reading a specific file
 * 		- Optimize code
 * 
 * 		Others: 
 * 
 * 		- Reads only files called RHS_matname
 */

int readrhs(double *b, char* matname, int size){
	
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int index = 0;
	double num = 0;
	char buff[256];
	int count = 0;
	
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
				if ((int) num == size) {printf("\nOK\n");}
				else {return 1;}
			}
			else if (count == 1){}
			else if (count == 2){
				index = (int) num;
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


/* 		Function name: printfigdatax
 * 		Purpose: Prints files with Data and output X
 * 
 * 		Inputs: Matrix name, X, Initial X, all measurement variables
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 		Others: 
 */

void printfigdatax(char *matname, int patternpower, int rhs, mat_t *G, mat_t *A, int percentpattern, int bniter, double btottime, double btimeiter, double btimemult, double btimetransp, long long bdcmnorm, long long bxdcmnorm, long long bdifmult, long long bdcmtransp, long long bxdcmtransp, long long bdiftransp, int brep, double *x, double *initx){
	
	char buffig[256];
	snprintf(buffig, sizeof buffig, "../Outputs/cg/DataCG/Fig_DATA_%s_%i_%i", matname, patternpower, rhs);
	FILE *outmnfig = fopen(buffig, "a");
	if (outmnfig == NULL){
		printf("Could not open writing file.");
	}
	
	fprintf(outmnfig, "%.2lf\t", 100.0 * ((double) G->nnz) / (double) (((A->nnz - A->size)/2) + A->size) - 100.0);
	fprintf(outmnfig, "%u\t", percentpattern);
	fprintf(outmnfig, "%u\t", G->nnz);
	
	fprintf(outmnfig, "%i\t", bniter);
	fprintf(outmnfig, "%lf\t", btottime);
	fprintf(outmnfig, "%lf\t", btimeiter);
	fprintf(outmnfig, "%lf\t", btimemult);
	fprintf(outmnfig, "%lf\t", btimetransp);
	fprintf(outmnfig, "%lli\t", bdcmnorm);
	fprintf(outmnfig, "%lli\t", bxdcmnorm);
	fprintf(outmnfig, "%lli\t", bdifmult);
	fprintf(outmnfig, "%lli\t", bdcmtransp);
	fprintf(outmnfig, "%lli\t", bxdcmtransp);
	fprintf(outmnfig, "%lli\t", bdiftransp);
	
	fprintf(outmnfig, "%i\t", A->nnz);
	fprintf(outmnfig, "%i\t", A->size);
	fprintf(outmnfig, "%i\n", brep);
	
	fclose(outmnfig);
	
	
	char buf_x[256];
	snprintf(buf_x, sizeof buf_x, "../Outputs/cg/DataCG/x_%s.log", matname);
	FILE *logt = fopen(buf_x, "w");
	for ( int i = 0; i < A->size; i++ ) {
		fprintf(logt, "%i\t%.8e\t%.8e\n", i, x[i], initx[i]);
	}
	fclose(logt);
}



/* 		Function name: checkbesttime
 * 		Purpose: Compares best time with repetition time. Updates best time and measurements of best time
 * 
 * 		Inputs: Measurement variables, best time
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- POINTERS!!
 * 
 * 		Others: 
 * 
 * 			- Not available
 */

void checkbesttime(double end_usec, double start_usec, double tmin, int bniter, int i, double btottime, double btimeiter, double btimemult, double btimetransp, long long bdcmnorm, long long bxdcmnorm, long long bdifmult, long long bdcmtransp, long long bxdcmtransp, long long bdiftransp, int brep, int m, double tmult, double ttransp, long long dcmnorm, long long xdcmnorm, long long dcmtransp, long long xdcmtransp){
	
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
}


/* 		Function name: printmat
 * 		Purpose: Prints in stdout a given matrix. If too large prints only 1o first initial rows and cols
 * 
 * 		Inputs: Matrix pointer
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 		Others: 
 */

void printmat(mat_t *A){

	int lim = A->size;
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



/* 		Function name: printarr
 * 		Purpose: Prints a given array. If too large prints 3 small sections
 * 
 * 		Inputs: Array pointer
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 		Others: 
 */

void printarr(double *arr, int dim){
	
	
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



/* 		Function name: papiinit
 * 		Purpose: Initialize papi functions
 * 
 * 		Inputs: retval, EventSet
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 *
 * 		- IFDEF or move to .h
 * 
 * 		Others: 
 */

void papiinit(int retval, int EventSet){
	retval = PAPI_library_init(PAPI_VER_CURRENT);
	if (retval != PAPI_VER_CURRENT){
		fprintf(stderr, "PAPI library init error!\n");
		exit(1);
		
	}
	if (PAPI_thread_init((long unsigned int (*)(void)) omp_get_thread_num) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	if (PAPI_create_eventset(&EventSet) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	if (PAPI_add_event(EventSet, PAPI_L1_DCM) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	if (PAPI_start(EventSet) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
}



/* 		Function name: papiend
 * 		Purpose: Stops papi
 * 
 * 		Inputs: EventSet, Pointer to array
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 *
 * 		Others: 
 */

void papiend(int EventSet, long_long *values){
	if (PAPI_stop(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
}



/* 		Function name: papicount
 * 		Purpose: Gets papi counters for SpMV operations
 * 
 * 		Inputs: Matrix pointers (mat_t), measurement variables, operating array pointers
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 *
 * 		Others: 
 */

void papicount(int dim, mat_t *G, mat_t *GtranspCOO, double *tmp, long long dcmnorm, long long xdcmnorm, long long dcmtransp, long long xdcmtransp, long_long *values, double *dumm, int *limits, int nthreads, int EventSet, double *r, double *s){
	
	pzerovect(tmp, dim);
	PAPI_reset(EventSet);
	pmultMatVect_DUMM(tmp, G);
	if (PAPI_read(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	xdcmnorm = values[0];
	pzerovect(dumm, dim);
	PAPI_reset(EventSet);
	pmultMatVectCOO_DUMM(dumm, GtranspCOO, limits, nthreads);
	if (PAPI_read(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	xdcmtransp = values[0];
	
	pzerovect(tmp, dim);
	PAPI_reset(EventSet);
	pmultMatVect(tmp, r, G);
	if (PAPI_read(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	dcmnorm = values[0];
	pzerovect(s, dim);
	PAPI_reset(EventSet);
	pmultMatVectCOO(s, tmp, GtranspCOO, limits, nthreads);
	if (PAPI_read(EventSet, values) != PAPI_OK) {printf("\nPAPI ERROR!\n");}
	dcmtransp = values[0];
}



/* 		Function name: inverseD
 * 		Purpose: Gets diagonal values of a matrix and inverts them
 * 
 * 		Inputs: Matrix pointer (mat_t)
 * 
 * 		Outputs:
 * 			- Matrix pointer (mat_t)
 * 
 * 		TODO:
 *
 * 		Others: 
 */

mat_t *inverseD(mat_t *mat){
	
	int i = 0;
	mat_t *D;
	
	D = malloc(sizeof(mat_t));
	D->size = mat->size;
	D->nnz = D->size;
	D->values = malloc(D->nnz*sizeof(double));
	D->cols = malloc(D->nnz*sizeof(int));
	D->rows = malloc((D->size + 1)*sizeof(int));
	
	for (i = 0; i < D->size; i++){
		D->values[i] = 1.0/getvalue_mat(i, i, mat);
		D->cols[i] = i;
		D->rows[i + 1] = i + 1;
	}
	return D;
}



/* 		Function name: cg_construct
 * 		Purpose: Solves Ax=b
 * 
 * 		Inputs: Matrix pointer (mat_t), array pointers (x, b), maximum number of iterations, tolerance 
 * 
 * 		Outputs:
 * 			- Pointer to resulting X
 * 
 * 		TODO:
 *
 * 		Others: 
 */

double *cg_construct(mat_t *A, double *x, double *b, int maxiter, double error){

	int i = 0;
	double *d;
	double d_new = 0.0;
	double *q;
	double *r;
	double alpha = 0.0;
	double d_old = 0.0;
	double beta = 0.0;
	double *tmp;
	int dim = A->size;
	
	r = malloc(dim*sizeof(double));
	d = malloc(dim*sizeof(double));
	q = malloc(dim*sizeof(double));
	tmp = malloc(dim*sizeof(double));
	
	pmultMatVect(tmp, x, A);
	psubtVects(b, tmp, r, dim);
	pequalvects(d, r, dim);
	d_new = pmultVectVect(r, r, dim);
	double norm_b = pmultVectVect(b, b, dim);
	double norm_r = sqrt(pmultVectVect(r, r, dim));
	while ((i < maxiter)&&(norm_r/norm_b > error)){
		pmultMatVect(q, d, A);
		alpha = d_new / pmultVectVect(d, q, dim);
		pscaleVect(alpha, d, tmp, dim);
		paddToVect(x, tmp, dim);
		pscaleVect(alpha, q, tmp, dim);
		psubToVect(r, tmp, dim);
		d_old = d_new;
		d_new = pmultVectVect(r, r, dim);
		beta = d_new/d_old;
		pscaleVect(beta, d, tmp, dim);
		paddVects(r, tmp, d, dim);
		norm_r = sqrt(pmultVectVect(r, r, dim));
		i++;
	}
	free(r); free(d); free(q); free(tmp);
	return x;
}



/* 		Function name: LUPDecompose
 * 		Purpose: LU Decomposition
 * 
 * 		Inputs: Pointer of pointers to input matrix, size, tolerance, unit permutation matrix
 * 
 * 		Outputs:
 * 			- 0 on failure
 * 			- 1 on success
 * 
 * 		TODO:
 *
 * 		Others: 
 */

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



/* 		Function name: LUPSolve
 * 		Purpose: Solves Ax=b with LU Decomposed matrix
 * 
 * 		Inputs: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * 
 * 		Outputs:
 * 			- x - solution vector of A*x=b
 * 
 * 		TODO:
 *
 * 		Others: 
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



/* 		Function name: precond
 * 		Purpose: Obtains an approximate inverse G from a given matrix
 * 
 * 		Inputs: Pointer to matrix (mat_t), Pointer to inverse (mat_t), Pointer to X array from general problem
 * 
 * 		Outputs:
 * 			- Inverse matrix of A (G) with specific pattern
 * 
 * 		TODO:
 * 
 * 			- Clean code and maybe redesign
 * 			- Work out way to filter G or not depending on user selection
 *
 * 		Others: Makes use of Lapack library, but could use implemented LUPD solver
 */

void precond(mat_t *mat, mat_t *G, double *xfinal){
	
	int i = 0, j = 0, k = 0;
	
	pat_t *pattern = malloc(sizeof(pat_t));
	pattern->rows = calloc((G->size + 1), sizeof(int));
	
	pat_t *expanded_patt = malloc(sizeof(pat_t));
	expanded_patt->rows = calloc((G->size + 1), sizeof(int));
	
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
	G->cols = malloc(G->nnz * sizeof(int));
	G->rows = malloc((G->size + 1) * sizeof(int));
	
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
		
		int rowelems = expanded_patt->rows[i + 1] - expanded_patt->rows[i];
		int *arrayindex = calloc(rowelems, sizeof(int));
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
		
		int poscount = 0;
		
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
		
		dgesv( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
		
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
	int diagcounter = 0, zerocounter = 0;
	diag = calloc(G->size, sizeof(double));
	
	double *diagvalue = calloc(G->size, sizeof(double));
	
	for (i = 0; i < G->size; i++){
		diagvalue[i] = 1E-03*getvalue_mat(i, i, G);
	}
	
	double zero = 1E-50;
	
	for (i = 0; i < G->size; i++){
		for (j = G->rows[i]; j < G->rows[i + 1]; j++){
			if (G->cols[j] == i){diag[i] = sqrt(fabs(G->values[j]));}
		}
		if (fabs(diag[i]) < zero) {
			diagcounter++;
		}
	}
	
	
	for (i = 0; i < G->size; i++){
		zero = diagvalue[i];
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
	gnoz->cols = calloc(gnoz->nnz, sizeof(int));
	gnoz->rows = calloc((gnoz->size + 1), sizeof(int));

	
	for (i = 0; i < G->size; i++){
		zero = diagvalue[i];
		for (j = G->rows[i]; j < G->rows[i + 1]; j++){
			if (fabs(G->values[j]) > zero){
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

	G->nnz = gnoz->nnz;
	G->values = gnoz->values;
	G->cols = gnoz->cols;
	G->rows = gnoz->rows;
	
	free(diagvalue);
	free(diag); free(pattern); free(expanded_patt); 
	free(errorrow); free(residrow); free(errrow);
}
	
	
	
/* 		Function name: transpose
 * 		Purpose: Transposes a given matrix in CSR format
 * 
 * 		Inputs: Pointer to matrix (mat_t), pointer to transposed matrix (mat_t)
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- No nees to pass 2 pointers, return pointer
 *
 * 		Others: 
 */

void transpose(mat_t *G, mat_t *Gtransp){
	
	int i = 0, j = 0, k = 0;
	int *contc = calloc(G->size, sizeof(int));
	
	Gtransp->nnz = G->nnz;
	Gtransp->values = calloc(Gtransp->nnz, sizeof(double));
	Gtransp->cols = calloc(Gtransp->nnz, sizeof(int));
	Gtransp->rows = calloc((Gtransp->size + 1), sizeof(int));
	
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
	
	
	
/* 		Function name: transposeCSC
 * 		Purpose: Transposes a given matrix in CSC format
 * 
 * 		Inputs: Pointer to matrix (mat_t), pointer to transposed matrix (mat_t)
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- No nees to pass 2 pointers, return pointer
 *
 * 		Others: 
 */

void transposeCSC(mat_t *G, mat_t *Gtransp){
	
	int i = 0;
	
	Gtransp->nnz = G->nnz;
	Gtransp->values = calloc(Gtransp->nnz, sizeof(double));
	Gtransp->rows = calloc(Gtransp->nnz, sizeof(int));
	Gtransp->cols = calloc((Gtransp->size + 1), sizeof(int));
	
	for (i = 0; i < Gtransp->nnz; i++){
		Gtransp->values = G->values;
		Gtransp->rows = G->cols;
	}
	for (i = 0; i < Gtransp->size + 1; i++){
		Gtransp->cols = G->rows;
	}
	
}
	
	
	
/* 		Function name: CSCtoCOO
 * 		Purpose: Converts a CSC formatted matrix to COO block format to prepare it for Matrix Vector parallel multiplications. 
 * 
 * 		Inputs: Pointer to matrix (mat_t), pointer to COO matrix (mat_t), pointer to array of limits, threads executing code
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- Check limits to exactly distribute load
 *
 * 		Others: 
 */

void CSCtoCOO(mat_t *Gtransp, mat_t *GtranspCOO, int *limits, int nthreads){
	
	int i, j, k;
	int counter = 0;
	int sum = 0;
	int *contr = calloc(Gtransp->size, sizeof(int));
	int *rowlim = calloc(nthreads + 1, sizeof(int));
	int *counterrows = calloc(nthreads, sizeof(int));
	int elems = (int) ceil(((double) Gtransp->nnz) / ((double) nthreads));
	
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




/* 		Function name: filterG
 * 		Purpose: Limits G to cache lines of A
 * 
 * 		Inputs: Pointer to matrix G (mat_t), pointer to matrix A (mat_t), size, Pointer to X array from general problem
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- No need to pass 2 pointers, return pointer
 *
 * 		Others: 
 * 
 * 			- We saw that doing this yields to no convergence, since we remove information
 */

void filterG(mat_t *G, mat_t *A, int dim, double *r){
	
	int i = 0, j = 0, k = 0, l = 0;
	int cachecol = 0;
	long long cacheline = 0;
	int testcol = 0;
	int counter = 0;
	
	
	mat_t *Gfinal = malloc(sizeof(mat_t));
	Gfinal->size = dim;
	Gfinal->nnz = G->nnz;
	Gfinal->rows = calloc(dim + 1, sizeof(int));
	Gfinal->cols = calloc(Gfinal->nnz, sizeof(int));
	
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
	Gfinal->cols = realloc(Gfinal->cols, Gfinal->nnz*sizeof(int));
	
	printf("Gvalues: %i\nGfinalvalues: %i\n\n", G->nnz, counter);
	
	G->nnz = Gfinal->nnz;
	G->cols = Gfinal->cols;
	G->rows = Gfinal->rows;
	
}














/* FILTERING STUFF */



/* 		Function name: sumupdown
 * 		Purpose: Gets row positions of new filtered matrix combining upper and lower rows
 * 
 * 		Inputs: Pointer to matrix G (mat_t), pointer to array column positions (mat_t), row
 * 
 * 		Outputs:
 * 			- elements of new row
 * 
 * 		TODO:
 * 
 * 			- 
 *
 * 		Others: 
 * 
 * 			- Not used. Yields to a non-symmetric matrix
 */

int sumupdown(mat_t *A, int *col, int ract){
	
	int i = 0; 
	int dim = A->size;
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
	
	int *tmpcol = calloc(counter, sizeof(int));
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



/* 		Function name: sumLLt
 * 		Purpose: Gets row positions of new filtered matrix combining right side elements, lower row element and right side lower row element
 * 
 * 		Inputs: Pointer to matrix G (mat_t), pointer to array column positions (mat_t), row
 * 
 * 		Outputs:
 * 			- elements of new row
 * 
 * 		TODO:
 * 
 * 			- 
 *
 * 		Others: 
 * 
 * 			- Not used. Yields to a non-symmetric matrix
 */

int sumLLt(mat_t *A, int *col, int ract){
	
	int i = 0; 
	int dim = A->size;
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
	
	int *tmpcol = calloc(counter, sizeof(int));
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



/* 		Function name: smoothmatLLt
 * 		Purpose: Smooths matrix using sumLLt filter
 * 
 * 		Inputs: Pointer to matrix A (mat_t) to be filtered
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- 
 *
 * 		Others: 
 * 
 * 			- 
 */

void smoothmatLLt(mat_t *A){
	
	int dim = A->size;
	
	int maxAsizemult = (A->nnz/A->size + 1) * 10;
	double *storevals = calloc(maxAsizemult*A->nnz, sizeof(double));
	int *storecols = calloc(maxAsizemult*A->nnz, sizeof(int));
	int *storerows = calloc(dim + 1, sizeof(int));
	int pos = 0;
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
	A->cols = realloc(A->cols, pos*sizeof(int));
	
	
	for (i = 0; i < pos; i++){
		A->values[i] = storevals[i];
		A->cols[i] = (int) storecols[i];
	}
	
	for (i = 0; i < (dim + 1); i++) A->rows[i] = (int) storerows[i];
	
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




/* 		Function name: smoothTOmatLLt
 * 		Purpose: Smooths matrix using sumLLt filter and stores it to a filtered matrix
 * 
 * 		Inputs: Pointer to matrix A (mat_t) to be filtered, Pointer to matrix Af (mat_t)
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- 
 *
 * 		Others: 
 * 
 * 			- 
 */

void smoothTOmatLLt(mat_t *A, mat_t *Af){
	
	int dim = A->size;
	
	int maxAsizemult = (A->nnz/A->size + 1) * 10;
	double *storevals = calloc(maxAsizemult*A->nnz, sizeof(double));
	int *storecols = calloc(maxAsizemult*A->nnz, sizeof(int));
	int *storerows = calloc(dim + 1, sizeof(int));
	int pos = 0;
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
	Af->cols = calloc(pos, sizeof(int));
	Af->rows = calloc(dim + 1, sizeof(int));
	Af->size = dim;
	
	for (i = 0; i < pos; i++){
		Af->values[i] = storevals[i];
		Af->cols[i] = (int) storecols[i];
	}
	
	for (i = 0; i < (dim + 1); i++) Af->rows[i] = (int) storerows[i];
	
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




/* 		Function name: sumL
 * 		Purpose: Gets filtered positions using lower diagonal positions
 * 
 * 		Inputs: Pointer to matrix G (mat_t), pointer to array column positions (mat_t), row
 * 
 * 		Outputs:
 * 			- elements of new row
 * 
 * 		TODO:
 * 
 * 			- 
 *
 * 		Others: 
 * 
 * 			- Not used. Yields to a non-symmetric matrix
 */


int sumL(mat_t *A, int *col, int ract){
	
	int i = 0; 
	int dim = A->size;
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
	
	int *tmpcol = calloc(counter, sizeof(int));
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



/* 		Function name: smoothmatL
 * 		Purpose: Smooths matrix using sumL filter
 * 
 * 		Inputs: Pointer to matrix A (mat_t) to be filtered
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- 
 *
 * 		Others: 
 * 
 * 			- 
 */

void smoothmatL(mat_t *A){
	
	int dim = A->size;
	
	int maxAsizemult = (A->nnz/A->size + 1) * 10;
	double *storevals = calloc(maxAsizemult*A->nnz, sizeof(double));
	int *storecols = calloc(maxAsizemult*A->nnz, sizeof(int));
	int *storerows = calloc(dim + 1, sizeof(int));
	int pos = 0;
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
	A->cols = realloc(A->cols, pos*sizeof(int));
	
	
	for (i = 0; i < pos; i++){
		A->values[i] = storevals[i];
		A->cols[i] = (int) storecols[i];
	}
	
	for (i = 0; i < (dim + 1); i++) A->rows[i] = (int) storerows[i];
	
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




/* 		Function name: smoothTOmatL
 * 		Purpose: Smooths matrix using sumL filter and stores it to a filtered matrix
 * 
 * 		Inputs: Pointer to matrix A (mat_t) to be filtered, Pointer to matrix Af (mat_t)
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- 
 *
 * 		Others: 
 * 
 * 			- 
 */

void smoothTOmatL(mat_t *A, mat_t *Af){
	
	int dim = A->size;
	
	int maxAsizemult = (A->nnz/A->size + 1) * 10;
	double *storevals = calloc(maxAsizemult*A->nnz, sizeof(double));
	int *storecols = calloc(maxAsizemult*A->nnz, sizeof(int));
	int *storerows = calloc(dim + 1, sizeof(int));
	int pos = 0;
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
	Af->cols = calloc(pos, sizeof(int));
	Af->rows = calloc(dim + 1, sizeof(int));
	Af->size = dim;
	
	for (i = 0; i < pos; i++){
		Af->values[i] = storevals[i];
		Af->cols[i] = (int) storecols[i];
	}
	
	for (i = 0; i < (dim + 1); i++) Af->rows[i] = (int) storerows[i];
	
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




/* 		Function name: smootharrL
 * 		Purpose: Smooths array using sumL filter and stores it to a filtered matrix
 * 
 * 		Inputs: Pointer to array b to be filtered, array size
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- 
 *
 * 		Others: 
 * 
 * 			- 
 */


void smootharrL(double *b, int dim){
	
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



/* 		Function name: backsmootharrLt
 * 		Purpose: Reverse smoothing process for arrays
 * 
 * 		Inputs: Pointer to array b to be filtered, array size
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- 
 *
 * 		Others: 
 * 
 * 			- 
 */

void backsmootharrLt(double *b, int dim){
	
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




/* 		Function name: backsmootharrLtV2
 * 		Purpose: Reverse smoothing process for arrays
 * 
 * 		Inputs: Pointer to array b to be filtered, array size
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- 
 *
 * 		Others: 
 * 
 * 			- To use only in filtered preconditioner function ()
 */

void backsmootharrLtV2(double *b, int dim, int *arrayindex){
	
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




/* 		Function name: solveLU
 * 		Purpose: LU solver 
 * 
 * 		Inputs: Pointer to matrix (mat_t), size, array b (Ax=b)
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- 
 *
 * 		Others: 
 * 
 * 			- 
 */

void solveLU(mat_t *A, int dim, double *b){
	
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
	
	dgesv( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
	
}




/* 		Function name: impLU
 * 		Purpose: Implemented LU solver: Matrix conversion
 * 
 * 		Inputs: Pointer to matrix (mat_t), size, array b (Ax=b)
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- 
 *
 * 		Others: 
 * 
 * 			- 
 */

void impLU(mat_t *A, int dim, double *b){
	
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




/* 		Function name: precond_vf
 * 		Purpose: Preconditioner for filtered cases. Exactly the same as precond but with filtering options
 * 
 * 		Inputs: Pointer to matrix (mat_t), size, array b (Ax=b)
 * 
 * 		Outputs:
 * 			- void
 * 
 * 		TODO:
 * 
 * 			- Adapt to use with another variable (not rhs)
 *
 * 		Others: 
 * 
 * 			- (Not in use)
 */

void precond_vf(mat_t *mat, mat_t *G, double *xfinal){
	
	int i = 0, j = 0, k = 0;
	
	pat_t *pattern = malloc(sizeof(pat_t));
	pattern->rows = calloc((G->size + 1), sizeof(int));
	
	pat_t *expanded_patt = malloc(sizeof(pat_t));
	expanded_patt->rows = calloc((G->size + 1), sizeof(int));
	
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
	G->cols = malloc(G->nnz * sizeof(int));
	G->rows = malloc((G->size + 1) * sizeof(int));
	
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
		
		int rowelems = expanded_patt->rows[i + 1] - expanded_patt->rows[i];
		int *arrayindex = calloc(rowelems, sizeof(int));
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
		
		int poscount = 0;
		
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
		
		dgesv( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
		
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
	int diagcounter = 0, zerocounter = 0;
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
	gnoz->cols = calloc(gnoz->nnz, sizeof(int));
	gnoz->rows = calloc((gnoz->size + 1), sizeof(int));

	
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
	
	free(diag); free(pattern); free(expanded_patt); 
	free(errorrow); free(residrow); free(errrow);
}
