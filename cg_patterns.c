
#include "cg_patterns.h"

#define cachesize 8

extern int percentpattern;
extern int patternpower;
extern FILE *outmn;
extern char matname[256];
extern mat_t *powmat;



/* 		Function name: compd
 * 		Purpose: compare function for doubles
 * 
 * 		Inputs: Compare values
 * 
 * 		Outputs:
 * 			- -1 on a>b
 * 			- 1 on b>a
 * 			- 0 a=b
 * 		TODO:
 *
 * 		- 
 * 
 * 		Others: 
 * 
 * 		- 
 */

int compd(const void * a, const void * b){
  if (*(double*)a > *(double*)b) return 1;
  else if (*(double*)a < *(double*)b) return -1;
  else return 0;  
}




/* 		Function name: comp
 * 		Purpose: compare function for integers
 * 
 * 		Inputs: Compare values
 * 
 * 		Outputs:
 * 			- Integer difference 
 * 
 * 		TODO:
 *
 * 		- 
 * 
 * 		Others: 
 * 
 * 		- 
 */

int comp (const void *a, const void *b){
	int *x = (int *) a;
	int *y = (int *) b;
	return *x - *y;
}




/* 		Function name: pexpandpattern
 * 		Purpose: Parallel pattern expanding
 * 
 * 		Inputs: Initial Pattern
 * 
 * 		Outputs: Final Pattern
 * 
 * 		TODO:
 *
 * 		- Check limits between threads
 * 
 * 		Others: 
 * 
 * 		- 
 */

void pexpandpattern(int dim, pat_t *pattern, pat_t *expanded_patt, double *xfinal){
	
	int i = 0, j = 0, k = 0, l = 0;
	int *poolcols = malloc(cachesize*pattern->nnz*sizeof(int));
	int *poolrows = calloc((dim + 1), sizeof(int));
	int numrowelems = 0;
	int cachecol = 0;
	int cacheline = 0;
	int cacheadded = 0;
	int cacheelements = 0;
	int tmpcol = 0;
	int calc = 0;
	int maxcacheelements = 0;
	int initk = 0;
	int addedelements = 0;
	int addedrow = 0;
	int nthreads = omp_get_max_threads(); 
	int isincacheline = 0;
	
	double maxtmp = (((double) pattern->nnz) * ((double) percentpattern)) / 100.0;
	
	int maxaddedelements = floor(maxtmp / nthreads);
	
	
	
	#pragma omp parallel for
	for (i = 0; i < cachesize*pattern->nnz; i++) poolcols[i] = -1;
	
	
	
	#pragma omp parallel private(j, numrowelems, cachecol, cacheline, cacheadded, cacheelements, k, tmpcol, l, calc, maxcacheelements, initk, addedelements, addedrow, isincacheline)
	{
		calc = 0;
		addedelements = 0;
		
		#pragma omp for
		for (i = 0; i < dim; i++){
			
			addedrow = 0;
			numrowelems = pattern->rows[i];
			
			for(j = pattern->rows[i]; j < pattern->rows[i + 1]; j++){
				
				cachecol = pattern->cols[j];
				cacheline = (int) (((long long) &xfinal[pattern->cols[j]]%64)/8);
				cacheadded = 0;
				cacheelements = 0;
				
				
				for (k = cacheline; k < cachesize; k++){
					tmpcol = cachecol + k;
					
					for (l = j; l < pattern->rows[i + 1]; l++){
						if (tmpcol == pattern->cols[l]) cacheelements++;
						if ((pattern->cols[l] - cachecol) >= cachesize) break;
					}
				}
				
				calc += cacheelements * ((int) percentpattern); // si > 75 -> 1, si > 175 -> 2
				maxcacheelements = 0;
				
				while (calc > 66){
					calc -= 100;
					maxcacheelements++;
				}
				
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
				initk = cacheline - maxcacheelements;
				if (initk < 0) initk = 0;
				
				
				
				
				for (k = initk; k < cachesize; k++){
					
					isincacheline = 0;
					tmpcol = cachecol - cacheline + k;
					
					for (l = j; l < pattern->rows[i + 1]; l++){
						if ((int) pattern->cols[l] == tmpcol){
							isincacheline = 1;
							break;
						}
					}
					
					if ((tmpcol >= 0) && (tmpcol <= i) && (isincacheline == 0) && (addedelements < maxaddedelements) && (cacheadded < maxcacheelements)){
						
						
						poolcols[numrowelems*cachesize + (j - numrowelems)*cachesize + k] = tmpcol;
						poolrows[i + 1]++;
						addedelements++;
						cacheadded++;
						addedrow++;
						
					}
					else if ((tmpcol >= 0) && (tmpcol <= i) && (isincacheline == 1)){
						
						
						poolcols[numrowelems*cachesize + (j - numrowelems)*cachesize + k] = tmpcol;
						poolrows[i + 1]++;
						if (tmpcol != cachecol) j++;
						
					}
					else{}
				}
				
				
				calc += (maxcacheelements - cacheadded) * 100;
				
			}
			
		}
	}
	
	for (i = 1; i < dim + 1; i++){poolrows[i] += poolrows[i - 1];}
	#pragma omp parallel for
	for (i = 0; i < dim + 1; i++){expanded_patt->rows[i] = (int) poolrows[i];}
	
	
	expanded_patt->nnz = expanded_patt->rows[dim];
	expanded_patt->cols = calloc(expanded_patt->nnz, sizeof(int));
	
	
	
	int count = 0;
	for (i = 0; i < cachesize*pattern->nnz; i++){
		if (poolcols[i] >= 0){
			expanded_patt->cols[count] = (int) poolcols[i];
			count++;
		}
	}
	
	
	printf("\n\n%u -> ^%.2lf%%\n\n", expanded_patt->nnz, 100 * (((double) expanded_patt->nnz - (double) pattern->nnz)/ ((double) pattern->nnz)));
	
	
	free(poolcols); free(poolrows);
	
}
	
	
	
/* 		Function name: optexpandpattern
 * 		Purpose: Parallel pattern expanding to positions that do not generate zeroes
 * 
 * 		Inputs: Initial Pattern
 * 
 * 		Outputs: Final Pattern
 * 
 * 		TODO:
 *
 * 		- Check limits between threads
 * 
 * 		Others: 
 * 
 * 		- When using it percent pattern is the number of times the procedure of adding 1 element to a row is performed
 */

void optexpandpattern(int dim, pat_t *pattern, pat_t *expanded_patt, double *xfinal, mat_t *A){
	
	
	int i, j, k, l, m;
	int cachecol = 0, maxloop = 0, isincacheline = 0, numrowelems = 0;
	int testcol = 0;
	long long cacheline = 0;
	int *poolcols = malloc(cachesize*pattern->nnz*sizeof(int));
	int *tmprows = calloc(dim + 1, sizeof(int));
	int addedelements = 0, totaladdedelements = 0;
	int counter = 0;
	int added = 0; 
	int maxaddedrow = 1;
	int addedrow = 0;
	int count = 0;
	int newsize = 0;
	
	int *rowelems = calloc(dim, sizeof(int));
	
	expanded_patt->nnz = pattern->nnz;
	expanded_patt->cols = calloc(pattern->nnz, sizeof(int));
	
	for (i = 0; i < expanded_patt->nnz; i++){expanded_patt->cols[i] = pattern->cols[i];}
	for (i = 0; i < dim + 1; i++){expanded_patt->rows[i] = pattern->rows[i];}
	
	
	#pragma omp parallel for
	for (i = 0; i < cachesize*expanded_patt->nnz; i++) poolcols[i] = -100;
	
	
	
	while (maxloop < percentpattern){
		
		newsize = cachesize*expanded_patt->nnz;
		poolcols = realloc(poolcols, newsize*sizeof(int));
		#pragma omp parallel for
		for (i = 0; i < cachesize*expanded_patt->nnz; i++) poolcols[i] = -100;
		#pragma omp parallel for
		for (i = 0; i < dim + 1; i++){tmprows[i] = 0;}
		
		addedelements = 0;
		counter = 0;
		
		#pragma omp parallel for private(numrowelems, cachecol, cacheline, j, k, l, m, testcol, isincacheline, added, addedrow) reduction(+:addedelements, counter)
		for (i = 0; i < dim; i++){
			
			addedrow = 0;
			numrowelems = expanded_patt->rows[i];
			
			for(j = expanded_patt->rows[i]; j < expanded_patt->rows[i + 1]; j++){
				
				cachecol = (int) expanded_patt->cols[j];
				cacheline = (((long long) &xfinal[expanded_patt->cols[j]]%64)/8);
				
				for (k = 0; k < cachesize; k++){
					
					isincacheline = 0;
					testcol = cachecol - (int) cacheline + k;
					
					for (l = expanded_patt->rows[i]; l < expanded_patt->rows[i + 1]; l++){
						if ((int) expanded_patt->cols[l] == testcol){
							isincacheline = 1;
							break;
						}
					}
					
					if ((testcol >= 0) && (testcol <= i) && (isincacheline == 0) && (addedrow < maxaddedrow)){
						added = 0;
						for (l = pattern->rows[testcol]; l < pattern->rows[testcol + 1]; l++){
							for(m = expanded_patt->rows[i]; m < expanded_patt->rows[i + 1]; m++){
								if ((pattern->cols[l] == expanded_patt->cols[m]) && (added == 0)){
									poolcols[numrowelems*cachesize + (j - numrowelems)*cachesize + k] = testcol;
									addedelements++;
									tmprows[i + 1]++;
									added = 1;
									addedrow++;
								}
							}
						}
					}
					else if ((testcol >= 0) && (testcol <= i) && (isincacheline == 1)){
						
						poolcols[numrowelems*cachesize + (j - numrowelems)*cachesize + k] = testcol;
						if (testcol != cachecol) j++;
						tmprows[i + 1]++;
						counter++;
					}
					else{/*poolcols[numrowelems*cachesize + (j - numrowelems)*cachesize + k] = -100;*/}
				}
			}
			
			rowelems[i] += addedrow;
		}
		
		totaladdedelements += addedelements;
		expanded_patt->nnz = counter + addedelements;
		expanded_patt->cols = realloc(expanded_patt->cols, expanded_patt->nnz*sizeof(int));
		
		count = 0;
		for (i = 0; i < newsize; i++){
			if (poolcols[i] >= 0){
				expanded_patt->cols[count] = (int) poolcols[i];
				count++;
			}
		}
		
		for (i = 1; i < dim + 1; i++){tmprows[i] += tmprows[i - 1];}
		
		
		
		#pragma omp parallel for
		for (i = 0; i < dim + 1; i++){expanded_patt->rows[i] = (int) tmprows[i];}
		
		
		char buf_p[256];
		snprintf(buf_p, sizeof buf_p, "../Outputs/cg/DataCG/addedelems/add_%s_%i_%i_%i.mtx", matname, percentpattern, patternpower, maxloop);
		FILE *tests= fopen(buf_p, "w");
		int sum = 0, pos = 0;
		double avg = 0.0;
		for(i = 0; i < dim; i++){
			sum += rowelems[i];
			if ((i%1000 == 0) && (i > 0)){
				avg = (double) sum / 1000.0;
				fprintf(tests, "%i %lf\n", pos, avg);
				pos++;
				sum = 0;
			}
		}
		fclose(tests);
		
		maxloop++;
	}
	
	
	int top = 0;
	
	for (i = 0; i < dim; i++) top+=rowelems[i];
	
	
	
	top = 0;
	for (i = 0; i < dim; i++){
		if (rowelems[i] > top) top = rowelems[i];
	}
	
	int *histogram = calloc(top + 1, sizeof(int));
	
	for (i = 0; i < dim; i++){
		histogram[rowelems[i]]++;
	}
	
	char buf_h[256];
	snprintf(buf_h, sizeof buf_h, "../Outputs/cg/DataCG/hist_%s_%i_%i", matname, patternpower, percentpattern);
	FILE *histograms= fopen(buf_h, "w");
	for(i = 0; i < top + 1; i++){
		fprintf(histograms, "%i %i\n", i, histogram[i]);
	}
	fclose(histograms);
	
	
	
	printf("INITIAL ELEMENTS:\t%i\n", pattern->nnz);
	printf("ADDED ELEMENTS:\t\t%i\n", addedelements);
	printf("SUM ADDED ELEMENTS:\t%i\n", totaladdedelements);
	printf("TOTAL ELEMENTS:\t\t%i\n", expanded_patt->nnz);
	printf("INCREASE ELEMENTS:\t%.2lf%%\n", ((double) expanded_patt->nnz)/((double) pattern->nnz)*100.0 - 100.0);
	
	free(histogram);
	free(poolcols); free(tmprows);
}



/* 		Function name: expandpattern
 * 		Purpose: Pattern expanding
 * 
 * 		Inputs: Initial Pattern
 * 
 * 		Outputs: Final Pattern
 * 
 * 		TODO:
 *
 * 
 * 		Others: 
 * 
 * 		- Deprecated
 */

void expandpattern(int dim, pat_t *pattern, pat_t *expanded_patt, double *xfinal){
	
	int i = 0, j = 0, k = 0, l = 0;
	int *elementstoadd = calloc(8, sizeof(int));
	int cacheelements = 0;
	int maxcacheelements = 0;
	int oldcachecol;
	int cacheadded = 0, oldcacheline = 0;
	int calc = 0;
	int tmpcol = 0;
	int cachecol = 0;
	int addedelements = 0, maxaddedelements = 0;
	int counterA = 0, counterB = 0;
	long long cacheline = 0;
	int initk = 0;
	
	double maxtmp = (((double) pattern->nnz) * ((double) percentpattern)) / 100.0;
	
	maxaddedelements = floor(maxtmp);
	expanded_patt->nnz = pattern->nnz + maxaddedelements;
	expanded_patt->cols = calloc(expanded_patt->nnz, sizeof(int));
	
	printf("%u -> %u\n", pattern->nnz, expanded_patt->nnz);
	
	for(i = 0; i < dim; i++){
		
		oldcachecol = -10;
		oldcacheline = 0;
		
		
		for(j = pattern->rows[i]; j < pattern->rows[i + 1]; j++){
			
			cachecol = pattern->cols[j];
			cacheline = (((long long) &xfinal[pattern->cols[j]]%64)/8);
			
			if (((cachecol - oldcachecol) > 7) || (((cachecol - oldcachecol) <= 7) && ((cacheline - oldcacheline) < 0))){
				
				oldcachecol = cachecol;
				oldcacheline = cacheline;
				cacheadded = 0;
				cacheelements = 1;
				
				for(k = 0; k < 8; k++){elementstoadd[k] = 0;}
				
				for (k = 1; k < (8 - (int) cacheline); k++){ // Right-most elements
					
					tmpcol = pattern->cols[j] + k;
					elementstoadd[(int) cacheline + k] = 1;
					
					for (l = pattern->rows[i]; l < pattern->rows[i + 1]; l++){
						
						if ((pattern->cols[l] == tmpcol) || (tmpcol > i)) {
							elementstoadd[(int) cacheline + k] = 0;
						}
						if (pattern->cols[l] == tmpcol) cacheelements++;
					}
				}
				for (k = (int) cacheline; k > 0 ; k--){		// Left-most elements
					
					tmpcol = pattern->cols[j] - k;
					elementstoadd[(int) cacheline - k] = 1;
					
					for (l = pattern->rows[i]; l < pattern->rows[i + 1]; l++){
						if ((pattern->cols[l] == tmpcol) || (tmpcol < 0)){
							elementstoadd[cacheline - k] = 0;
						}
						if (pattern->cols[l] == tmpcol) cacheelements++;
					}
				}
				
				
				
				calc += cacheelements * ((int) percentpattern); // si > 75 -> 1, si > 175 -> 2
				maxcacheelements = 0;
				
				while (calc > 66){
					calc -= 100;
					maxcacheelements++;
				}
				
				
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
				initk = (int) cacheline - (int) maxcacheelements;
				if (initk < 0) initk = 0;
				
				for (k = initk; k < 8; k++){
					if ((elementstoadd[k] == 1) && (addedelements < maxaddedelements) && (cacheadded < maxcacheelements) && ((k + pattern->cols[j] - (int) cacheline) < i) && ((k + pattern->cols[j] - (int) cacheline) >= 0)){
						expanded_patt->cols[counterB] = pattern->cols[j] + k - (int) cacheline;
						
						
						counterB++;
						addedelements++;
						cacheadded++;
					}
					else if ((elementstoadd[k] == 0) && ((k + pattern->cols[j] - (int) cacheline) <= i) && ((k + pattern->cols[j] - (int) cacheline) >= 0) && (pattern->cols[counterA] == (pattern->cols[j] + k - (int) cacheline))){
						expanded_patt->cols[counterB] = pattern->cols[counterA];
						
						
						counterA++;
						counterB++;
						
					}
				}
				calc += (maxcacheelements - cacheadded) * 100;
			}
		}
		expanded_patt->rows[i + 1] = pattern->rows[i + 1] + addedelements;
	}
	expanded_patt->nnz = expanded_patt->rows[dim];
	
	
	printf("\n\n%u -> ^%.2lf%%\n\n", expanded_patt->nnz, 100 * (((double) expanded_patt->nnz - (double) pattern->nnz)/ ((double) pattern->nnz)));
	
	
	free(elementstoadd);
}




/* 		Function name: ltp
 * 		Purpose: Gets lower triangular pattern of input matrix
 * 
 * 		Inputs: Initial Pattern
 * 
 * 		Outputs: Final Pattern
 * 
 * 		TODO:
 *
 * 		- 
 * 
 * 		Others: 
 * 
 * 		- 
 */

void ltp(mat_t *mat, pat_t *pattern, pat_t *expanded_patt){

	int i = 0, j = 0;
	int patterncount = 0;
	int counter = 0;

	for (i = 0; i < mat->size; i++){							// For every row in the matrix
		for (j = mat->rows[i]; j < mat->rows[i + 1]; j++){		// Gets only lower triangle values
			if (mat->cols[j] <= i){patterncount++;}
		}
		pattern->rows[i + 1] = patterncount;
	}
	pattern->nnz = patterncount;
	pattern->cols = malloc(pattern->nnz*sizeof(int));
	for (i = 0; i < mat->size; i++){							// For every row in the matrix
		for (j = mat->rows[i]; j < mat->rows[i + 1]; j++){		// Gets only lower triangle values
			if (mat->cols[j] <= i){
				pattern->cols[counter] = mat->cols[j];
				counter++;
			}
		}
	}

	expanded_patt->rows = pattern->rows;
	expanded_patt->cols = pattern->cols;
	expanded_patt->nnz = pattern->nnz;
}



/* 		Function name: flt
 * 		Purpose: Gets full lower triangular pattern
 * 
 * 		Inputs: Initial Pattern
 * 
 * 		Outputs: Final Pattern
 * 
 * 		TODO:
 *
 * 
 * 		Others: 
 * 
 * 		- 
 */

void flt(mat_t *mat, pat_t *pattern, pat_t *expanded_patt){

	int i = 0, j = 0;
	int counter = 0, counter2 = 0;

	for (i = 0; i < mat->size; i++){							// For every row in the matrix
		pattern->rows[i + 1] = pattern->rows[i] + i + 1;
	}
	pattern->nnz = pattern->rows[mat->size];
	pattern->cols = malloc(pattern->nnz*sizeof(int));
	for (i = 0; i < mat->size; i++){							// For every row in the matrix
		counter = 0;
		for (j = pattern->rows[i]; j < pattern->rows[i + 1]; j++){		// Gets only lower triangle values
			pattern->cols[counter2] = counter;
			counter++;
			counter2++;
		}
	}

	expanded_patt->rows = pattern->rows;
	expanded_patt->cols = pattern->cols;
	expanded_patt->nnz = pattern->nnz;

}



/* 		Function name: perclt
 * 		Purpose: Gets lower triangular pattern and expands it
 * 
 * 		Inputs: Initial Pattern
 * 
 * 		Outputs: Final Pattern
 * 
 * 		TODO:
 *
 * 		- 
 * 
 * 		Others: 
 * 
 * 		- 
 */

void perclt(mat_t *mat, pat_t *pattern, pat_t *expanded_patt, double *xfinal){

	int i = 0, j = 0;
	int counter = 0;
	int patterncount = 0;


	// WE GET INITIAL PATTERN PROPERLY STORED

	for (i = 0; i < mat->size; i++){	// First get initial pattern
		for (j = mat->rows[i]; j < mat->rows[i + 1]; j++){
			if (mat->cols[j] <= i){patterncount++;}
		}
		pattern->rows[i + 1] = patterncount;
	}

	pattern->nnz = patterncount;
	pattern->cols = malloc(pattern->nnz*sizeof(int));

	for (i = 0; i < mat->size; i++){
		for (j = mat->rows[i]; j < mat->rows[i + 1]; j++){
			if (mat->cols[j] <= i){
				pattern->cols[counter] = mat->cols[j];
				counter++;
			}
		}
	}
	
	pexpandpattern(mat->size, pattern, expanded_patt, xfinal);
	
}



/* 		Function name: powerA
 * 		Purpose: Powers a pattern to the desired power
 * 
 * 		Inputs: Initial Pattern
 * 
 * 		Outputs: Final Pattern
 * 
 * 		TODO:
 *
 * 		- Parallelize
 * 
 * 		Others: 
 * 
 * 		- 
 */

void powerA(mat_t *mat, pat_t *pattern, pat_t *expanded_patt, double *xfinal){

	int i = 0, j = 0, k = 0, l = 0;
	int counter = 0;
	int patterncount = 0;
	int elements = 0;
	int dim = mat->size;

	
	pat_t *apower;
	apower = malloc(sizeof(mat_t));
	apower->rows = calloc((dim + 1), sizeof(int));
	apower->cols = calloc((20 * patternpower * mat->nnz), sizeof(int));
	
	pat_t *apowertemp;
	apowertemp = malloc(sizeof(mat_t));
	apowertemp->rows = calloc((dim + 1), sizeof(int));
	apowertemp->cols = calloc((20 * patternpower * mat->nnz), sizeof(int));
	
	int *pospower = calloc(dim, sizeof(int));
	int *poscol = calloc(dim, sizeof(int));
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	apower->nnz = mat->nnz;
	
	for (i = 0; i < apower->nnz; i++) {
		apower->cols[i] = mat->cols[i];
	}
	
	for (i = 0; i < (dim + 1); i++) {
		apower->rows[i] = mat->rows[i];
	}
	
	/* POWER A to the patternpower */
	
	double start_usec = 0.0, end_usec = 0.0;
	double tsort = 0.0, tfor = 0.0, tend = 0.0;
	double ttotal = 0.0;
	
	ttotal = omp_get_wtime();
	for (int power = 1; power < patternpower; power++){
		for (j = 0; j < dim; j++){
			start_usec = omp_get_wtime();
			for (k = apower->rows[j]; k < apower->rows[j + 1]; k++){
				for (l = apower->rows[apower->cols[k]]; l < apower->rows[apower->cols[k] + 1]; l++){
					if (pospower[apower->cols[l]] == 1){
						continue;
					}
					else {
						pospower[apower->cols[l]] = 1;
						poscol[elements] = apower->cols[l];
						elements++;
					}
				}
			}
			end_usec = omp_get_wtime();
			tfor += end_usec - start_usec; 
			
			start_usec = omp_get_wtime();
			qsort(poscol, elements, sizeof(int), comp);
			end_usec = omp_get_wtime();
			tsort += end_usec - start_usec; 
			
			start_usec = omp_get_wtime();
			// Addelements to array
			for (i = 0; i < elements; i++){
				apowertemp->cols[counter] = poscol[i];
				counter++;
			}
			apowertemp->rows[j + 1] = counter;
			
			for (i = 0; i < elements; i++){
				pospower[poscol[i]] = 0;
				poscol[i] = 0;
			}
			
			end_usec = omp_get_wtime();
			tend += end_usec - start_usec; 
			elements = 0;
		}
		
		apowertemp->nnz = counter;
		counter = 0;
		for (j = 0; j < apowertemp->nnz; j++) {
			apower->cols[j] = apowertemp->cols[j];
		}
		for (j = 0; j < (dim + 1); j++) {
			apower->rows[j] = apowertemp->rows[j];
		}
		apower->nnz = apowertemp->nnz;
	}
	
	ttotal = omp_get_wtime() - ttotal;
	
// 	printf("For Time: %lf\tSorting Time: %lf\tEnd Time: %lf\tTotal Time: %lf", tfor, tsort, tend, ttotal);
	
	elements = 0;
	counter = 0;
	
	/* GET PATTERN OF POWERED A */
	
	for (i = 0; i < dim; i++){	// First get initial pattern
		for (j = apower->rows[i]; j < apower->rows[i + 1]; j++){
			if (apower->cols[j] <= i){patterncount++;}
		}
		pattern->rows[i + 1] = patterncount;
	}

	pattern->nnz = patterncount;
	pattern->cols = malloc(pattern->nnz*sizeof(int));

	for (i = 0; i < dim; i++){
		for (j = apower->rows[i]; j < apower->rows[i + 1]; j++){
			if (apower->cols[j] <= i){
				pattern->cols[counter] = apower->cols[j];
				counter++;
			}
		}
	}
	
	free(apower->cols);
	free(apower->rows);
	free(apower);
	free(apowertemp->cols);
	free(apowertemp->rows);
	free(apowertemp);
	free(pospower);
	free(poscol);
	
	
	
	
// 	optexpandpattern(mat->size, pattern, expanded_patt, xfinal, mat);
	
	
	double t1 = omp_get_wtime();
	pexpandpattern(mat->size, pattern, expanded_patt, xfinal);
	double t2 = omp_get_wtime();
	
	printf("TIME EXPAND PATTERN: %lf\n", t2 - t1);
	
	
// 	limitexpandpattern(mat->size, pattern, expanded_patt, xfinal, mat);
	
}




/* 		Function name: powerAf
 * 		Purpose: Powers a pattern to the desired power but for filtered matrices
 * 
 * 		Inputs: Initial Pattern
 * 
 * 		Outputs: Final Pattern
 * 
 * 		TODO:
 *
 * 		- Parallelize
 * 
 * 		Others: 
 * 
 * 		- 
 */

void powerAf(mat_t *mat, pat_t *pattern, pat_t *expanded_patt, double *xfinal){

	int i = 0, j = 0, k = 0, l = 0;
	int counter = 0;
	int patterncount = 0;
	int elements = 0;
	int dim = mat->size;

	
	pat_t *apower;
	apower = malloc(sizeof(mat_t));
	apower->rows = calloc((dim + 1), sizeof(int));
	apower->cols = calloc((20 * patternpower * mat->nnz), sizeof(int));
	
	pat_t *apowertemp;
	apowertemp = malloc(sizeof(mat_t));
	apowertemp->rows = calloc((dim + 1), sizeof(int));
	apowertemp->cols = calloc((20 * patternpower * mat->nnz), sizeof(int));
	
	int *pospower = calloc(dim, sizeof(int));
	int *poscol = calloc(dim, sizeof(int));
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	apower->nnz = mat->nnz;
	
	for (i = 0; i < apower->nnz; i++) {
		apower->cols[i] = mat->cols[i];
	}
	
	for (i = 0; i < (dim + 1); i++) {
		apower->rows[i] = mat->rows[i];
	}
	
	/* POWER A to the patternpower */
	
	double start_usec = 0.0, end_usec = 0.0;
	double tsort = 0.0, tfor = 0.0, tend = 0.0;
	double ttotal = 0.0;
	
	ttotal = omp_get_wtime();
	for (int power = 1; power < patternpower; power++){
		for (j = 0; j < dim; j++){
			start_usec = omp_get_wtime();
			for (k = apower->rows[j]; k < apower->rows[j + 1]; k++){
				for (l = apower->rows[apower->cols[k]]; l < apower->rows[apower->cols[k] + 1]; l++){
					if (pospower[apower->cols[l]] == 1){
						continue;
					}
					else {
						pospower[apower->cols[l]] = 1;
						poscol[elements] = apower->cols[l];
						elements++;
					}
				}
			}
			end_usec = omp_get_wtime();
			tfor += end_usec - start_usec; 
			
			start_usec = omp_get_wtime();
			qsort(poscol, elements, sizeof(int), comp);
			end_usec = omp_get_wtime();
			tsort += end_usec - start_usec; 
			
			start_usec = omp_get_wtime();
			// Addelements to array
			for (i = 0; i < elements; i++){
				apowertemp->cols[counter] = poscol[i];
				counter++;
			}
			apowertemp->rows[j + 1] = counter;
			
			for (i = 0; i < elements; i++){
				pospower[poscol[i]] = 0;
				poscol[i] = 0;
			}
			
			end_usec = omp_get_wtime();
			tend += end_usec - start_usec; 
			elements = 0;
		}
		
		apowertemp->nnz = counter;
		counter = 0;
		for (j = 0; j < apowertemp->nnz; j++) {
			apower->cols[j] = apowertemp->cols[j];
		}
		for (j = 0; j < (dim + 1); j++) {
			apower->rows[j] = apowertemp->rows[j];
		}
		apower->nnz = apowertemp->nnz;
	}
	
	ttotal = omp_get_wtime() - ttotal;
	
// 	printf("For Time: %lf\tSorting Time: %lf\tEnd Time: %lf\tTotal Time: %lf\n\n", tfor, tsort, tend, ttotal);
	
	elements = 0;
	counter = 0;
	
	/* GET PATTERN OF POWERED A */
	
	for (i = 0; i < dim; i++){	// First get initial pattern
		for (j = apower->rows[i]; j < apower->rows[i + 1]; j++){
			if (apower->cols[j] <= i){patterncount++;}
		}
		pattern->rows[i + 1] = patterncount;
	}

	pattern->nnz = patterncount;
	pattern->cols = malloc(pattern->nnz*sizeof(int));

	for (i = 0; i < dim; i++){
		for (j = apower->rows[i]; j < apower->rows[i + 1]; j++){
			if (apower->cols[j] <= i){
				pattern->cols[counter] = apower->cols[j];
				counter++;
			}
		}
	}
	
	free(apower->cols);
	free(apower->rows);
	free(apower);
	free(apowertemp->cols);
	free(apowertemp->rows);
	free(apowertemp);
	free(pospower);
	free(poscol);
	
	
	
	
// 	optexpandpattern(mat->size, pattern, expanded_patt, xfinal, mat);
	
	pexpandpattern(mat->size, pattern, expanded_patt, xfinal);
	
// 	limitexpandpattern(mat->size, pattern, expanded_patt, xfinal, mat);
	
}




/* 		Function name: OptpowerA
 * 		Purpose: Powers a pattern to the desired power but limitin it to lower triangular pattern. Ensures no zeroes are generated
 * 
 * 		Inputs: Initial Pattern
 * 
 * 		Outputs: Final Pattern
 * 
 * 		TODO:
 *
 * 		- Parallelize
 * 
 * 		Others: 
 * 
 * 		- 
 */

void OptpowerA(mat_t *mat, pat_t *pattern, pat_t *expanded_patt, double *xfinal){
	
	
	int i, j, k, l;
	int counter = 0;
	int elements = 0;
	int dim = mat->size;
	int *pospower = calloc(dim, sizeof(int));
	int *poscol = calloc(dim, sizeof(int));
	
	pat_t *apower;
	apower = malloc(sizeof(mat_t));
	apower->nnz = (mat->nnz - dim)/2 + dim;
	apower->rows = calloc((dim + 1), sizeof(int));
	apower->cols = calloc(apower->nnz, sizeof(int));
	
	pat_t *apowertmp;
	apowertmp = malloc(sizeof(mat_t));
	apowertmp->nnz = (mat->nnz - dim)/2 + dim;
	apowertmp->rows = calloc((dim + 1), sizeof(int));
	apowertmp->cols = calloc(10*patternpower*apowertmp->nnz, sizeof(int));
	
	
	// TEMPORARY
	pattern->nnz = apower->nnz;
	pattern->cols = calloc(pattern->nnz, sizeof(int));
	// TEMPORARY
	
	for (i = 0; i < dim; i++){
		for (j = mat->rows[i]; j < mat->rows[i + 1]; j++){
			if (mat->cols[j] <= i){
				pattern->cols[counter] = mat->cols[j];
				apower->cols[counter] = mat->cols[j];
				counter++;
			}
		}
		pattern->rows[i + 1] = counter;
		apower->rows[i + 1] = counter;
	}
	
	
	
	
	
	for (int power = 1; power < patternpower; power++){
		
		counter = 0;
		
		for (j = 0; j < dim; j++){
			for (k = apower->rows[j]; k < apower->rows[j + 1]; k++){
				for (l = apower->rows[apower->cols[k]]; l < apower->rows[apower->cols[k] + 1]; l++){
					if (pospower[apower->cols[l]] == 1){
						continue;
					}
					else {
						pospower[apower->cols[l]] = 1;
						poscol[elements] = apower->cols[l];
						elements++;
					}
				}
			}
			
			qsort(poscol, elements, sizeof(int), comp);
			
			for (i = 0; i < elements; i++){
				apowertmp->cols[counter] = poscol[i];
				counter++;
			}
			apowertmp->rows[j + 1] = counter;
			
			for (i = 0; i < elements; i++){
				pospower[poscol[i]] = 0;
				poscol[i] = 0;
			}
			
			elements = 0;
		}
		
		apowertmp->nnz = counter;
		apower->nnz = apowertmp->nnz;
		apower->cols = realloc(apower->cols, apower->nnz*sizeof(int));
		
		for (j = 0; j < apowertmp->nnz; j++) {
			apower->cols[j] = apowertmp->cols[j];
		}
		for (j = 0; j < (dim + 1); j++) {
			apower->rows[j] = apowertmp->rows[j];
		}
	}
	
	
	
	free(pattern->cols);
	free(pattern->rows);
	
	
	pattern->nnz = apower->nnz;
	pattern->rows = apower->rows;
	pattern->cols = apower->cols;
	
	
	
// 	free(apower);
	free(apowertmp);
	free(poscol);
	free(pospower);
	
	
	// AT THIS POINT I MUST HAVE A GOOD pattern
	
	optexpandpattern(mat->size, pattern, expanded_patt, xfinal, mat);
	
	
}




/* 		Function name: poweredA
 * 		Purpose: Powers a pattern to the desired power and calculates values
 * 
 * 		Inputs: Initial Pattern
 * 
 * 		Outputs: Final Pattern
 * 
 * 		TODO:
 *
 * 		- Parallelize
 * 
 * 		Others: 
 * 
 * 		- 
 */

void poweredA(mat_t *mat, int dim){

	int i = 0, j = 0, k = 0, l = 0;
	int counter = 0;
	int elements = 0;
	
	int cont = 0;
	
	mat_t *apower;
	apower = malloc(sizeof(mat_t));
	apower->rows = calloc((dim + 1), sizeof(int));
	apower->cols = calloc((10 * mat->nnz), sizeof(int));
	apower->values = calloc((10 * mat->nnz), sizeof(int));
	
	mat_t *apowertemp;
	apowertemp = malloc(sizeof(mat_t));
	apowertemp->rows = calloc((dim + 1), sizeof(int));
	apowertemp->cols = calloc((10 * mat->nnz), sizeof(int));
	apowertemp->values = calloc((10 * mat->nnz), sizeof(int));
	
	int *pospower = calloc(dim, sizeof(int));
	int *poscol = calloc(dim, sizeof(int));
	
	double value1 = 1.0, value2 = 1.0;
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	apower->nnz = mat->nnz;
	
	for (i = 0; i < apower->nnz; i++) {
		apower->cols[i] = mat->cols[i];
		apower->values[i] = mat->values[i];
	}
	
	for (i = 0; i < (dim + 1); i++) {
		apower->rows[i] = mat->rows[i];
	}
	
	
	
	/* POWER A to the patternpower */
	
	
	for (int power = 0; power < 1; power++){
		
		for (j = 0; j < dim; j++){
			for (k = apower->rows[j]; k < apower->rows[j + 1]; k++){
				for (l = apower->rows[apower->cols[k]]; l < apower->rows[apower->cols[k] + 1]; l++){
					if (pospower[apower->cols[l]] == 1){
						continue;
					}
					else {
						pospower[apower->cols[l]] = 1;
						poscol[elements] = apower->cols[l];
						elements++;
					}
				}
			}
			
			qsort(poscol, elements, sizeof(int), comp);
			
			// Addelements to array
			for (i = 0; i < elements; i++){
				apowertemp->cols[counter] = poscol[i];
				counter++;
			}
			apowertemp->rows[j + 1] = counter;
			
			
			for (i = 0; i < elements; i++){
				pospower[poscol[i]] = 0;
				poscol[i] = 0;
			}
			
			elements = 0;
		}
		
		apowertemp->nnz = counter;
		
		
		
// 		#pragma omp parallel for private(k, l, poscol, cont, value1, value2, pospower)
		for (j = 0; j < dim; j++){
			for (k = apowertemp->rows[j]; k < apowertemp->rows[j + 1]; k++){
				for (l = apowertemp->rows[j]; l < apowertemp->rows[j + 1]; l++){
					pospower[apowertemp->cols[l]] = 1;
					poscol[cont] = apowertemp->cols[l];
					cont++;
				}
				for (l = apowertemp->rows[apowertemp->cols[k]]; l < apowertemp->rows[apowertemp->cols[k] + 1]; l++){
					if (pospower[apowertemp->cols[l]] == 1){
						continue;
					}
					else{
						pospower[apowertemp->cols[l]] = 1;
						poscol[cont] = apowertemp->cols[l];
						cont++;
					}
				}
				
				qsort(poscol, cont, sizeof(int), comp);
				
				for (l = 0; l < cont; l++){
					
					value1 = getvalue_mat(j, poscol[l], apower);
					value2 = getvalue_mat(apowertemp->cols[k], poscol[l], mat);
					
					apowertemp->values[k] += (value1*value2);
				}
				
				for (l = 0; l < cont; l++){
					pospower[poscol[l]] = 0;
					poscol[l] = 0;
				}
				cont = 0;
			}
		}
		
		counter = 0;
		for (j = 0; j < apowertemp->nnz; j++) {
			apower->cols[j] = apowertemp->cols[j];
			apower->values[j] = apowertemp->values[j];
		}
		for (j = 0; j < (dim + 1); j++) {
			apower->rows[j] = apowertemp->rows[j];
		}
		apower->nnz = apowertemp->nnz;
		
	}
	
	
	powmat = malloc(sizeof(mat_t));
	powmat->nnz = apowertemp->nnz;
	powmat->cols = malloc(powmat->nnz*sizeof(int));
	powmat->values = malloc(powmat->nnz*sizeof(double));
	powmat->size = dim;
	powmat->rows = malloc((dim + 1) * sizeof(int));
	
	
	for (j = 0; j < apowertemp->nnz; j++) {
		powmat->cols[j] = apowertemp->cols[j];
		powmat->values[j] = apowertemp->values[j];
	}
	for (j = 0; j < (dim + 1); j++) {
		powmat->rows[j] = apowertemp->rows[j];
	}
	
	
	free(apower->cols);
	free(apower->rows);
	free(apower);
	free(apowertemp->cols);
	free(apowertemp->rows);
	free(apowertemp);
	free(pospower);
	free(poscol);
	
}

