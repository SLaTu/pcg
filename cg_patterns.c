
#include "cg_patterns.h"

extern unsigned int percentpattern;
extern unsigned int patternpower;
extern FILE *outmn;
extern char matname[256];

int compd(const void * a, const void * b){
  if (*(double*)a > *(double*)b) return 1;
  else if (*(double*)a < *(double*)b) return -1;
  else return 0;  
}

int comp (const void *a, const void *b){
	int *x = (int *) a;
	int *y = (int *) b;
	return *x - *y;
}


void optexpandpattern(unsigned int dim, pat_t *pattern, pat_t *expanded_patt, double *xfinal, mat_t *A){
	
	int cachesize = 8;
	
	int i, j, k, l, m;
	int cachecol = 0, maxloop = 0, testcol = 0, isincacheline = 0, numrowelems = 0;
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
	expanded_patt->cols = calloc(pattern->nnz, sizeof(unsigned int));
	
	for (i = 0; i < expanded_patt->nnz; i++){expanded_patt->cols[i] = pattern->cols[i];}
	for (i = 0; i < dim + 1; i++){expanded_patt->rows[i] = pattern->rows[i];}
	
	
	#pragma omp parallel for
	for (i = 0; i < cachesize*expanded_patt->nnz; i++) poolcols[i] = -100;
	
	
	
	while (maxloop < percentpattern){
		
		newsize = cachesize*expanded_patt->nnz;
		poolcols = realloc(poolcols, newsize*sizeof(unsigned int));
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
		expanded_patt->cols = realloc(expanded_patt->cols, expanded_patt->nnz*sizeof(unsigned int));
		
		count = 0;
		for (i = 0; i < newsize; i++){
			if (poolcols[i] >= 0){
				expanded_patt->cols[count] = (unsigned int) poolcols[i];
				count++;
			}
		}
		
		for (i = 1; i < dim + 1; i++){tmprows[i] += tmprows[i - 1];}
		
		
		
		#pragma omp parallel for
		for (i = 0; i < dim + 1; i++){expanded_patt->rows[i] = (unsigned int) tmprows[i];}
		
		
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

void expandpattern(unsigned int dim, pat_t *pattern, pat_t *expanded_patt, double *xfinal){
	
	unsigned int i = 0, j = 0, k = 0, l = 0;
	unsigned int *elementstoadd = calloc(8, sizeof(unsigned int));
	int cacheelements = 0;
	unsigned int maxcacheelements = 0;
	int oldcachecol;
	unsigned int cacheadded = 0, oldcacheline = 0;
	int calc = 0;
	int tmpcol = 0;
	unsigned int cachecol = 0;
	int addedelements = 0, maxaddedelements = 0;
	unsigned int counterA = 0, counterB = 0;
	long long cacheline = 0;
	int initk = 0;
	
	double maxtmp = (((double) pattern->nnz) * ((double) percentpattern)) / 100.0;
	
	maxaddedelements = floor(maxtmp);
	expanded_patt->nnz = pattern->nnz + maxaddedelements;
	expanded_patt->cols = calloc(expanded_patt->nnz, sizeof(unsigned int));
	
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
				
				for (k = 1; k < (8 - (unsigned int) cacheline); k++){ // Right-most elements
					
					tmpcol = pattern->cols[j] + k;
					elementstoadd[(unsigned int) cacheline + k] = 1;
					
					for (l = pattern->rows[i]; l < pattern->rows[i + 1]; l++){
						
						if ((pattern->cols[l] == tmpcol) || (tmpcol > i)) {
							elementstoadd[(unsigned int) cacheline + k] = 0;
						}
						if (pattern->cols[l] == tmpcol) cacheelements++;
					}
				}
				for (k = (unsigned int) cacheline; k > 0 ; k--){		// Left-most elements
					
					tmpcol = pattern->cols[j] - k;
					elementstoadd[(unsigned int) cacheline - k] = 1;
					
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
					if ((elementstoadd[k] == 1) && (addedelements < maxaddedelements) && (cacheadded < maxcacheelements) && ((k + pattern->cols[j] - (int) cacheline) < i) && ((k + pattern->cols[j] - (unsigned int) cacheline) >= 0)){
						expanded_patt->cols[counterB] = pattern->cols[j] + k - (unsigned int) cacheline;
						
						
						counterB++;
						addedelements++;
						cacheadded++;
					}
					else if ((elementstoadd[k] == 0) && ((k + pattern->cols[j] - (unsigned int) cacheline) <= i) && ((k + pattern->cols[j] - (unsigned int) cacheline) >= 0) && (pattern->cols[counterA] == (pattern->cols[j] + k - (unsigned int) cacheline))){
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


void ltp(mat_t *mat, pat_t *pattern, pat_t *expanded_patt){

	unsigned int i = 0, j = 0;
	unsigned int patterncount = 0;
	unsigned int counter = 0;

	for (i = 0; i < mat->size; i++){							// For every row in the matrix
		for (j = mat->rows[i]; j < mat->rows[i + 1]; j++){		// Gets only lower triangle values
			if (mat->cols[j] <= i){patterncount++;}
		}
		pattern->rows[i + 1] = patterncount;
	}
	pattern->nnz = patterncount;
	pattern->cols = malloc(pattern->nnz*sizeof(unsigned int));
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


void flt(mat_t *mat, pat_t *pattern, pat_t *expanded_patt){

	unsigned int i = 0, j = 0;
	unsigned int counter = 0, counter2 = 0;

	for (i = 0; i < mat->size; i++){							// For every row in the matrix
		pattern->rows[i + 1] = pattern->rows[i] + i + 1;
	}
	pattern->nnz = pattern->rows[mat->size];
	pattern->cols = malloc(pattern->nnz*sizeof(unsigned int));
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


void powerA(mat_t *mat, pat_t *pattern, pat_t *expanded_patt, double *xfinal){

	unsigned int i = 0, j = 0, k = 0, l = 0;
	unsigned int counter = 0;
	unsigned int patterncount = 0;
	unsigned int elements = 0;
	unsigned int dim = mat->size;

	
	pat_t *apower;
	apower = malloc(sizeof(mat_t));
	apower->rows = calloc((dim + 1), sizeof(unsigned int));
	apower->cols = calloc((20 * patternpower * mat->nnz), sizeof(unsigned int));
	
	pat_t *apowertemp;
	apowertemp = malloc(sizeof(mat_t));
	apowertemp->rows = calloc((dim + 1), sizeof(unsigned int));
	apowertemp->cols = calloc((20 * patternpower * mat->nnz), sizeof(unsigned int));
	
	unsigned int *pospower = calloc(dim, sizeof(unsigned int));
	unsigned int *poscol = calloc(dim, sizeof(unsigned int));
	
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
			qsort(poscol, elements, sizeof(unsigned int), comp);
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
	pattern->cols = malloc(pattern->nnz*sizeof(unsigned int));

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
	
	
	
	
	optexpandpattern(mat->size, pattern, expanded_patt, xfinal, mat);
	
// 	expandpattern(mat->size, pattern, expanded_patt, xfinal);
	
	
	
}


void perclt(mat_t *mat, pat_t *pattern, pat_t *expanded_patt, double *xfinal){

	unsigned int i = 0, j = 0;
	unsigned int counter = 0;
	unsigned int patterncount = 0;


	// WE GET INITIAL PATTERN PROPERLY STORED

	for (i = 0; i < mat->size; i++){	// First get initial pattern
		for (j = mat->rows[i]; j < mat->rows[i + 1]; j++){
			if (mat->cols[j] <= i){patterncount++;}
		}
		pattern->rows[i + 1] = patterncount;
	}

	pattern->nnz = patterncount;
	pattern->cols = malloc(pattern->nnz*sizeof(unsigned int));

	for (i = 0; i < mat->size; i++){
		for (j = mat->rows[i]; j < mat->rows[i + 1]; j++){
			if (mat->cols[j] <= i){
				pattern->cols[counter] = mat->cols[j];
				counter++;
			}
		}
	}
	
	expandpattern(mat->size, pattern, expanded_patt, xfinal);
	
}
