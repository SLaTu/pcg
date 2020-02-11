
#include "cg_patterns.h"

extern unsigned int percentpattern;
extern unsigned int patternpower;
extern FILE *outmn;

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

void expandpattern(unsigned int dim, pat_t *pattern, pat_t *expanded_patt, double *xfinal){
	
	unsigned int i = 0, j = 0, k = 0, l = 0, m = 0;
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
	
	int lim1 = 100;
	int lim2 = 110;

// 	for (i = lim1; i < lim2; i++) printf("%u, ", pattern->rows[i] - pattern->rows[i - 1]);
// 	printf("\n");

	for(i = 0; i < dim; i++){

		oldcachecol = -10;
		oldcacheline = 0;
		
// 		if ((i > lim1) && (i < lim2)) printf("\n");
		
		for(j = pattern->rows[i]; j < pattern->rows[i + 1]; j++){
			
			cachecol = pattern->cols[j];
			cacheline = (((long long) &xfinal[pattern->cols[j]]%64)/8);
// 			if ((i > lim1) && (i < lim2)) printf("\n");
// 			if ((i > lim1) && (i < lim2)) printf("CACHELINE: %lli\tCOLUMN: %u\t\t ", cacheline, cachecol);
			
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
				
				
// 				if ((i > lim1) && (i < lim2)) for (m = 0; m < 8; m++) printf("%u ", elementstoadd[m]);
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
// 				if ((i > lim1) && (i < lim2)) printf("\t");
				
				calc += cacheelements * ((int) percentpattern); // si > 75 -> 1, si > 175 -> 2
				maxcacheelements = 0;
				
				while (calc > 66){
					calc -= 100;
					maxcacheelements++;
				}
				
// 				if ((maxcacheelements + cacheelements) > 8) {
// 					calc += ((maxcacheelements + cacheelements) - 8) * 100;
// 				}
				
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
				initk = (int) cacheline - (int) maxcacheelements;
				if (initk < 0) initk = 0;
				
// 				for (k = 0; k < 8; k++){
				for (k = initk; k < 8; k++){
					if ((elementstoadd[k] == 1) && (addedelements < maxaddedelements) && (cacheadded < maxcacheelements) && ((k + pattern->cols[j] - (int) cacheline) < i) && ((k + pattern->cols[j] - (unsigned int) cacheline) >= 0)){
						expanded_patt->cols[counterB] = pattern->cols[j] + k - (unsigned int) cacheline;
						
// 						if ((i > lim1) && (i < lim2)) printf("%u\t", expanded_patt->cols[counterB]);
						
						counterB++;
						addedelements++;
						cacheadded++;
					}
					else if ((elementstoadd[k] == 0) && ((k + pattern->cols[j] - (unsigned int) cacheline) <= i) && ((k + pattern->cols[j] - (unsigned int) cacheline) >= 0) && (pattern->cols[counterA] == (pattern->cols[j] + k - (unsigned int) cacheline))){
						expanded_patt->cols[counterB] = pattern->cols[counterA];
						
// 						if ((i > lim1) && (i < lim2)) printf("|%u|\t", expanded_patt->cols[counterB]);
						
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
	
// 	printf("\n");
// 	printf("\n");
// 	for (i = lim1 + 1; i < lim2; i++){
// 		for (j = expanded_patt->rows[i]; j < expanded_patt->rows[i + 1]; j++){
// 			printf("%u,  ", expanded_patt->cols[j]);
// 		}
// 		printf("\n");
// 	}
	
	printf("\n\n%u -> ^%.2lf%%\n\n", expanded_patt->nnz, 100 * (((double) expanded_patt->nnz - (double) pattern->nnz)/ ((double) pattern->nnz)));
	
	
// 	fprintf(outmn, "%u\t", expanded_patt->nnz);
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
	
// 	int limone = 9990;
// 	int limtwo = 10000;
// 	
// 	
// 	for (i = limone; i < limtwo; i++){
// 		for (j = apower->rows[i]; j < apower->rows[i + 1]; j++){
// 			printf("%u, ", apower->cols[j]);
// 		}
// 		printf("\t%u\n", apower->rows[i + 1]);
// 	}
	
	
	
	
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
	
	printf("For Time: %lf\tSorting Time: %lf\tEnd Time: %lf\tTotal Time: %lf\n\n", tfor, tsort, tend, ttotal);
	
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
	
	
// 	for (i = limone; i < limtwo; i++){
// 		for (j = pattern->rows[i]; j < pattern->rows[i + 1]; j++){
// 			printf("%u, ", pattern->cols[j]);
// 		}
// 		printf("\t%u\n", pattern->rows[i + 1]);
// 	}
	
	
	
	expandpattern(mat->size, pattern, expanded_patt, xfinal);
	
	
// 	for (i = limone; i < limtwo; i++){
// 		for (j = expanded_patt->rows[i]; j < expanded_patt->rows[i + 1]; j++){
// 			printf("%u, ", expanded_patt->cols[j]);
// 		}
// 		printf("\t%u\n", expanded_patt->rows[i + 1]);
// 	}
	
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
