 

#include "cg_basics.h"

void multMatVect(double *tmp, double *vect, mat_t *mat){
	int i = 0, j = 0;
	for (i = 0; i < mat->size; i++){
		tmp[i] = 0.0;
		for (j = mat->rows[i]; j < mat->rows[i + 1]; j++)
		{
			tmp[i] += mat->values[j] * vect[mat->cols[j]];
		}
	}
}

void multMatVectCSC(double *tmp, double *vect, mat_t *mat){
	int i = 0, j = 0;
	for(i = 0; i < mat->size; i++) tmp[i] = 0.0;
	for (i = 0; i < mat->size; i++){
		for (j = mat->cols[i]; j < mat->cols[i + 1]; j++)
		{
			tmp[mat->rows[j]] += mat->values[j] * vect[i];
		}
	}
}

void zerovect(double *vect, unsigned int max){
	unsigned int i = 0;
	for (i = 0; i < max; i++){
		vect[i] = 0.0;
	}
}

double multVectVect(double *V1, double *V2, unsigned int max){
	
	double result = 0.0; int i;
	for (i = 0; i < max; i++) result += V1[i] * V2[i];
	return result;
}

void scaleVect(double scalar, double *V1, double *tmp, unsigned int max){	
	int i = 0;
	for (i = 0; i < max; i++) {tmp[i] = 0.0; tmp[i] = scalar * V1[i];}
}

void addToVect(double *V1, double *V2, unsigned int max){
	int i = 0;
	for (i = 0; i < max; i++){V1[i] += V2[i];}
}

void subToVect(double *V1, double *V2, unsigned int max){
	int i = 0;
	for (i = 0; i < max; i++){V1[i] -= V2[i];}
}

void addVects(double *V1, double *V2, double *Vres, unsigned int max){
	int i = 0;
	for (i = 0; i < max; i++){Vres[i] = 0.0; Vres[i] = V1[i] + V2[i];}
}

void subtVects(double *V1, double *V2, double *Vres, unsigned int max){
	int i = 0;
	for (i = 0; i < max; i++){Vres[i] = 0.0; Vres[i] = V1[i] - V2[i];}
}

double getvalue_mat(unsigned int row, unsigned int col, mat_t *mat){

	double value = 0.0;
	unsigned int i = 0;

	for (i = mat->rows[row]; i < mat->rows[row+1]; i++){
		if (mat->cols[i] == col) value = mat->values[i];
	}
	return value;
}























/* PARALLEL VERSIONS */


void pmultMatVect(double *tmp, double *vect, mat_t *mat){
	int i = 0, j = 0;
	double aux;
	#pragma omp parallel for schedule(runtime) private(j, aux)
	for (i = 0; i < mat->size; i++){
		aux = 0.0;
		for (j = mat->rows[i]; j < mat->rows[i + 1]; j++)
		{
			aux += mat->values[j] * vect[mat->cols[j]];
		}
		tmp[i] = aux;
	}
}

void pmultMatVectCSC(double *tmp, double *vect, mat_t *mat){
	int i = 0, j = 0;
	double temp = 0.0;
	int ind = 0;
	#pragma omp parallel for schedule(runtime) private(j, temp, ind)
	for (i = 0; i < mat->size; i++){
		double vecti = vect[i];
		for (j = mat->cols[i]; j < mat->cols[i + 1]; j++)
		{
			ind = mat->rows[j];
			temp = mat->values[j] * vecti;
			#pragma omp atomic
			tmp[ind] += temp;
		}
	}
}



void pmultMatVect_DUMM(double *tmp, mat_t *mat){
	int i = 0, j = 0;
	#pragma omp parallel for schedule(runtime) private(j)
	for (i = 0; i < mat->size; i++){
		for (j = mat->rows[i]; j < mat->rows[i + 1]; j++)
		{
			tmp[i] += mat->values[j] * mat->cols[j];
		}
	}
}

void pmultMatVectCSC_DUMM(double *tmp, mat_t *mat){
	int i = 0, j = 0;
	double temp = 0.0;
	int ind = 0;
	#pragma omp parallel for schedule(runtime) private(j, temp, ind)
	for (i = 0; i < mat->size; i++){
		for (j = mat->cols[i]; j < mat->cols[i + 1]; j++)
		{
			ind = mat->rows[j];
			temp = mat->values[j] * i;
			#pragma omp atomic
			tmp[ind] += temp;
		}
	}
}

void pzerovect(double *vect, unsigned int max){
	unsigned int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++){
		vect[i] = 0.0;
	}
}

double pmultVectVect(double *V1, double *V2, unsigned int max){
	
	double result = 0.0; int i;
	#pragma omp parallel for reduction(+:result)
	for (i = 0; i < max; i++) result += V1[i] * V2[i];
	return result;
}

void pscaleVect(double scalar, double *V1, double *tmp, unsigned int max){	
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++) {tmp[i] = 0.0; tmp[i] = scalar * V1[i];}
}

void paddToVect(double *V1, double *V2, unsigned int max){
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++){V1[i] += V2[i];}
}

void psubToVect(double *V1, double *V2, unsigned int max){
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++){V1[i] -= V2[i];}
}

void paddVects(double *V1, double *V2, double *Vres, unsigned int max){
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++){Vres[i] = 0.0; Vres[i] = V1[i] + V2[i];}
}

void psubtVects(double *V1, double *V2, double *Vres, unsigned int max){
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++){Vres[i] = 0.0; Vres[i] = V1[i] - V2[i];}
}

void pequalvects(double *V1, double *V2, unsigned int max){
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++) {V1[i] = V2[i];}
}









