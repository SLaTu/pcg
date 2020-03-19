 #include "cg_basics.h"



double getvalue_mat(int row, int col, mat_t *mat){

	double value = 0.0;
	int i = 0;

	for (i = mat->rows[row]; i < mat->rows[row+1]; i++){
		if (mat->cols[i] == col) value = mat->values[i];
	}
	return value;
}


double getvalue_matCSC(int row, int col, mat_t *mat){

	double value = 0.0;
	int i = 0;

	for (i = mat->cols[col]; i < mat->cols[col + 1]; i++){
		if (mat->rows[i] == row) value = mat->values[i];
	}
	return value;
}


void pmultMatVect(double *tmp, double *vect, mat_t *mat){
	int i = 0, j = 0;
	double aux;
	#pragma omp parallel for private(j, aux)
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
	#pragma omp parallel for private(j, temp, ind)
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

void pmultMatVectCOO(double *tmp, double *vect, mat_t *mat, int *limits, int nthreads){
	int i = 0, j = 0;
	#pragma omp parallel for private(j, i)
	for (i = 0; i < nthreads; i++){
		for (j = limits[i]; j < limits[i + 1]; j++){
			tmp[mat->rows[j]] += mat->values[j] * vect[mat->cols[j]];
		}
	}
}

void pmultMatVect_DUMM(double *tmp, mat_t *mat){
	int i = 0, j = 0;
	#pragma omp parallel for private(j)
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
	#pragma omp parallel for private(j, temp, ind)
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

void pmultMatVectCOO_DUMM(double *tmp, mat_t *mat, int *limits, int nthreads){
	int i = 0, j = 0;
	#pragma omp parallel for private(j, i)
	for (i = 0; i < nthreads; i++){
		for (j = limits[i]; j < limits[i + 1]; j++){
			tmp[mat->rows[j]] += mat->values[j] * mat->cols[j];
		}
	}
}

void pcopyvect(double *IN, double *OUT, int max){
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++) OUT[i] = IN[i];
}

void pzerovect(double *vect, int max){
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++){
		vect[i] = 0.0;
	}
}

double pmultVectVect(double *V1, double *V2, int max){
	
	double result = 0.0; int i;
	#pragma omp parallel for reduction(+:result)
	for (i = 0; i < max; i++) result += V1[i] * V2[i];
	return result;
}

void pscaleVect(double scalar, double *V1, double *tmp, int max){	
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++) {tmp[i] = 0.0; tmp[i] = scalar * V1[i];}
}

void paddToVect(double *V1, double *V2, int max){
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++){V1[i] += V2[i];}
}

void psubToVect(double *V1, double *V2, int max){
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++){V1[i] -= V2[i];}
}

void paddVects(double *V1, double *V2, double *Vres, int max){
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++){Vres[i] = V1[i] + V2[i];}
}

void psubtVects(double *V1, double *V2, double *Vres, int max){
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++){Vres[i] = V1[i] - V2[i];}
}

void pequalvects(double *V1, double *V2, int max){
	int i = 0;
	#pragma omp parallel for
	for (i = 0; i < max; i++) {V1[i] = V2[i];}
}


void cleararrays(double *r, double *rold, double *d, double *q, double *s, double *tmp, double *x, int dim){
	#pragma omp parallel for
	for (int i = 0; i < dim; i++){
		r[i] = 0.0;
		rold[i] = 0.0;
		d[i] = 0.0;
		q[i] = 0.0;
		s[i] = 0.0;
		tmp[i] = 0.0;
		x[i] = 0.0;
	}
}


void clearlogs(double *elapses, double *residuals, int imax){
	#pragma omp parallel for
	for (int i = 0; i < imax; i++){
		elapses[i] = 0.0;
		residuals[i] = 0.0;
	}
}




