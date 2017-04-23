#include "head.h"

extern float *a_tam, *b_tam, *c_tam;
extern float *c_new_tam, *d_new_tam;

void tam(float *A, float *d, float *x, int n){
	/*
	float *a, *b, *c;
	float *c_new, *d_new;

	size_t size;
	size = n*sizeof(float);

	a = (float*)malloc(size);
	b = (float*)malloc(size);
	c = (float*)malloc(size);
	c_new = (float*)malloc(size);
	d_new = (float*)malloc(size);
	*/

	#pragma omp single
	{
	a_tam[0] = 0.0;
	c_tam[n-1] = 0.0;
	}

	int i;
	#pragma omp for private(i) nowait
	for(i=0;i<n-1;i++){
		c_tam[i] = A[i*n+i+1];
		a_tam[i+1] = A[(i+1)*n+i];
	}

	#pragma omp for private(i)
	for(i=0;i<n;i++){
		b_tam[i] = A[i*n+i];
	}

	#pragma omp single
	{
	c_new_tam[0] = c_tam[0]/b_tam[0];
	}

	#pragma omp for private(i)
	for(i=1;i<n-1;i++){
		c_new_tam[i] = c_tam[i]/(b_tam[i]-a_tam[i]*c_new_tam[i-1]);
	}

	#pragma omp single
	{
	d_new_tam[0] = d[0]/b_tam[0];
	for(i=1;i<n;i++){
		d_new_tam[i] = (d[i]-a_tam[i]*d_new_tam[i-1])/(b_tam[i]-a_tam[i]*c_new_tam[i-1]);
	}

	x[n-1] = d_new_tam[n-1];
	for(i=n-2;i>-1;i--){
		x[i] = d_new_tam[i]-c_new_tam[i]*x[i+1];
	}
	}

	/*
	free(a);free(b);free(c);
	free(c_new);free(d_new);
	*/
}
