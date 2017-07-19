#include "head.h"

void tam(float *A, float *d, float *x, int n){

	float *a, *b, *c;
	float *c_new, *d_new;

	size_t size;
	size = n*sizeof(float);

	a = (float*)malloc(size);
	b = (float*)malloc(size);
	c = (float*)malloc(size);
	c_new = (float*)malloc(size);
	d_new = (float*)malloc(size);

	a[0] = 0.0;
	c[n-1] = 0.0;

	int i;
	for(i=0;i<n-1;i++){
		c[i] = A[i*n+i+1];
		a[i+1] = A[(i+1)*n+i];
	}

	for(i=0;i<n;i++){
		b[i] = A[i*n+i];
	}

	c_new[0] = c[0]/b[0];
	for(i=1;i<n-1;i++){
		c_new[i] = c[i]/(b[i]-a[i]*c_new[i-1]);
	}

	d_new[0] = d[0]/b[0];
	for(i=1;i<n;i++){
		d_new[i] = (d[i]-a[i]*d_new[i-1])/(b[i]-a[i]*c_new[i-1]);
	}

	x[n-1] = d_new[n-1];
	for(i=n-2;i>-1;i--){
		x[i] = d_new[i]-c_new[i]*x[i+1];
	}

	free(a);free(b);free(c);
	free(c_new);free(d_new);
}
