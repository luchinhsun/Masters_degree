#include "head.h"

int k_t = (tfinal/dt);

float *V1, *V2, *F, *ue;
float *x, *y;
float *t;
float *d_V1, *d_V2, *d_F;
float *d_x, *d_y;
float *d_t;

void Allocate(){

	size_t size;
	size = Np*Np*sizeof(float);

	V1 = (float*)malloc(size);
	V2 = (float*)malloc(size);
	F = (float*)malloc(size);
	ue = (float*)malloc(size);
	cudaMalloc((void**)&d_V1, size);
	cudaMalloc((void**)&d_V2, size);
	cudaMalloc((void**)&d_F, size);


	size = Np*sizeof(float);

	x = (float*)malloc(size);
	y = (float*)malloc(size);
	cudaMalloc((void**)&d_x, size);
	cudaMalloc((void**)&d_y, size);

	size = sizeof(float);
	t = (float*)malloc(size);
        cudaMalloc((void**)&d_t, size);
}

void Save_Result(){

	FILE *pFile;
	int i,j;
	int index;
	int n;
	n = Np;
	pFile = fopen("V1.txt","w+");
	// Save the matrix V1
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			index = i*n + j;
			fprintf(pFile, "%g", V1[index]);
			if (j == (n-1)) {
				fprintf(pFile, "\n");
			}else{
				fprintf(pFile, "\t");
			}
		}
	}
	fclose(pFile);
}

void Free(){
	free(V1);free(V2);free(F);free(ue);
	free(x);free(y);
	free(t);
	cudaFree(d_V1);cudaFree(d_V2);cudaFree(d_F);
	cudaFree(d_x);cudaFree(d_y);
	cudaFree(d_t);
}

void Send_to_Device(){
	cudaError_t Error;
	size_t size;
        size = Np*Np*sizeof(float);

	Error = cudaMemcpy(d_V1, V1, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy V1->d_V1) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_F, F, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy F->d_F) = %s\n",cudaGetErrorString(Error));

	size = Np*sizeof(float);
        Error = cudaMemcpy(d_x, x, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy x->d_x) = %s\n",cudaGetErrorString(Error));
        Error = cudaMemcpy(d_y, y, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy y->d_y) = %s\n",cudaGetErrorString(Error));

	size = 1*sizeof(float);
        Error = cudaMemcpy(d_t, t, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy t->d_t) = %s\n",cudaGetErrorString(Error));
}

void Send_to_Host(){
	cudaError_t Error;
        size_t size;
        size = Np*Np*sizeof(float);

	Error = cudaMemcpy(V1, d_V1, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_V1->V1) = %s\n",cudaGetErrorString(Error));

	size = 1*sizeof(float);
        Error = cudaMemcpy(t, d_t, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_t->t) = %s\n",cudaGetErrorString(Error));
}
