#include "head.h"

int k_t = (tfinal/dt);

float *Ax, *Ay, *V1, *V2, *W, *F, *V_tmp, *W_tmp, *ue;
float *b, *x, *y, *ut, *Vt;
float *t;

//GPU variable
float *d_Ax, *d_Ay, *d_V1, *d_V2;
float *d_b, *d_ut, *d_Vt;

//GPU tam variable
float *d_a_tam, *d_b_tam,*d_c_tam;
float *d_c_new_tam,*d_d_new_tam;

//GPU RK variable
float *d_F, *d_x, *d_y, *d_V_tmp;
float *d_t;

clock_t start;
clock_t end;
float time_used;

void Allocate(){

	size_t size;
	size = Np*Np*sizeof(float);

	Ax = (float*)malloc(size);
	Ay = (float*)malloc(size);
	V1 = (float*)malloc(size);
	V2 = (float*)malloc(size);
	W = (float*)malloc(size);
	F = (float*)malloc(size);
	V_tmp = (float*)malloc(size);
	W_tmp = (float*)malloc(size);
	ue = (float*)malloc(size);

	cudaError_t Error;

        Error = cudaMalloc((void**)&d_Ax, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_Ax) = %s\n", cudaGetErrorString(Error));
        Error = cudaMalloc((void**)&d_Ay, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_Ay) = %s\n", cudaGetErrorString(Error));
        Error = cudaMalloc((void**)&d_V1, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_V1) = %s\n", cudaGetErrorString(Error));
        Error = cudaMalloc((void**)&d_V2, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_V2) = %s\n", cudaGetErrorString(Error));
        Error = cudaMalloc((void**)&d_V_tmp, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_V_tmp) = %s\n", cudaGetErrorString(Error));
        Error = cudaMalloc((void**)&d_F, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_F) = %s\n", cudaGetErrorString(Error));

	size = Np*sizeof(float);

	b = (float*)malloc(size);
	x = (float*)malloc(size);
	y = (float*)malloc(size);
	ut = (float*)malloc(size);
	Vt = (float*)malloc(size);

        Error = cudaMalloc((void**)&d_b, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_b) = %s\n", cudaGetErrorString(Error));
        Error = cudaMalloc((void**)&d_x, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_x) = %s\n", cudaGetErrorString(Error));
        Error = cudaMalloc((void**)&d_y, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_y) = %s\n", cudaGetErrorString(Error));
        Error = cudaMalloc((void**)&d_ut, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_ut) = %s\n", cudaGetErrorString(Error));
        Error = cudaMalloc((void**)&d_Vt, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_Vt) = %s\n", cudaGetErrorString(Error));

	Error = cudaMalloc((void**)&d_a_tam, size);
	if(Error != cudaSuccess)
	printf("CUDA error(malloc d_a_tam) = %s\n", cudaGetErrorString(Error));
        Error = cudaMalloc((void**)&d_b_tam, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_b_tam) = %s\n", cudaGetErrorString(Error));
        Error = cudaMalloc((void**)&d_c_tam, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_c_tam) = %s\n", cudaGetErrorString(Error));
        Error = cudaMalloc((void**)&d_c_new_tam, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_c_new_tam) = %s\n", cudaGetErrorString(Error));
        Error = cudaMalloc((void**)&d_d_new_tam, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_d_new_tam) = %s\n", cudaGetErrorString(Error));

	size = 1*sizeof(float);
	t = (float*)malloc(size);
        Error = cudaMalloc((void**)&d_t, size);
        if(Error != cudaSuccess)
        printf("CUDA error(malloc d_t) = %s\n", cudaGetErrorString(Error));
}

void Save_Result(){

	FILE *pFile;
	int i,j;
	int index;
	int n;
	n = Np;
	pFile = fopen("V1.txt","w+");
	// Save the matrix A
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
/*
    pFile = fopen("F.txt","w+");
    // Save the matrix A
    for (i = 0; i < n; i++) {
    	for (j = 0; j < n; j++) {
        	index = i*n + j;
            fprintf(pFile, "%g", F[index]);
            if (j == (n-1)) {
        		fprintf(pFile, "\n");
            }else{
            	fprintf(pFile, "\t");
            }
    	}
	}
    fclose(pFile);
    pFile = fopen("V_tmp.txt","w+");
    // Save the matrix A
    for (i = 0; i < n; i++) {
    	for (j = 0; j < n; j++) {
        	index = i*n + j;
            fprintf(pFile, "%g", V_tmp[index]);
            if (j == (n-1)) {
        		fprintf(pFile, "\n");
            }else{
            	fprintf(pFile, "\t");
            }
    	}
	}
    fclose(pFile);
    pFile = fopen("b.txt","w+");
    for (i = 0; i < n; i++) {
        fprintf(pFile, "%g", b[i]);
        fprintf(pFile, "\t");
    }
    fclose(pFile);
    pFile = fopen("x.txt","w+");
    for (i = 0; i < n; i++) {
        fprintf(pFile, "%g", x[i]);
        fprintf(pFile, "\t");
    }
    fclose(pFile);
	*/
}

void Free(){
	free(Ax);free(Ay);free(V1);free(V2);
	free(W);free(F);free(V_tmp);free(W_tmp);free(ue);
	free(b);free(x);free(y);free(ut);free(Vt);
	free(t);

	cudaFree(d_Ax);cudaFree(d_Ay);cudaFree(d_V1);cudaFree(d_V2);
	cudaFree(d_V_tmp);cudaFree(d_F);
	cudaFree(d_t);
	cudaFree(d_b);cudaFree(d_x);cudaFree(d_y);cudaFree(d_Vt);cudaFree(d_ut);
	cudaFree(d_a_tam);cudaFree(d_b_tam);cudaFree(d_c_tam);
	cudaFree(d_c_new_tam);cudaFree(d_d_new_tam);
}

void Send_to_Device(){
	start = clock();

	cudaError_t Error;
	size_t size;
	size = Np*Np*sizeof(float);
	Error = cudaMemcpy(d_V1, V1, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy V1->d_V1) = %s\n",cudaGetErrorString(Error));
/*
        Error = cudaMemcpy(d_V2, V2, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy V2->d_V2) = %s\n",cudaGetErrorString(Error));
*/
	Error = cudaMemcpy(d_Ax, Ax, size, cudaMemcpyHostToDevice);
	if (Error != cudaSuccess)
	printf("CUDA error(copy Ax->d_Ax) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_Ay, Ay, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy Ay->d_Ay) = %s\n",cudaGetErrorString(Error));

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

/*
	size = Np*sizeof(float);
        Error = cudaMemcpy(d_b, b, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy b->d_b) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_ut, ut, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy ut->d_ut) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_Vt, Vt, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy Vt->d_Vt) = %s\n",cudaGetErrorString(Error));
*/
}

void Send_to_Host(){
	cudaError_t Error;
        size_t size;
	size = Np*Np*sizeof(float);
	Error = cudaMemcpy(V1, d_V1, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_V1->V1) = %s\n",cudaGetErrorString(Error));
/*
	Error = cudaMemcpy(V2, d_V2, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_V2->V2) = %s\n",cudaGetErrorString(Error));

	Error = cudaMemcpy(V_tmp, d_V_tmp, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_V_tmp->V_tmp) = %s\n",cudaGetErrorString(Error));

        Error = cudaMemcpy(Ax, d_Ax, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_Ax->Ax) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(Ay, d_Ay, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_Ay->Ay) = %s\n",cudaGetErrorString(Error));
        size = Np*sizeof(float);
        Error = cudaMemcpy(b, d_b, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_b->b) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(ut, d_ut, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_ut->ut) = %s\n",cudaGetErrorString(Error));
        Error = cudaMemcpy(Vt, d_Vt, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_Vt->Vt) = %s\n",cudaGetErrorString(Error));
*/

	end = clock();
        time_used = (float)(end - start)/ CLOCKS_PER_SEC;
        printf("\ntime in cu = %f\n",time_used);
}

