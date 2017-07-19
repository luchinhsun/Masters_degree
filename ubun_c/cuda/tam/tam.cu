#include "head.h"

extern float *d_a_tam, *d_b_tam,*d_c_tam;
extern float *d_c_new_tam,*d_d_new_tam;
extern float *d_Ax, *d_Ay;
extern float *d_b, *d_ut, *d_Vt;

extern float *d_V1, *d_V2;
//extern float *d_b, *d_ut, *d_Vt;

__global__ void GPU_adi_x(float *d_V1, float *d_b, int j){
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        
	if(i<Np){
        if(j==0){
                d_b[i] = d_V1[i*Np+j] + r*(-d_V1[i*Np+j] + d_V1[i*Np+j+1]);
        }else if(j==Np-1){
                d_b[i] = d_V1[i*Np+j] + r*(d_V1[i*Np+j-1] - d_V1[i*Np+j]);
        }else{
                d_b[i] = d_V1[i*Np+j] + (r/2)*(d_V1[i*Np+j-1] - 2*d_V1[i*Np+j] + d_V1[i*Np+j+1]);
        }
	}

}

__global__ void GPU_getV2(float *d_V2, float *d_Vt, int j){
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if(i<Np){
		d_V2[i*Np+j] = d_Vt[i];
	}

}

__global__ void GPU_adi_y(float *d_V2, float *d_b, int i){
        int j = blockDim.x * blockIdx.x + threadIdx.x;

	if(j<Np){
	if(i==0){
        	d_b[j] = d_V2[i*Np+j] + (r/2)*(-2*d_V2[i*Np+j] + 2*d_V2[(i+1)*Np+j]);
        }else if(i==Np-1){
                d_b[j] = d_V2[i*Np+j] + (r/2)*(2*d_V2[(i-1)*Np+j] - 2*d_V2[i*Np+j]);
        }else{
                d_b[j] = d_V2[i*Np+j] + (r/2)*(d_V2[(i-1)*Np+j] - 2*d_V2[i*Np+j] + d_V2[(i+1)*Np+j]);
        }
	}
}

__global__ void GPU_getV1(float *d_V1, float *d_ut, int i){
        int j = blockDim.x * blockIdx.x + threadIdx.x;

        if(j<Np){
                d_V1[i*Np+j] = d_ut[j];
        }

}

__global__ void GPU_tam(float *A, float *d_a_tam, float *d_b_tam, float *d_c_tam, int n){
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if(i==0){
		d_a_tam[0] = 0.0;
		d_c_tam[n-1] = 0.0;
	}
	__syncthreads();

	if(i<n-1){
		d_c_tam[i] = A[i*n+i+1];
		d_a_tam[i+1] = A[(i+1)*n+i];
	}

	if(i<n){
		d_b_tam[i] = A[i*n+i];
	}
}

__global__ void GPU_tam2(float *d_a_tam, float *d_b_tam, float *d_c_tam, float *d_c_new_tam, float *d_d_new_tam, float *d, float *x, int n){
	int i;
	d_c_new_tam[0] = d_c_tam[0]/d_b_tam[0];
        for(i=1;i<n-1;i++){
                d_c_new_tam[i] = d_c_tam[i]/(d_b_tam[i]-d_a_tam[i]*d_c_new_tam[i-1]);
        }

        d_d_new_tam[0] = d[0]/d_b_tam[0];
        for(i=1;i<n;i++){
                d_d_new_tam[i] = (d[i]-d_a_tam[i]*d_d_new_tam[i-1])/(d_b_tam[i]-d_a_tam[i]*d_c_new_tam[i-1]);
        }

        x[n-1] = d_d_new_tam[n-1];
        for(i=n-2;i>-1;i--){
                x[i] = d_d_new_tam[i]-d_c_new_tam[i]*x[i+1];
        }

}

void tam(){
	int tpb = 256;
	int bpg = (Np*Np+tpb-1)/tpb;
	int bpg1 = (Np+tpb-1)/tpb;

	//GPU_tam<<<bpg, tpb>>>(d_Ax, d_a_tam, d_b_tam, d_c_tam, Np);
	int j;
	for(j=0;j<Np;j++){
		GPU_adi_x<<<bpg, tpb>>>(d_V1, d_b, j);
		GPU_tam<<<bpg, tpb>>>(d_Ax, d_a_tam, d_b_tam, d_c_tam, Np);
		GPU_tam2<<<1, 1>>>(d_a_tam, d_b_tam, d_c_tam, d_c_new_tam, d_d_new_tam, d_b, d_Vt, Np);
		GPU_getV2<<<bpg1, tpb>>>(d_V2, d_Vt, j);
	}
}

void tam2(){
        int tpb = 256;
        int bpg = (Np*Np+tpb-1)/tpb;
	int bpg1 = (Np+tpb-1)/tpb;

        //GPU_tam<<<bpg, tpb>>>(d_Ay, d_a_tam, d_b_tam, d_c_tam, Np);
	int i;
	for(i=0;i<Np;i++){
		GPU_adi_y<<<bpg, tpb>>>(d_V2, d_b, i);
		GPU_tam<<<bpg, tpb>>>(d_Ay, d_a_tam, d_b_tam, d_c_tam, Np);
        	GPU_tam2<<<1, 1>>>(d_a_tam, d_b_tam, d_c_tam, d_c_new_tam, d_d_new_tam, d_b, d_ut, Np);
		GPU_getV1<<<bpg1, tpb>>>(d_V1, d_ut, i);
	}
}

