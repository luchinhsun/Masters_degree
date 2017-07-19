#include "head.h"

extern float *b, *ut, *Vt;
extern float *d_V1, *d_V2;
extern float *d_b, *d_ut, *d_Vt;

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
void ADI1(int j){
        int tpb = 256;
        int bpg = (Np+tpb-1)/tpb;

        GPU_adi_x<<<bpg, tpb>>>(d_V1, d_b, j);
}
void ADI1_2(int j){
	int tpb = 256;
        int bpg = (Np+tpb-1)/tpb;

        GPU_getV2<<<bpg, tpb>>>(d_V2, d_Vt, j);
}

void ADI2(int i){
        int tpb = 256;
        int bpg = (Np+tpb-1)/tpb;

        GPU_adi_y<<<bpg, tpb>>>(d_V2, d_b, i);
}

void ADI2_2(int i){
	int tpb = 256;
        int bpg = (Np+tpb-1)/tpb;

        GPU_getV1<<<bpg, tpb>>>(d_V1, d_ut, i);
}

void Send_to_D(){
	size_t size;
	size = Np*sizeof(float);
	cudaMemcpy(d_Vt, Vt, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_ut, ut, size, cudaMemcpyHostToDevice);
}

void Send_to_H(){
	size_t size;
        size = Np*sizeof(float);
	cudaMemcpy(b, d_b, size, cudaMemcpyDeviceToHost);
}
