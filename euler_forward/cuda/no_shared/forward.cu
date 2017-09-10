#include "head.h"

extern float *d_F, *d_x, *d_y;
extern float *d_V1, *d_V2;
extern float *d_t;

__global__ void GPU_fsource(float *d_V2, float *d_t, float *d_x, float *d_y, float *d_F){
        int k = blockDim.x * blockIdx.x + threadIdx.x;
	int i, j;

        if(k<Np*Np){
		i = k/Np;
		j = k%Np;
		d_F[i*Np+j] = exp(-2.0*d_t[0])*cos(M_PI*d_x[i])*
                                                cos(M_PI*d_y[j])*(2.0*M_PI*M_PI)
                                                        -2.0*d_V2[i*Np+j];
        }
}

__global__ void forward(float *d_V1, float *d_V2, float *d_F){
	int k = blockDim.x * blockIdx.x + threadIdx.x;
	int i = k/Np;
        int j = k%Np;

	if(k<Np*Np){
		d_V2[i*Np+j] = d_V1[(i+1)*(Np+2)+j+1] + dt/h1*
                                (d_V1[i*(Np+2)+j+1]+d_V1[(i+2)*(Np+2)+j+1]+d_V1[(i+1)*(Np+2)+j]+d_V1[(i+1)*(Np+2)+j+2]
                                        -4*d_V1[(i+1)*(Np+2)+j+1]) + dt*d_F[i*Np+j];
	}
}

__global__ void update_V1(float *d_V1, float *d_V2, float *d_t){
	int k = blockDim.x * blockIdx.x + threadIdx.x;
	int i, j;

	if(k<Np*Np){
		i = k/Np;
		j = k%Np;
		d_V1[(i+1)*(Np+2)+j+1] = d_V2[i*Np+j];
	}

	if(k<Np){
		d_V1[k+1] = d_V2[(Np-1)*Np+k];
                d_V1[(Np+1)*(Np+2)+k+1] = d_V2[k];
                d_V1[(k+1)*(Np+2)] = d_V2[k*Np+Np-1];
                d_V1[(k+1)*(Np+2)+Np+1] = d_V2[k*Np];
	}

	if(k==0){
                d_t[0] = d_t[0]+dt;
        }

}

void forward_euler(){
	int tpb = 256;
        int bpg = (Np*Np+tpb-1)/tpb;

	update_V1<<<bpg, tpb>>>(d_V1, d_V2, d_t);
	
	forward<<<bpg, tpb>>>(d_V1, d_V2, d_F);

	GPU_fsource<<<bpg, tpb>>>(d_V2, d_t, d_x, d_y, d_F);
	//cudaDeviceSynchronize();
}
