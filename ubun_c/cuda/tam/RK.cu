#include "head.h"

extern float *d_F, *d_x, *d_y;
extern float *d_V2;
extern float *d_V_tmp;
extern float *d_t;

__global__ void GPU_fsource(float *d_V_tmp, float *d_t, float *d_x, float *d_y, float *d_F){
        int k = blockDim.x * blockIdx.x + threadIdx.x;
        int j;

        if(k<Np){
                for(j=0;j<Np;j++){
                        d_F[k*Np+j] = exp(-2.0*d_t[0])*cos(M_PI*d_x[k])*
                                                cos(M_PI*d_y[j])*(2.0*M_PI*M_PI)
                                                        -2.0*d_V_tmp[k*Np+j];
                }
        }
}

__global__ void GPU_RKa(float *d_V_tmp, float *d_V2, float *d_F, float *d_t){
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if(i<Np*Np){
		d_V_tmp[i] = d_V2[i] + (1.0/2.0)*dt*d_F[i];
	}

	if(i==0)	d_t[0] = d_t[0]+dt/2.0;
}

__global__ void GPU_RKb(float *d_V2, float *d_F, float *d_t){
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if(i<Np*Np){
		d_V2[i] = d_V2[i] + dt*d_F[i];
        }

	if(i==0)	d_t[0] = d_t[0]+dt/2.0;
}

void RK(){
	int tpb = 256;
        int bpg = (Np*Np+tpb-1)/tpb;

	GPU_RKa<<<bpg, tpb>>>(d_V_tmp, d_V2, d_F, d_t);
        GPU_fsource<<<bpg, tpb>>>(d_V_tmp, d_t, d_x, d_y, d_F);
	GPU_RKb<<<bpg, tpb>>>(d_V2, d_F, d_t);
}

