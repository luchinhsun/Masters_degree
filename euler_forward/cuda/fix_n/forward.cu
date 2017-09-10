#include "head.h"

#define block_dim 32

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
	int yid = threadIdx.y + blockIdx.y * (blockDim.y-2);
	int xid = threadIdx.x + blockIdx.x * (blockDim.x-2);
	int totalx = (gridDim.x-1) * (blockDim.x-2) + blockDim.x;
	int id = xid + totalx * yid;
	int V2id = (yid-1)* (totalx-2) + xid - 1;
	int j = threadIdx.x;
	int i = threadIdx.y;

	__shared__ float V1[block_dim][block_dim];
	V1[i][j] = d_V1[id];
	__syncthreads();

	if(i>0&&i<(block_dim-1)&&j>0&&j<(block_dim-1)){
		d_V2[V2id] = V1[i][j] + dt/h1*(V1[i-1][j]+V1[i+1][j]+V1[i][j-1]+V1[i][j+1]
                                        -4*V1[i][j]) + dt*d_F[V2id];
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

	int twodtpb = block_dim;
	int twodbpg = (Np/(twodtpb-2));
	dim3 threads(twodtpb, twodtpb, 1);
	dim3 grid(twodbpg, twodbpg);
	forward<<<grid, threads>>>(d_V1, d_V2, d_F);

	GPU_fsource<<<bpg, tpb>>>(d_V2, d_t, d_x, d_y, d_F);
	//cudaDeviceSynchronize();
}
