#include "head.h"

#define block_dim 32

extern float *d_F, *d_x, *d_y;
extern float *d_V1, *d_V2;
extern float *d_t;

__global__ void GPU_fsource(float *d_V1, float *d_t, float *d_x, float *d_y, float *d_F){

        int k = blockDim.x * blockIdx.x + threadIdx.x;
        int i, j;

        if(k<Np*Np){
		i = k/Np;
                j = k%Np;
                d_F[i*Np+j] = exp(-2.0*d_t[0])*cos(M_PI*d_x[i])*
                                                cos(M_PI*d_y[j])*(2.0*M_PI*M_PI)
                                                        -2.0*d_V1[i*Np+j];
        }
}

__global__ void forward(float *d_V1, float *d_V2, float *d_F){
	int y = gridDim.x * blockDim.x * (threadIdx.y + blockIdx.y * blockDim.y);
        int x = threadIdx.x + blockIdx.x * blockDim.x;
        int id = x + y;
	int newid = id - (threadIdx.y + blockIdx.y * blockDim.y)*(gridDim.x * blockDim.x - Np);
	int j = threadIdx.x;
	int i = threadIdx.y;
	int a = newid/Np;
	int b = newid%Np;

	__shared__ float V1[block_dim][block_dim];

	V1[i][j] = d_V1[newid];
	__syncthreads();

	if (newid<Np*Np){

		if(i>0&&i<(block_dim-1)&&j>0&&j<(block_dim-1)&&a>0&&a<Np-1&&b>0&&b<Np-1){
			d_V2[newid] = V1[i][j]+dt/h1*(V1[i-1][j]+V1[i+1][j]+V1[i][j-1]+V1[i][j+1]
					-4.0*V1[i][j])+dt*d_F[newid];
		}
		else if(a>0&&a<Np-1&&b>0&&b<Np-1){
			d_V2[a*Np+b] = d_V1[a*Np+b] + dt/h1*(d_V1[(a+1)*Np+b]+d_V1[(a-1)*Np+b]
					+d_V1[a*Np+b+1]+d_V1[a*Np+b-1]-4.0*d_V1[a*Np+b])
						+dt*d_F[a*Np+b];
		}
	}
}

__global__ void forward_bound(float *d_V1, float *d_V2, float *d_F){
	int k = blockDim.x * blockIdx.x + threadIdx.x;
	int i;

	if(k<(Np-2)){
		i = k+1;
		d_V2[0*Np+i] = d_V1[0*Np+i] + dt/h1*
                                (d_V1[(0+1)*Np+i]+d_V1[(Np-1)*Np+i]+d_V1[0*Np+i+1]+d_V1[0*Np+i-1]
                                        -4*d_V1[0*Np+i]) + dt*d_F[0*Np+i];
		d_V2[(Np-1)*Np+i] = d_V1[(Np-1)*Np+i] + dt/h1*
                                (d_V1[0*Np+i]+d_V1[(Np-1-1)*Np+i]+d_V1[(Np-1)*Np+i+1]+d_V1[(Np-1)*Np+i-1]
                                        -4*d_V1[(Np-1)*Np+i]) + dt*d_F[(Np-1)*Np+i];
		d_V2[i*Np+0] = d_V1[i*Np+0] + dt/h1*
                                (d_V1[(i+1)*Np+0]+d_V1[(i-1)*Np+0]+d_V1[i*Np+0+1]+d_V1[i*Np+Np-1]
                                        -4*d_V1[i*Np+0]) + dt*d_F[i*Np+0];
		d_V2[i*Np+(Np-1)] = d_V1[i*Np+(Np-1)] + dt/h1*
                                (d_V1[(i+1)*Np+(Np-1)]+d_V1[(i-1)*Np+(Np-1)]+d_V1[i*Np+(Np-1)-1]+d_V1[i*Np+0]
                                        -4*d_V1[i*Np+(Np-1)]) + dt*d_F[i*Np+(Np-1)];
	}
	if(k==(Np-2)){
		d_V2[0*Np+0] = d_V1[0*Np+0] + dt/h1*
                                (d_V1[(0+1)*Np+0]+d_V1[(Np-1)*Np+0]+d_V1[0*Np+0+1]+d_V1[0*Np+Np-1]
					-4*d_V1[0*Np+0]) + dt*d_F[0*Np+0];
		d_V2[(Np-1)*Np+0] = d_V1[(Np-1)*Np+0] + dt/h1*
                                (d_V1[0*Np+0]+d_V1[(Np-1-1)*Np+0]+d_V1[(Np-1)*Np+0+1]+d_V1[(Np-1)*Np+Np-1]
                                        -4*d_V1[(Np-1)*Np+0]) + dt*d_F[(Np-1)*Np+0];
		d_V2[0*Np+(Np-1)] = d_V1[0*Np+(Np-1)] + dt/h1*
                                (d_V1[(0+1)*Np+(Np-1)]+d_V1[(Np-1)*Np+(Np-1)]+d_V1[0*Np+0]+d_V1[0*Np+(Np-1)-1]
                                        -4*d_V1[0*Np+(Np-1)]) + dt*d_F[0*Np+(Np-1)];
		d_V2[(Np-1)*Np+(Np-1)] = d_V1[(Np-1)*Np+(Np-1)] + dt/h1*
                                (d_V1[(Np-1-1)*Np+(Np-1)]+d_V1[0*Np+(Np-1)]+d_V1[(Np-1)*Np+0]+d_V1[(Np-1)*Np+(Np-1)-1]
                                        -4*d_V1[(Np-1)*Np+(Np-1)]) + dt*d_F[(Np-1)*Np+(Np-1)];
	}
}

__global__ void update_V1(float *d_V1, float *d_V2, float *d_t){
	int k = blockDim.x * blockIdx.x + threadIdx.x;

	if(k<(Np*Np)){
		d_V1[k] = d_V2[k];
	}

	if(k==0){
                d_t[0] = d_t[0]+dt;
        }

}

void forward_euler(){
        int tpb = block_dim;
        int bpg = (Np+tpb-1)/tpb;
	dim3 threads(tpb, tpb, 1);
        dim3 grid(bpg, bpg);
	forward<<<grid, threads>>>(d_V1, d_V2, d_F);

	tpb = 256;
	bpg = ((Np-1)+tpb-1)/tpb;
	forward_bound<<<bpg, tpb>>>(d_V1, d_V2, d_F);
	bpg = (Np*Np+tpb-1)/tpb;
	update_V1<<<bpg, tpb>>>(d_V1, d_V2, d_t);

	bpg = (Np*Np+tpb-1)/tpb;
	GPU_fsource<<<bpg, tpb>>>(d_V1, d_t, d_x, d_y, d_F);
	//cudaDeviceSynchronize();
}
