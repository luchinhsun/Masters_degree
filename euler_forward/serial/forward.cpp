#include "head.h"

void forward(float *d_V1, float *d_V2, float *d_F){

	int i, j;
	for(i=0;i<Np;i++){
		for(j=0;j<Np;j++){
			d_V2[i*Np+j] = d_V1[(i+1)*(Np+2)+j+1] + dt/h1*
		                (d_V1[i*(Np+2)+j+1]+d_V1[(i+2)*(Np+2)+j+1]+d_V1[(i+1)*(Np+2)+j]+d_V1[(i+1)*(Np+2)+j+2]
                                        -4*d_V1[(i+1)*(Np+2)+j+1]) + dt*d_F[i*Np+j];
		}
	}
}

void update_V1(float *d_V1, float *d_V2){
	int i, j;

	for(i=0;i<Np;i++){
		for(j=0;j<Np;j++){
			d_V1[(i+1)*(Np+2)+j+1] = d_V2[i*Np+j];
		}
		d_V1[i+1] = d_V2[(Np-1)*Np+i];
		d_V1[(Np+1)*(Np+2)+i+1] = d_V2[i];
		d_V1[(i+1)*(Np+2)] = d_V2[i*Np+Np-1];
		d_V1[(i+1)*(Np+2)+Np+1] = d_V2[i*Np];
	}
}
