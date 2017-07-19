#include "head.h"

void fsource(float *V1, float t, float *x, float *y, float *F){

	int i,j;

	for(i=0;i<Np;i++){
                for(j=0;j<Np;j++){
                        F[i*Np+j] = exp(-2*t)*cos(M_PI*x[i])*cos(M_PI*y[j])*(2*M_PI*M_PI)-2*V1[i*Np+j];
                }
        }
}
