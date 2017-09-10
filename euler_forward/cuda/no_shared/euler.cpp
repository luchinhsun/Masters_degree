/*
	creat by Chinhsun.Lu
*/

#include "head.h"

extern int k_t;

extern float *V1, *V2, *F, *ue;
extern float *x, *y;
extern float *t;

int main(){

	//Allocate and initial variable
	Allocate();
	Init();

	int i, j, k;

	//Time t
	t[0] = 0.0;

	struct timeb start, end;
        int diff;

	ftime(&start);
	Send_to_Device();
	//---------- Big loop for time t --------------------------------------
	for(k=0;k<k_t;k++){
		forward_euler();
	}
	Send_to_Host();
	ftime(&end);
        diff = (int)(1000.0*(end.time-start.time)+(end.millitm-start.millitm));
        printf("\nTime = %d ms\n", diff);

	for(i=0;i<Np;i++){
		for(j=0;j<Np;j++){
			ue[i*Np+j] = exp(-2.0*tfinal)*cos(M_PI*x[i])*cos(M_PI*y[j]);
		}
	}

	float e = 0.0;;
	for(i=0;i<Np*Np;i++){
		if(fabs(ue[i]-V2[i])> e)	e = fabs(ue[i]-V2[i]);
	}

	printf("\nerr = %g\n", e);
        printf("r = %g\n", r);

	Save_Result();
	Free();

	printf("complete\n\n");
	return 0;
}

void Init(){

	int i,j;

	for(i=0;i<Np;i++){
		x[i] = a_c + i*h;
		y[i] = c_c + i*h;
	}

	for(i=0;i<Np;i++){
		for(j=0;j<Np;j++){
			V2[i*Np+j] = exp(-2*0)*cos(M_PI*x[i])*cos(M_PI*y[j]);
		}
	}

	fsource(V2, 0, x, y, F);
}
