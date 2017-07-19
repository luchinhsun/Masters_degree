/*
	creat by Chinhsun.Lu
*/

#include "head.h"
/*
#define a_c 0.0
#define b_c 2.0
#define c_c 0.0
#define d_c 2.0
#define tfinal 0.5
#define N 80
#define M N
#define Mp (M+1)
#define Np (N+1)
#define h ((b_c-a_c)/N)
#define h1 (h*h)
#define dt 0.01
*/
extern int k_t;

//#define r (dt/h1)

extern float *Ax, *Ay, *V1, *V2, *W, *F, *V_tmp, *W_tmp, *ue;
extern float *b, *x, *y, *ut, *Vt;

int main(){

	Allocate();
	Init();

	int i, j, k;

	float t = 0.0;

	struct timeb start, end;
	int diff;

	ftime(&start);
	//---------- Big loop for time t --------------------------------------
	for(k=0;k<k_t;k++){

		//--- sweep in x-direction --------------------------------------
		for(j=0;j<Np;j++){
			for(i=0;i<Np;i++){
				if(j==0){
					b[i] = V1[i*Np+j] + r*(-V1[i*Np+j] + V1[i*Np+j+1]);
				}else if(j==Np-1){
					b[i] = V1[i*Np+j] + r*(V1[i*Np+j-1] - V1[i*Np+j]);
				}else{
					b[i] = V1[i*Np+j] + (r/2)*(V1[i*Np+j-1] - 2*V1[i*Np+j] + V1[i*Np+j+1]);
				}
			}
			//tam(Ax, b, Vt, Np);
			CG(Ax, Vt, b, Np);
			for (i=0;i<Np;i++){
				V2[i*Np+j] = Vt[i];
			}
		}

		//-------- RK2 for ODE ---------------------------
		for(i=0;i<Np;i++){
			for(j=0;j<Np;j++){
				V_tmp[i*Np+j] = V2[i*Np+j] + (1.0/2.0)*dt*F[i*Np+j];
				//W_tmp[i*Np+j] = W[i*Np+j] + (1/2)*dt*(-2*V2[i*Np+j]);
			}
		}

		t = t+dt/2;

		fsource(V_tmp, t, x, y, F);

		for(i=0;i<Np;i++){
			for(j=0;j<Np;j++){
				V2[i*Np+j] = V2[i*Np+j] + dt*F[i*Np+j];
				//W[i*Np+j] = W[i*Np+j] + dt*(-2*V_tmp[i*Np+j]);
			}
		}

		t = t+dt/2;

		//-------------- loop in y -direction --------------------------------
		for(i=0;i<Np;i++){
			for(j=0;j<Np;j++){
				if(i==0){
					b[j] = V2[i*Np+j] + (r/2)*(-2*V2[i*Np+j] + 2*V2[(i+1)*Np+j]);
				}else if(i==Np-1){
					b[j] = V2[i*Np+j] + (r/2)*(2*V2[(i-1)*Np+j] - 2*V2[i*Np+j]);
				}else{
					b[j] = V2[i*Np+j] + (r/2)*(V2[(i-1)*Np+j] - 2*V2[i*Np+j] + V2[(i+1)*Np+j]);
				}
			}
			//tam(Ay, b, ut, Np);
			CG(Ay, ut, b, Np);
			for (j=0;j<Np;j++){
				V1[i*Np+j] = ut[j];
			}
		}
	}
	ftime(&end);
	diff = (int)(1000.0*(end.time-start.time)+(end.millitm-start.millitm));
	printf("\nTime = %d ms\n", diff);

	//-------------- check error between exact solution  --------------------------------
	for(i=0;i<Np;i++){
		for(j=0;j<Np;j++){
			ue[i*Np+j] = exp(-2.0*tfinal)*cos(M_PI*x[i])*cos(M_PI*y[j]);
		}
	}

	float e;
	for(i=0;i<Np*Np;i++){
		if(fabs(ue[i]-V1[i])>0.0)	e = fabs(ue[i]-V1[i]);
	}

	printf("\nerr = %g\n", e);

	//Save_Result();
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
			V1[i*Np+j] = exp(-2*0)*cos(M_PI*x[i])*cos(M_PI*y[j]);
			//W[i*Np+j] = V1[i*Np+j];
		}
	}

	fsource(V1, 0, x, y, F);

	for(i=0;i<Np;i++){
			if(i == 0){
				Ax[i*Np+i+1] = -r;
				Ay[i*Np+i+1] = -r;
			}else if(i == Np-1){
				Ax[i*Np+i-1] = -r;
				Ay[i*Np+i-1] = -r;
			}else{
				Ax[i*Np+i+1] = -r/2;
				Ax[i*Np+i-1] = -r/2;
				Ay[i*Np+i+1] = -r/2;
				Ay[i*Np+i-1] = -r/2;
			}
			Ax[i*Np+i] = 1+r;
			Ay[i*Np+i] = 1+r;
	}
}
