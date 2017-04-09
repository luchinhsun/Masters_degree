/*
	creat by Chinhsun.Lu
*/

#include <stdio.h> 
#include <malloc.h>
#include <math.h>

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

int k_t = (tfinal/dt);

#define r (dt/h1)

float *Ax, *Ay, *V1, *V2, *W, *F, *V_tmp, *W_tmp;
float *b, *x, *y, *ut, *Vt;

void Allocate();
void Init();
void tam(float *A, float *d, float *x, int n);
void Save_Result();
void Free();

int main(){
	
	Allocate();
	Init();
	
	int i, j, k;
	
	//---------- Big loop for time t --------------------------------------
	for(k=0;k<k_t;k++){
		
		//--- sweep in x-direction --------------------------------------
		for(j=0;j<Np;j++){
			/*
			for(i=0;i<Np;i++){
				b[i] = 0.0;
			}
			*/
			for(i=0;i<Np;i++){
				if(j==0){
					b[i] = V1[i*Np+j] + r*(-V1[i*Np+j] + V1[i*Np+j+1]);
				}else if(j==Np-1){
					b[i] = V1[i*Np+j] + r*(V1[i*Np+j-1] - V1[i*Np+j]);
				}else{
					b[i] = V1[i*Np+j] + (r/2)*(V1[i*Np+j-1] - 2*V1[i*Np+j] + V1[i*Np+j+1]);
				}
			}
			tam(Ax, b, Vt, Np);
			for (i=0;i<Np;i++){
				V2[i*Np+j] = Vt[i];
			}
		}
		
		//-------- RK2 for ODE ---------------------------
		for(i=0;i<Np;i++){
			for(j=0;j<Np;j++){
				V_tmp[i*Np+j] = V2[i*Np+j] + (1/2)*dt*F[i*Np+j];
				W_tmp[i*Np+j] = W[i*Np+j] + (1/2)*dt*(-2*V2[i*Np+j]);
			}
		}
		for(i=0;i<Np;i++){
			for(j=0;j<Np;j++){
				V2[i*Np+j] = V2[i*Np+j] + dt*F[i*Np+j];
				W[i*Np+j] = W[i*Np+j] + dt*(-2*V_tmp[i*Np+j]);
			}
		}
		
		//-------------- loop in y -direction --------------------------------
		for(i=0;i<Np;i++){
			/*
			for(i=0;i<Np;i++){
				b[i] = 0.0;
			}
			*/
			for(j=0;j<Np;j++){
				if(i==0){
					b[j] = V2[i*Np+j] + (r/2)*(-2*V2[i*Np+j] + 2*V2[(i+1)*Np+j]);
				}else if(i==Np-1){
					b[j] = V2[i*Np+j] + (r/2)*(2*V2[(i-1)*Np+j] - 2*V2[i*Np+j]);
				}else{
					b[j] = V2[i*Np+j] + (r/2)*(V2[(i-1)*Np+j] - 2*V2[i*Np+j] + V2[(i+1)*Np+j]);
				}
			}
			tam(Ay, b, ut, Np);
			for (j=0;j<Np;j++){
				V1[i*Np+j] = ut[j];
			}
		}
	}
	
	Save_Result();
	Free();
	
	printf("\ncomplete\n\n");
	return 0;
}

void Allocate(){
	
	size_t size;
	size = Np*Np*sizeof(float);
	
	Ax = (float*)malloc(size);
	Ay = (float*)malloc(size);
	V1 = (float*)malloc(size);
	V2 = (float*)malloc(size);
	W = (float*)malloc(size);
	F = (float*)malloc(size);
	V_tmp = (float*)malloc(size);
	W_tmp = (float*)malloc(size);
	
	size = Np*sizeof(float);
	
	b = (float*)malloc(size);
	x = (float*)malloc(size);
	y = (float*)malloc(size);
	ut = (float*)malloc(size);
	Vt = (float*)malloc(size);
}

void Init(){
	
	int i,j;
	
	for(i=0;i<Np;i++){
		x[i] = a_c + i*h;
		y[i] = c_c + i*h;
	} 
	
	for(i=0;i<Np;i++){
		for(j=0;j<Np;j++){
			V1[i*Np+j] = exp(-2*M_PI*M_PI*0)*cos(M_PI*x[i])*cos(M_PI*y[j]);
			W[i*Np+j] = V1[i*Np+j];
		}
	}
	
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

void tam(float *A, float *d, float *x, int n){
	
	float *a, *b, *c;
	float *c_new, *d_new;
	
	size_t size;
	size = n*sizeof(float);
	
	a = (float*)malloc(size);
	b = (float*)malloc(size);
	c = (float*)malloc(size);
	c_new = (float*)malloc(size);
	d_new = (float*)malloc(size);
	
	a[0] = 0.0;
	c[n-1] = 0.0;

	int i;
	for(i=0;i<n-1;i++){
		c[i] = A[i*n+i+1];
		a[i+1] = A[(i+1)*n+i];
	}
	
	for(i=0;i<n;i++){
		b[i] = A[i*n+i];
	}
	
	c_new[0] = c[0]/b[0];
	for(i=1;i<n-1;i++){
		c_new[i] = c[i]/(b[i]-a[i]*c_new[i-1]);
	}
	
	d_new[0] = d[0]/b[0];
	for(i=1;i<n;i++){
		d_new[i] = (d[i]-a[i]*d_new[i-1])/(b[i]-a[i]*c_new[i-1]);	
	}
	
	x[n-1] = d_new[n-1];
	for(i=n-2;i>-1;i--){
		x[i] = d_new[i]-c_new[i]*x[i+1];
	}
	
	free(a);free(b);free(c);
	free(c_new);free(d_new);
}

void Save_Result(){

	FILE *pFile;
    int i,j;
    int index;
    int n;
    n = Np;
    pFile = fopen("V1.txt","w+");
    // Save the matrix A
    for (i = 0; i < n; i++) {
    	for (j = 0; j < n; j++) {
        	index = i*n + j;
            fprintf(pFile, "%g", V1[index]);
            if (j == (n-1)) {
        		fprintf(pFile, "\n");
            }else{
            	fprintf(pFile, "\t");
            }
    	}
	}
    fclose(pFile);
    /*
    pFile = fopen("Ax.txt","w+");
    // Save the matrix A
    for (i = 0; i < n; i++) {
    	for (j = 0; j < n; j++) {
        	index = i*n + j;
            fprintf(pFile, "%g", Ax[index]);
            if (j == (n-1)) {
        		fprintf(pFile, "\n");
            }else{
            	fprintf(pFile, "\t");
            }
    	}
	}
    fclose(pFile);
    pFile = fopen("Ay.txt","w+");
    // Save the matrix A
    for (i = 0; i < n; i++) {
    	for (j = 0; j < n; j++) {
        	index = i*n + j;
            fprintf(pFile, "%g", Ay[index]);
            if (j == (n-1)) {
        		fprintf(pFile, "\n");
            }else{
            	fprintf(pFile, "\t");
            }
    	}
	}
    fclose(pFile);
    pFile = fopen("b.txt","w+");
    for (i = 0; i < n; i++) {
        fprintf(pFile, "%g", b[i]);
        fprintf(pFile, "\t");
    }
    fclose(pFile);
    pFile = fopen("x.txt","w+");
    for (i = 0; i < n; i++) {
        fprintf(pFile, "%g", x[i]);
        fprintf(pFile, "\t");
    }
    fclose(pFile);
	*/
}

void Free(){
	free(Ax);free(Ay);free(V1);free(V2);
	free(W);free(F);free(V_tmp);free(W_tmp);
	free(b);free(x);free(y);free(ut);free(Vt);
}
