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

//int k_t = (tfinal/dt);

#define r (dt/h1)

//float *Ax, *Ay, *V1, *V2, *W, *F, *V_tmp, *W_tmp;
//float *b, *x, *y, *ut, *Vt;

void Allocate();
void Init();
void fsource(float *V1,float t, float *x, float *y, float *F);
void tam(float *A, float *d, float *x, int n);
void Save_Result();
void Free();
