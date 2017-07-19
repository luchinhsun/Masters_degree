/*
	creat by Chinhsun.Lu
*/

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <sys/timeb.h>
#include <cusparse.h>

#define a_c 0.0
#define b_c 2.0
#define c_c 0.0
#define d_c 2.0
#define tfinal 0.5
#define N 2000//0~511
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

#define nnz (3*(Np-2)+4)

void Allocate();
void Init();
void fsource(float *V1,float t, float *x, float *y, float *F);
//void tam(float *A, float *d, float *x, int n);
//void tam();
//void tam2();
void RK();
void Save_Result();
void Free();

void Send_to_Device();
void Send_to_Host();

void Send_to_D();
void Send_to_H();

void sparse_Allocate_Memory();
void sparse_Send_To_Device();
//void ADI1(int j);
//void ADI1_2(int j);
//void ADI2(int i);
//void ADI2_2(int i);
void ADI1();
void ADI2();
void sparse_Free_Memory();
void cusparse_analysis();
