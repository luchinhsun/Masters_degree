/*
	creat by Chinhsun.Lu
*/

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <sys/timeb.h>

#define a_c 0.0
#define b_c 2.0
#define c_c 0.0
#define d_c 2.0
#define tfinal 0.5
#define N ((32-2)*7+32-3)
#define M N
#define Mp (M+1)
#define Np (N+1)
#define h ((b_c-a_c)/N)
#define h1 (h*h)
#define dt 1e-5
#define r (dt/h1)

void Allocate();
void Init();
void fsource(float *V2,float t, float *x, float *y, float *F);
void Save_Result();
void Free();

void Send_to_Device();
void Send_to_Host();
void forward_euler();
