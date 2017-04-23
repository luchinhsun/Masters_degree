#include "head.h"

int k_t = (tfinal/dt);

float *Ax, *Ay, *V1, *V2, *W, *F, *V_tmp, *W_tmp, *ue;
float *b, *x, *y, *ut, *Vt;
float *a_tam, *b_tam, *c_tam;
float *c_new_tam, *d_new_tam;

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
	ue = (float*)malloc(size);

	size = Np*sizeof(float);

	b = (float*)malloc(size);
	x = (float*)malloc(size);
	y = (float*)malloc(size);
	ut = (float*)malloc(size);
	Vt = (float*)malloc(size);

	size = Np*sizeof(float);

	a_tam = (float*)malloc(size);
        b_tam = (float*)malloc(size);
        c_tam = (float*)malloc(size);
        c_new_tam = (float*)malloc(size);
        d_new_tam = (float*)malloc(size);

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
    pFile = fopen("F.txt","w+");
    // Save the matrix A
    for (i = 0; i < n; i++) {
    	for (j = 0; j < n; j++) {
        	index = i*n + j;
            fprintf(pFile, "%g", F[index]);
            if (j == (n-1)) {
        		fprintf(pFile, "\n");
            }else{
            	fprintf(pFile, "\t");
            }
    	}
	}
    fclose(pFile);
    pFile = fopen("V_tmp.txt","w+");
    // Save the matrix A
    for (i = 0; i < n; i++) {
    	for (j = 0; j < n; j++) {
        	index = i*n + j;
            fprintf(pFile, "%g", V_tmp[index]);
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
	free(W);free(F);free(V_tmp);free(W_tmp);free(ue);
	free(b);free(x);free(y);free(ut);free(Vt);

	free(a_tam);free(b_tam);free(c_tam);
        free(c_new_tam);free(d_new_tam);
}
