#include "head.h"

int k_t = (tfinal/dt);

float *V1, *V2, *F, *ue;
float *x, *y;
float *t;

void Allocate(){

	size_t size;
	size = (Np+2)*(Np+2)*sizeof(float);

	V1 = (float*)malloc(size);

	size = Np*Np*sizeof(float);

	V2 = (float*)malloc(size);
	F = (float*)malloc(size);
	ue = (float*)malloc(size);

	size = Np*sizeof(float);

	x = (float*)malloc(size);
	y = (float*)malloc(size);

	size = sizeof(float);
	t = (float*)malloc(size);
}

void Save_Result(){

	FILE *pFile;
	int i,j;
	int index;
	int n;
	n = Np;
	pFile = fopen("V1.txt","w+");
	// Save the matrix V1
	for (i = 0; i < n+2; i++) {
		for (j = 0; j < n+2; j++) {
			index = i*(n+2) + j;
			fprintf(pFile, "%g", V1[index]);
			if (j == n+1) {
				fprintf(pFile, "\n");
			}else{
				fprintf(pFile, "\t");
			}
		}
	}
	fclose(pFile);

	pFile = fopen("V2.txt","w+");
        // Save the matrix V2
        for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                        index = i*n + j;
                        fprintf(pFile, "%g", V2[index]);
                        if (j == n-1) {
                                fprintf(pFile, "\n");
                        }else{
                                fprintf(pFile, "\t");
                        }
                }
        }
        fclose(pFile);
}

void Free(){
	free(V1);free(V2);free(F);free(ue);
	free(x);free(y);
	free(t);
}
