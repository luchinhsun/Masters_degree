// CG (Conjugate Gradient) Demonstration code
// Prof. Matthew Smith, NCKU, 2016
// All together in one file to make
// it easier to build

// This code will solve a steady 2D heat transfer problem.
// Here we use 2D central differences to create an Ax = B
// matrix. The init() function creates A and B - A is a 
// NxN matrix while B is a 1xN vector. Our solution is x.
/*
#include <stdio.h>
#include <malloc.h>

#define NX 25   // Solve a 25x25 matrix
#define NY 25
#define N (NX*NY)			// Number of cells - our matrix is N*N	
#define L 1.0
#define H 1.0
#define DX (L/NX)
#define DY (H/NY)
#define K 0.1				// Conductivity
#define PHI_X (K/(DX*DX))
#define PHI_Y (K/(DY*DY))	
//#define float double

// Define our variables of interest
float *h_A;	// A matrix
float *h_B;	// B vector
*/
#include "head.h"

//float *h_x;	// x (solution) vector
float *h_R;	// Residual vector
float *h_P;
float *h_AP;
float *h_scalars;  // Scalars


void Allocate_Memory() {

	size_t size;

	// Our N*N variable (A)
	//size = N*N*sizeof(float);
	//h_A = (float*)malloc(size); 
	
	// Our 1D (N) variables
	size = Np*sizeof(float);
	//h_B = (float*)malloc(size); 
	//h_x = (float*)malloc(size);
	h_R = (float*)malloc(size);
	h_P = (float*)malloc(size);
	h_AP = (float*)malloc(size);
	
	// Small array holding scalars
	size = 5*sizeof(float);
	h_scalars = (float*)malloc(size);

}


void Free_Memory() {

	// Now we better free the memory on the GPU
	//free(h_A); 
	//free(h_B); 
	//free(h_x); 
	free(h_R); 
	free(h_P);
	free(h_AP);
	free(h_scalars);

}

void SetUp_CG(float *B, float *R, float *P, int n) {

	for (int i = 0; i < n; i++) {
		R[i] = B[i];
		P[i] = B[i];
	}

}


float DotProduct(float *x, float *y, int n) {

	float temp = 0;
	for (int i = 0; i < n; i++) {
		temp = temp + x[i]*y[i];
	}
	return temp; 
}


void MatrixVectorProduct(float *x, float *y, float *z, int n) {

	int i ,j; 
	int index;
	
	for (i = 0; i < n; i++) {
		z[i] = 0.0;
		index = i*n; // Jump to the start
		for (j = 0; j < n; j++) {
			z[i] += x[index+j]*y[j];
		}
	}
}


void Update_x(float *x, float *P, float alpha, int n) {

	for (int i = 0; i < n; i++) {
		x[i] = x[i] + alpha*P[i];
	}

}

void Update_P(float *P, float *R, float beta, int n) {

	for (int i = 0; i < n; i++) {
		P[i] = R[i] + beta*P[i];
	}

}

void Update_R(float *R, float *AP, float alpha, int n) {

	for (int i = 0; i < n; i++) {
		R[i] = R[i] - alpha*AP[i];
	}

}

/*
void Init() {

	int i, j;
	int index;
	int x_cell, y_cell;
	// Set up our A and B matrix

	// Full 2D Heat Transfer test
	for (i = 0; i < N; i++) {

		y_cell = (int)(i/NX);
		x_cell = i - y_cell*NX;

		// Find the diagonal		
		index = i*N + i;
		h_A[index] = 2.0*PHI_X + 2.0*PHI_Y;

		// Modify A for left and right				
		if (x_cell > 0) {
			h_A[index-1] = -PHI_X;
		}

		if (x_cell < (NX-1)) {
			h_A[index+1] = -PHI_X;
		}

		// Modifiy for up and down
		if (y_cell > 0) {
			h_A[index-NX] = -PHI_Y;
		}
		if (y_cell < (NY-1)) {
			h_A[index+NX] = -PHI_Y;
		}


		// Set B now
		if (y_cell == 0) {
			h_B[i] = PHI_Y;
		} else {
			h_B[i] = 0.0;
		}

		// And our initial guess x
		h_x[i] = 0.0; 
	}
		

}


void Save_Result() {

	FILE *pFile;
	int i,j;
	int index;
	pFile = fopen("A.txt","w");
	// Save the matrix A
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			index = i*N + j;
			fprintf(pFile, "%g", h_A[index]);
			if (j == (N-1)) {
				fprintf(pFile, "\n");
			} else {
				fprintf(pFile, "\t");
			}
		}
	}
	fclose(pFile);

	pFile = fopen("B.txt","w");
	// Save the vector B
	for (i = 0; i < N; i++) {
		fprintf(pFile, "%g\n", h_B[i]);
	}
	fclose(pFile);

	pFile = fopen("X.txt","w");
	// Save the vector B
	for (i = 0; i < N; i++) {
		fprintf(pFile, "%g\n", h_x[i]);
	}
	fclose(pFile);

	pFile = fopen("R.txt","w");
	// Save the vector R
	for (i = 0; i < N; i++) {
		fprintf(pFile, "%g\n", h_R[i]);
	}
	fclose(pFile);


}
*/
void CG(float *h_A, float *h_x, float *h_B, int n) {
	
	int i;
	float alpha;
	float beta;
	float PTAP, RTR, RTR_NEW;
	// Allocate the memory
	Allocate_Memory();

	//Init();

	SetUp_CG(h_B, h_R, h_P, n);	

    // Here I will fix the number of iterations to 33. 
    // Theoretically we need 625 (worst case) but practically speaking
    // we know CG will converge with sqrt(N) ~ 25 iterations.
	for (i = 0; i < 33; i++) {

		// Compute AP
		MatrixVectorProduct(h_A, h_P, h_AP, n); // Send it to R so we can get it back and have a look
		
		// Now to compute PTAP 
		PTAP = DotProduct(h_P, h_AP, n);   

		// Compute RTR
		RTR = DotProduct(h_R, h_R, n); 
		
		// Update X (and compute alpha)
		alpha = RTR/PTAP;
		Update_x(h_x, h_P, alpha, n);
		Update_R(h_R, h_AP, alpha, n);

		// Compute the new residual
		RTR_NEW = DotProduct(h_R, h_R, n);

		// Update P (and compute beta)
		beta = RTR_NEW/RTR;
		Update_P(h_P, h_R, beta, n);
		
		//printf("Iteration %i = RTR = %g, RTR_new = %g, PTAP = %g, ALPHA = %g, BETA = %g\n", i, RTR, RTR_NEW, PTAP, alpha, beta);
		if(RTR_NEW<1e-4)        break;
	}	

	// Save result
	//Save_Result();

	// Free memory
	//Free_Memory();

	//printf("Complete\n");

}







