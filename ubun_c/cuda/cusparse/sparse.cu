#include "head.h"
//#include <time.h>
/*
clock_t start_cu;
clock_t end_cu;
float time_used_cu;
*/
//variable for cusparse
cusparseStatus_t status;
cusparseHandle_t handle=0;
cusparseMatDescr_t descr=0;
cusparseMatDescr_t descrL=0;
cusparseMatDescr_t descrU=0;
cusparseSolveAnalysisInfo_t infoA=0;
cusparseSolveAnalysisInfo_t info_u=0;
int *cooRowIndexHostPtr;
int * cooColIndexHostPtr;
float * cooValHostPtr;
int *cooRowIndex;
int * cooColIndex;
float * cooVal;
float * cooValLU;
float * temp;
float * d_V1_t;
int * csrRowPtr;

float done =1.0;

extern float *d_b, *d_ut, *d_Vt;
extern float *d_V1, *d_V2;

void sparse_Allocate_Memory(){
	//cusparse
	size_t size = nnz*sizeof(int);
	cooRowIndexHostPtr = (int *) malloc(size);
	cooColIndexHostPtr = (int *) malloc(size);
	cooValHostPtr = (float *)malloc(nnz*sizeof(float));

	cooRowIndexHostPtr[0] = 0;cooColIndexHostPtr[0]=0;cooValHostPtr[0]=r+1;
	cooRowIndexHostPtr[1] = 0;cooColIndexHostPtr[1]=1;cooValHostPtr[1]=-r;

	cooRowIndexHostPtr[2] = 1;cooColIndexHostPtr[2]=0;cooValHostPtr[2]=-r/2;
	cooRowIndexHostPtr[3] = 1;cooColIndexHostPtr[3]=1;cooValHostPtr[3]=r+1;
	cooRowIndexHostPtr[4] = 1;cooColIndexHostPtr[4]=2;cooValHostPtr[4]=-r/2;
	int i;
	for(i=5;i<(nnz-3);i=i+3){
		cooRowIndexHostPtr[i] = cooRowIndexHostPtr[i-3]+1;	
		cooColIndexHostPtr[i] = cooColIndexHostPtr[i-3]+1;
		cooRowIndexHostPtr[i+1] = cooRowIndexHostPtr[i];	
		cooColIndexHostPtr[i+1] = cooColIndexHostPtr[i]+1;
		cooRowIndexHostPtr[i+2] = cooRowIndexHostPtr[i+1];	
		cooColIndexHostPtr[i+2] = cooColIndexHostPtr[i+1]+1;
		cooValHostPtr[i]=-r/2;
		cooValHostPtr[i+1]=r+1;
		cooValHostPtr[i+2]=-r/2;
	}
	cooRowIndexHostPtr[nnz-2] = Np-1;cooColIndexHostPtr[nnz-2]=Np-2;cooValHostPtr[nnz-2]=-r;
        cooRowIndexHostPtr[nnz-1] = Np-1;cooColIndexHostPtr[nnz-1]=Np-1;cooValHostPtr[nnz-1]=r+1;

	cudaError_t Error;

	Error = cudaMalloc((void**)&cooRowIndex, size);
	printf("CUDA error(malloc RowIndex) = %s\n",cudaGetErrorString(Error));
	Error = cudaMalloc((void**)&cooColIndex, size);
	printf("CUDA error(malloc ColIndex) = %s\n",cudaGetErrorString(Error));
	Error = cudaMalloc((void**)&cooVal, nnz*sizeof(float));
	printf("CUDA error(malloc Val) = %s\n",cudaGetErrorString(Error));
	Error = cudaMalloc((void**)&cooValLU, nnz*sizeof(float));
        printf("CUDA error(malloc Val) = %s\n",cudaGetErrorString(Error));

	Error = cudaMalloc((void**)&temp, Np*Np*sizeof(float));
        printf("CUDA error(malloc temp) = %s\n",cudaGetErrorString(Error));
	Error = cudaMalloc((void**)&d_V1_t, Np*Np*sizeof(float));
        printf("CUDA error(malloc d_V1_t) = %s\n",cudaGetErrorString(Error));
	Error = cudaMalloc((void**)&csrRowPtr,(Np+1)*sizeof(int));
        printf("CUDA error(malloc csrRowPtr) = %s\n",cudaGetErrorString(Error));

	status= cusparseCreate(&handle);
	status= cusparseCreateMatDescr(&descr);

	cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

	status = cusparseCreateMatDescr(&descrL);
	cusparseSetMatType(descrL,CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrL,CUSPARSE_INDEX_BASE_ZERO);
	cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
    	cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT);

    	status = cusparseCreateMatDescr(&descrU);
    	cusparseSetMatType(descrU,CUSPARSE_MATRIX_TYPE_GENERAL);
    	cusparseSetMatIndexBase(descrU,CUSPARSE_INDEX_BASE_ZERO);
    	cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
    	cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);

        status = cusparseCreateSolveAnalysisInfo(&infoA);
        status = cusparseCreateSolveAnalysisInfo(&info_u);
}

void sparse_Send_To_Device(){
	cudaError_t Error;
	size_t size = nnz*sizeof(int);
	Error = cudaMemcpy(cooRowIndex, cooRowIndexHostPtr, size, cudaMemcpyHostToDevice);
	printf("CUDA error(memcpy RowIndex) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(cooColIndex, cooColIndexHostPtr, size, cudaMemcpyHostToDevice);
	printf("CUDA error(memcpy ColIndex) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(cooVal, cooValHostPtr, (size_t)(nnz*sizeof(float)), cudaMemcpyHostToDevice);
	printf("CUDA error(memcpy Val) = %s\n",cudaGetErrorString(Error));
}

void cusparse_analysis(){

	status= cusparseXcoo2csr(handle,cooRowIndex,nnz,Np, csrRowPtr,CUSPARSE_INDEX_BASE_ZERO);
	/*
        if (status != CUSPARSE_STATUS_SUCCESS) {
                printf("shit1");
        }
	*/
	status= cusparseScsrsm_analysis(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, Np, nnz, descr,
                                                cooVal, csrRowPtr, cooColIndex, infoA);
        cudaMemcpy(cooValLU, cooVal, nnz*sizeof(float), cudaMemcpyDeviceToDevice);
        status = cusparseScsrilu0(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, Np, descr,
                                                cooValLU, csrRowPtr, cooColIndex, infoA);
        status = cusparseScsrsm_analysis(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, Np, nnz, descrU,
                                                cooVal, csrRowPtr, cooColIndex, info_u);
}

void sparse(){
	status = cusparseScsrsm_solve(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, Np, Np, &done, descrL,
                                              cooValLU, csrRowPtr, cooColIndex, infoA, d_V2, Np,
						temp, Np);
	status = cusparseScsrsm_solve(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, Np, Np, &done, descrU,
                                              cooValLU, csrRowPtr, cooColIndex, info_u, temp, Np,
						d_V2, Np);
	/*
	if (status != CUSPARSE_STATUS_SUCCESS) {
                printf("shit2");
        }
	*/
}

void sparse2(){
        status = cusparseScsrsm_solve(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, Np, Np, &done, descrL,
                                              cooValLU, csrRowPtr, cooColIndex, infoA, d_V1_t, Np,
						temp, Np);
        status = cusparseScsrsm_solve(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, Np, Np, &done, descrU,
                                              cooValLU, csrRowPtr, cooColIndex, info_u, temp, Np,
						d_V1_t, Np);
	/*
        if (status != CUSPARSE_STATUS_SUCCESS) {
                printf("shit4");
        }
	*/
}

void sparse_Free_Memory(){

	status = cusparseDestroyMatDescr(descr); descr = 0;
	status = cusparseDestroy(handle); handle = 0;
	status = cusparseDestroyMatDescr(descrL); descrL = 0;
        status = cusparseDestroyMatDescr(descrU); descrU = 0;

        if (cooRowIndexHostPtr) free(cooRowIndexHostPtr);
        if (cooColIndexHostPtr) free(cooColIndexHostPtr);
        if (cooValHostPtr) free(cooValHostPtr);
	if (temp) cudaFree(temp);
	if (d_V1_t) cudaFree(d_V1_t);
        if (csrRowPtr) cudaFree(csrRowPtr);
        if (cooRowIndex) cudaFree(cooRowIndex);
        if (cooColIndex) cudaFree(cooColIndex);
        if (cooVal) cudaFree(cooVal);
	if (cooValLU) cudaFree(cooValLU);
        if (descr) cusparseDestroyMatDescr(descr);
	if (descrL) cusparseDestroyMatDescr(descrL);
	if (descrU) cusparseDestroyMatDescr(descrU);
        if (handle) cusparseDestroy(handle);


	cusparseDestroySolveAnalysisInfo(infoA);
        cusparseDestroySolveAnalysisInfo(info_u);

}

__global__ void GPU_adi_x(float *d_V1, float *d_V2){
        int i = blockDim.x * blockIdx.x + threadIdx.x;

	int j;
        if(i<Np){
	for(j=0;j<Np;++j){
        if(j==0){
                d_V2[i*Np+j] = d_V1[i*Np+j] + r*(-d_V1[i*Np+j] + d_V1[i*Np+j+1]);
        }else if(j==Np-1){
                d_V2[i*Np+j] = d_V1[i*Np+j] + r*(d_V1[i*Np+j-1] - d_V1[i*Np+j]);
        }else{
                d_V2[i*Np+j] = d_V1[i*Np+j] + (r/2)*(d_V1[i*Np+j-1] - 2*d_V1[i*Np+j] + d_V1[i*Np+j+1]);
        }
        }
	}

}

__global__ void GPU_adi_y(float *d_V2, float *d_V1){
        int j = blockDim.x * blockIdx.x + threadIdx.x;

	int i;
        if(j<Np){
	for(i=0;i<Np;++i){
        if(i==0){
                d_V1[j*Np+i] = d_V2[i*Np+j] + (r/2)*(-2*d_V2[i*Np+j] + 2*d_V2[(i+1)*Np+j]);
        }else if(i==Np-1){
                d_V1[j*Np+i] = d_V2[i*Np+j] + (r/2)*(2*d_V2[(i-1)*Np+j] - 2*d_V2[i*Np+j]);
        }else{
                d_V1[j*Np+i] = d_V2[i*Np+j] + (r/2)*(d_V2[(i-1)*Np+j] - 2*d_V2[i*Np+j] + d_V2[(i+1)*Np+j]);
        }
	}
        }
}

__global__ void GPU_trans(float *d_V1, float *d_V1_t){
        int k = blockDim.x * blockIdx.x + threadIdx.x;
        int i;
        int j;

        if(k<Np*Np){
                i = k/Np;
                j = k%Np;
                d_V1[i*Np+j] = d_V1_t[j*Np+i];
        }
}

void ADI1(){
        int tpb = 256;
        int bpg = (Np+tpb-1)/tpb;

        GPU_adi_x<<<bpg, tpb>>>(d_V1, d_V2);
	//start_cu = clock();
        sparse();
	//end_cu = clock();
	//time_used_cu = (float)(end_cu - start_cu)/ CLOCKS_PER_SEC;
	//printf("\n time in sparse kernal = %f\n",time_used_cu);
}

void ADI2(){
        int tpb = 256;
        int bpg = (Np+tpb-1)/tpb;

        GPU_adi_y<<<bpg, tpb>>>(d_V2, d_V1_t);
        sparse2();
	bpg = (Np*Np+tpb-1)/tpb;
        GPU_trans<<<bpg, tpb>>>(d_V1, d_V1_t);
}
