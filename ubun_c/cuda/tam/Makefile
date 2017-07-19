CUDA_LIB :=/usr/local/cuda/lib64 -lcuda -lcudart

All: A B C D E
	g++ ADI_RK.o memory.o tam.o fsource.o RK.o -o test.run -L $(CUDA_LIB)
A:
	nvcc --ptxas-options=-v -arch=sm_61 memory.cu -c
B:
	nvcc --ptxas-options=-v -arch=sm_61 tam.cu -c
C:
	g++ ADI_RK.cpp -c
D:
	g++ fsource.cpp -c
E:
	nvcc --ptxas-options=-v -arch=sm_61 RK.cu -c
clean:
	rm *.o *.txt *.run
