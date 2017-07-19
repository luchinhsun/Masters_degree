CUDA_LIB :=/usr/local/cuda/lib64 -lcuda -lcudart -lcusparse_static -lculibos

CUSPARSE_LIB :=/usr/local/cuda/include

All: A B C D E
	g++ ADI_RK.o memory.o sparse.o fsource.o RK.o -o test.run -L $(CUDA_LIB) -lpthread -ldl
A:
	nvcc -arch=sm_61 memory.cu -c
B:
	nvcc -arch=sm_61 sparse.cu -c -I $(CUSPARSE_LIB)
C:
	g++ ADI_RK.cpp -c -I $(CUSPARSE_LIB)
D:
	g++ fsource.cpp -c -I $(CUSPARSE_LIB)
E:
	nvcc -arch=sm_61 RK.cu -c
clean:
	rm *.o *.txt *.run
