CUDA_LIB :=/usr/local/cuda/lib64 -lcuda -lcudart

All: A B C D
	g++ euler.o memory.o forward.o fsource.o -o test.run -L $(CUDA_LIB)
A:
	nvcc -arch=sm_61 --ptxas-options=-v memory.cu -c
B:
	nvcc -arch=sm_61 --ptxas-options=-v forward.cu -c
C:
	g++ euler.cpp -c
D:
	g++ fsource.cpp -c
clean:
	rm *.o *.txt *.run
