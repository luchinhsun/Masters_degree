All: A B C D
	g++ ADI_RK.o memory.o tam.o fsource.o -fopenmp -o test.run
A:
	g++ memory.cpp -c
B:
	g++ tam.cpp -c -fopenmp
C:
	g++ ADI_RK.cpp -c -fopenmp
D:
	g++ fsource.cpp -c -fopenmp
clean:
	rm *.o *.txt *.run
