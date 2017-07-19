All: A B C D
	g++ ADI_RK.o memory.o CG.o fsource.o -o test.run
A:
	g++ memory.cpp -c
B:
	g++ CG.cpp -c
C:
	g++ ADI_RK.cpp -c
D:
	g++ fsource.cpp -c
clean:
	rm *.o *.txt *.run
