all: Gauss

Gauss: main.o ErrorNorm.o MatrixNorm.o MatrixInput.o MatrixOutput.o Function.o Gauss.o
	g++ main.o ErrorNorm.o MatrixNorm.o MatrixInput.o MatrixOutput.o Function.o Gauss.o -lpthread -o Gauss

main.o: main.cpp
	g++ -c main.cpp -lpthread -O3

ErrorNorm.o: ErrorNorm.cpp
	g++ -c ErrorNorm.cpp -lpthread -O3

MatrixNorm.o: MatrixNorm.cpp
	g++ -c MatrixNorm.cpp -lpthread -O3

MatrixInput.o: MatrixInput.cpp
	g++ -c MatrixInput.cpp -lpthread -O3

MatrixOutput.o: MatrixOutput.cpp
	g++ -c MatrixOutput.cpp -lpthread -O3

Function.o: Function.cpp
	g++ -c Function.cpp -lpthread -O3

Gauss.o: Gauss.cpp
	g++ -c Gauss.cpp -lpthread-O3

clean: 
	rm -rf *.o Gauss
