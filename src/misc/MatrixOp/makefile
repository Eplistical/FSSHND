CXX = g++ 
OPT = -O3 -fPIC
SRC = .

all : matmat.o eig.o logm.o determinant.o

matmat.o : matmat.cpp
	$(CXX) $(OPT) -c $^ -o $@

eig.o : eig.cpp
	$(CXX) $(OPT) -c $^ -o $@

logm.o : logm.cpp
	$(CXX) $(OPT) -c $^ -o $@

determinant.o : determinant.cpp
	$(CXX) $(OPT) -c $^ -o $@

clean: 
	rm *.o
