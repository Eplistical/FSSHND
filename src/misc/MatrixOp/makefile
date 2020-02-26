CXX = g++ 
OPT = -O3 -fPIC
SRC = .

all : matmat.o eig.o

matmat.o : matmat.cpp
	$(CXX) $(OPT) -c $^ -o $@

eig.o : eig.cpp
	$(CXX) $(OPT) -c $^ -o $@

clean: 
	rm *.o
