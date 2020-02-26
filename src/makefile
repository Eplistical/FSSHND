CXX = g++ 
MPICXX = mpicxx
OPT = -O0

LIBS += -L./misc -Wl,--rpath=./misc -lmisc
LIBS += -lboost_program_options -lopenblas -llapack -lpthread -lgfortran 
SRC = .

all : main

main : main.cpp hamiltonian.o tully1_hamiltonian.o
	$(CXX) -o $@ $^ $(LIBS)

hamiltonian.o : $(SRC)/hamiltonian.cpp
	$(CXX) $(OPT) -c $^ -o $@

tully1_hamiltonian.o : $(SRC)/tully1_hamiltonian.cpp
	$(CXX) $(OPT) -c $^ -o $@

misc : 
	$(MAKE) -C $@

clean: 
	rm *.o

.PHONY: all clean misc
