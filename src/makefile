CXX = g++ 
MPICXX = mpicxx
OPT = -O0

LIBS += -L./misc
LIBS += -lopenblas -llapack -lpthread -lgfortran 
LIBS += -lmisc
SRC = .

all : main

main : hamiltonian.o tully1_hamiltonian.o fssh_trajectory.o main.o
	$(CXX) -o $@ $^ $(LIBS)

fssh_trajectory.o : $(SRC)/fssh_trajectory.cpp
	$(CXX) $(OPT) -c $^ -o $@

hamiltonian.o : $(SRC)/hamiltonian.cpp
	$(CXX) $(OPT) -c $^ -o $@

tully1_hamiltonian.o : $(SRC)/tully1_hamiltonian.cpp
	$(CXX) $(OPT) -c $^ -o $@

misc : 
	$(MAKE) -C $@

clean: 
	rm *.o

.PHONY: all clean misc
