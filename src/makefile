ROOTDIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

CXX = g++ 
MPICXX = mpicxx
OPT += -std=c++14
#OPT += -O0
OPT += -O3 

LIBS += -L$(ROOTDIR)/misc -Wl,--rpath=$(ROOTDIR)/misc -lmisc
LIBS += -lboost_program_options -lopenblas -llapack -lpthread -lgfortran 


all : tully1 tully1_2d multistate nuscatter

tully1 : hamiltonian.o tully1_hamiltonian.o tully1.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

tully1_2d : hamiltonian.o tully1_2d_hamiltonian.o tully1_2d.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

multistate : hamiltonian.o multistate_hamiltonian.o multistate.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

nuscatter : hamiltonian.o nuscatter_hamiltonian.o nuscatter.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

hamiltonian.o : hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

tully1_hamiltonian.o : tully1_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

tully1_2d_hamiltonian.o : tully1_2d_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

multistate_hamiltonian.o : multistate_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

nuscatter_hamiltonian.o : nuscatter_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

misc : 
	$(MAKE) -C $@

clean: 
	rm *.o

.PHONY: all clean misc
