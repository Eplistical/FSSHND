ROOTDIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

CXX = g++ 
MPICXX = mpicxx
OPT += -std=c++14
#OPT += -O0
OPT += -O3 

LIBS += -L$(ROOTDIR)/misc -Wl,--rpath=$(ROOTDIR)/misc -lmisc
LIBS += -lboost_program_options -lopenblas -llapack -lpthread -lgfortran


all : tully1 tully1_2d flat_2d multistate conner nuscatter three_state three_state_1d test1d2s yanze_straight

tully1 : hamiltonian.o adjust_evt_phase.o tully1_hamiltonian.o tully1.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

tully1_2d : hamiltonian.o adjust_evt_phase.o tully1_2d_hamiltonian.o tully1_2d.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

flat_2d : hamiltonian.o adjust_evt_phase.o flat_2d_hamiltonian.o flat_2d.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

multistate : hamiltonian.o adjust_evt_phase.o multistate_hamiltonian.o multistate.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

three_state : hamiltonian.o adjust_evt_phase.o three_state_hamiltonian.o three_state.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

three_state_1d : hamiltonian.o adjust_evt_phase.o three_state_1d_hamiltonian.o three_state_1d.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

conner : hamiltonian.o adjust_evt_phase.o conner_hamiltonian.o conner.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

nuscatter : hamiltonian.o adjust_evt_phase.o nuscatter_hamiltonian.o nuscatter.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

test1d2s : hamiltonian.o adjust_evt_phase.o test1d2s_hamiltonian.o test1d2s.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

yanze_straight : hamiltonian.o adjust_evt_phase.o yanze_straight_hamiltonian.o yanze_straight.cpp
	$(MPICXX) $(OPT) -o $@ $^ $(LIBS)

hamiltonian.o : hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

adjust_evt_phase.o : adjust_evt_phase.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

tully1_hamiltonian.o : tully1_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

tully1_2d_hamiltonian.o : tully1_2d_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

flat_2d_hamiltonian.o : flat_2d_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

multistate_hamiltonian.o : multistate_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

three_state_hamiltonian.o : three_state_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

three_state_1d_hamiltonian.o : three_state_1d_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

conner_hamiltonian.o : conner_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

nuscatter_hamiltonian.o : nuscatter_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

test1d2s_hamiltonian.o : test1d2s_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

yanze_straight_hamiltonian.o : yanze_straight_hamiltonian.cpp
	$(CXX) $(OPT) -o $@ -c $^ 

misc : 
	$(MAKE) -C $@

clean: 
	rm *.o

.PHONY: all clean misc
