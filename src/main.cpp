#include <iostream>
#include <complex>
#include <vector>
#include "tully1_hamiltonian.hpp"
#include "fssh_trajectory.hpp"
#include "misc/matrixop.hpp"

using namespace mqc;
using namespace std;

using hamiltonian_t = Tully1_Hamiltonian;
using trajectory_t = FSSH_Trajectory<hamiltonian_t>;

int ndim = 2;
int edim = 2;
double mass = 2000.0;
int init_s = 0;
vector<double> init_r { -5.0, 0.0 };
vector<double> init_v { 20.0 / mass, 0.0 / mass };
int N = 1;
hamiltonian_t hami;

vector<trajectory_t> gen_swarm(int N) {
    /*
     * generate a swarm of trajectories
     */
    vector<complex<double>> init_c(edim, matrixop::ZEROZ);
    init_c[init_s] = matrixop::ONEZ;
    vector<trajectory_t> swarm;
    for (int i(0); i < N; ++i) {
        swarm.emplace_back(hami);
        swarm.back().setup(mass, init_r, init_v, init_c, init_s);
    }
    return swarm;
}


int main(int argc, char** argv) {
    vector<trajectory_t> swarm = gen_swarm(1);

    for (auto& x : swarm[0].get_r()) {
        cout << x << " ";
    }
    cout << endl;
    for (auto& x : swarm[0].get_v()) {
        cout << x << " ";
    }
    cout << endl;
    for (auto& x : swarm[0].get_c()) {
        cout << x.real() << " ";
    }
    cout << endl;
    cout << swarm[0].get_s() << " ";
    cout << endl;
    return 0;
}
