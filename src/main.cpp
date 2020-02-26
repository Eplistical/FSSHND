#include <iostream>
#include <algorithm>
#include <complex>
#include <vector>
#include "tully1_hamiltonian.hpp"
#include "fssh_trajectory.hpp"
#include "misc/matrixop.hpp"
#include "boost/program_options.hpp"

using namespace mqc;
using namespace std;
namespace po = boost::program_options;

using hamiltonian_t = Tully1_Hamiltonian;
using trajectory_t = FSSH_Trajectory<hamiltonian_t>;

int ndim = 2;
int edim = 2;
double mass = 2000.0;
int init_s = 0;
vector<double> init_r { -5.0, 0.0 };
vector<double> init_p { 3.0, 0.0 };
vector<double> sigma_r { 0.0, 0.0 };
vector<double> sigma_p { 0.0, 0.0 };
int Ntraj = 1;
int Nstep = 100000;
int output_step = 100;
double dt = 0.1;
int Nrec = Nstep / output_step;
bool enable_hop = true;
vector<double> potential_params;
int seed = 42;
hamiltonian_t hami;


bool argparse(int argc, char** argv) 
{
    /*
     * parse cmd arguments
     */
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("Ntraj", po::value<decltype(Ntraj)>(&Ntraj), "# traj")
        ("Nstep", po::value<decltype(Nstep)>(&Nstep), "# step")
        ("output_step", po::value<decltype(output_step)>(&output_step), "# step for output")
        ("dt", po::value<decltype(dt)>(&dt), "single time step")
        ("mass", po::value<decltype(mass)>(&mass)->multitoken(), "mass vector")
        ("init_r", po::value<decltype(init_r)>(&init_r)->multitoken(), "init_r vector")
        ("init_p", po::value<decltype(init_p)>(&init_p)->multitoken(), "init_p vector")
        ("sigma_r", po::value<decltype(sigma_r)>(&sigma_r)->multitoken(), "sigma_r vector")
        ("sigma_p", po::value<decltype(sigma_p)>(&sigma_p)->multitoken(), "sigma_p vector")
        ("init_s", po::value<decltype(init_s)>(&init_s), "init_s")
        ("potential_params", po::value<decltype(potential_params)>(&potential_params)->multitoken(), "potential_params vector")
        ("seed", po::value<decltype(seed)>(&seed), "random seed")
        ;
    po::variables_map vm; 
    //po::store(po::parse_command_line(argc, argv, desc), vm);
    po::store(parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::notify(vm);    

    // check

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return false;
    }
    return true;
}



vector<trajectory_t> gen_swarm(int N) {
    /*
     * generate a swarm of trajectories
     */
    vector<complex<double>> init_c(edim, matrixop::ZEROZ);
    init_c[init_s] = matrixop::ONEZ;
    vector<trajectory_t> swarm;
    for (int i(0); i < N; ++i) {
        swarm.emplace_back(hami);
        swarm.back().setup(mass, init_r, init_p / mass, init_c, init_s);
    }
    return swarm;
}

bool check_end(const trajectory_t& traj) {
    /*
     * determine if a trajectory reaches the end of its simulation
     */
    return abs(traj.get_r()[0]) > 6.0;
    //return false;
}


void run() {
    /*
     * run simulation
     */
    vector<trajectory_t> swarm = gen_swarm(Ntraj);
    for (int istep(0); istep < Nstep; ++istep) {
        // propagation
        for (trajectory_t& traj : swarm) {
            if (not check_end(traj)) {
                traj.integrator(dt);
                if (enable_hop) {
                    traj.hopper();
                }
            }
        }
        // recording
        if (istep % output_step) {
            if (all_of(swarm.begin(), swarm.end(), check_end)) {
                // check if all trajectories get end
                break;
            }
            cout << swarm[0].get_r()[0] << endl;
        }
    }


    // --- collecting --- // 


    // --- output --- //


}


int main(int argc, char** argv) {
    run();
    return 0;
}
