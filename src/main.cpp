#include <iostream>
#include <algorithm>
#include <complex>
#include <vector>
#include "tully1_hamiltonian.hpp"
#include "fssh_trajectory.hpp"
#include "traj_recorder.hpp"
#include "misc/vector.hpp"
#include "misc/matrixop.hpp"
#include "misc/randomer.hpp"
#include "misc/ioer.hpp"
#include "misc/timer.hpp"
#include "boost/program_options.hpp"

using namespace mqc;
using namespace std;
namespace po = boost::program_options;

using hamiltonian_t = Tully1_Hamiltonian;
using trajectory_t = FSSH_Trajectory<hamiltonian_t>;
using recorder_t = traj_recorder<trajectory_t>;

int ndim = 2;
int edim = 2;
double mass = 2000.0;
int init_s = 0;
vector<double> init_r { -6.0, 0.0 };
vector<double> init_p { 20.0, 0.0 };
vector<double> sigma_r { 0.5, 0.5 };
vector<double> sigma_p { 1.0, 1.0 };
int Ntraj = 5;
int Nstep = 10000;
int output_step = 1000;
double dt = 0.1;
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



vector<trajectory_t> gen_swarm(int Ntraj) {
    /*
     * generate a swarm of trajectories
     */
    vector<complex<double>> init_c(edim, matrixop::ZEROZ);
    init_c[init_s] = matrixop::ONEZ;
    vector<trajectory_t> swarm;

    for (int itraj(0); itraj < Ntraj; ++itraj) {
        swarm.emplace_back(hami);
        vector<double> r(ndim);
        vector<double> v(ndim);
        for (int i(0); i < ndim; ++i) {
            r[i] = randomer::normal(init_r[i], sigma_r[i]);
            v[i] = randomer::normal(init_p[i], sigma_p[i]) / mass;
        }
        swarm.back().setup(mass, r, v, init_c, init_s);
    }
    return swarm;
}

bool check_end(const trajectory_t& traj) {
    /*
     * determine if a trajectory reaches the end of its simulation
     */
    return (traj.get_r()[0] > 8.0 and traj.get_v()[0] > 0.0) 
        or (traj.get_r()[0] < -8.0 and traj.get_v()[0] < 0.0);
    //return false;
}


void run() {
    /*
     * run simulation
     */
    recorder_t recorder;
    vector<trajectory_t> swarm = gen_swarm(Ntraj);
    for (int istep(0); istep < Nstep; ++istep) {
        // recording
        if (istep % output_step == 0) {
            recorder.stamp(swarm);
            if (all_of(swarm.begin(), swarm.end(), check_end)) {
                // check if all trajectories get end
                for (int irec(istep / output_step); irec < Nstep / output_step; ++irec) {
                    recorder.stamp(swarm);
                }
                break;
            }
        }
        // propagation
        for (trajectory_t& traj : swarm) {
            if (not check_end(traj)) {
                traj.integrator(dt);
                if (enable_hop) {
                    traj.hopper(dt);
                }
            }
        }
    }

    // --- collecting --- // 




    // --- output --- //

    // header
    hami.output_params();
    ioer::info("# simulation paras: ",  
    " ndim = ", ndim,
    " edim = ", edim,
    " mass = ", mass,
    " Ntraj = ", Ntraj,
    " Nstep = ", Nstep,
    " output_step = ", output_step,
    " dt = ", dt,
    " init_r = ", init_r,
    " init_p = ", init_p,
    " sigma_r = ", sigma_r,
    ""
    );
    ioer::tabout("# t", "n0T", "n0R", "n1T", "n1R");

    int Nrec = recorder.get_Nrec();
    auto tarr = recorder.get_t_by_traj(0);
    for (int irec(0); irec < Nrec; ++irec) {
        auto sarr = recorder.get_s_by_rec(irec);
        auto rarr = recorder.get_r_by_rec(irec);
        double n0T = 0.0;
        double n0R = 0.0;
        double n1T = 0.0;
        double n1R = 0.0;
        for (int itraj(0); itraj < Ntraj; ++itraj) {
            if (sarr[itraj] == 0) {
                if (rarr[0+itraj*ndim] < 0.0) {
                    n0R += 1.0;
                }
                else {
                    n0T += 1.0;
                }
            }
            else if (sarr[itraj] == 1) {
                if (rarr[0+itraj*ndim] < 0.0) {
                    n1R += 1.0;
                }
                else {
                    n1T += 1.0;
                }
            }
        }
        n0R /= Ntraj;
        n0T /= Ntraj;
        n1R /= Ntraj;
        n1T /= Ntraj;
        ioer::tabout(tarr[irec], n0T, n0R, n1T, n1R);
    }
}


int main(int argc, char** argv) {
    // parse args
    if (argparse(argc, argv) == false) {
        return -1;
    }
    randomer::seed(seed);
    timer::tic();
    run();
    ioer::info("# ", timer::toc());
    return 0;
}
