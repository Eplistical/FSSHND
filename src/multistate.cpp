#include <iostream>
#include <algorithm>
#include <complex>
#include <memory>
#include <vector>
#include "multistate_hamiltonian.hpp"
#include "fssh_trajectory.hpp"
#include "traj_recorder.hpp"
#include "misc/vector.hpp"
#include "misc/matrixop.hpp"
#include "misc/randomer.hpp"
#include "misc/ioer.hpp"
#include "misc/timer.hpp"
#include "misc/MPIer.hpp"
#include "boost/program_options.hpp"

using namespace mqc;
using namespace std;
namespace po = boost::program_options;

using hamiltonian_t = Multistate_Hamiltonian;
using trajectory_t = FSSH_Trajectory<hamiltonian_t>;
using recorder_t = traj_recorder<trajectory_t>;

int Nsite;
int ndim;
int edim;
double mass;
int init_s = 0;
vector<double> init_r;
vector<double> init_p;
vector<double> sigma_r;
vector<double> sigma_p;
int Ntraj = 2000;
int Nstep = 100000;
int output_step = 1000;
double dt = 0.1;
bool enable_hop = true;
vector<double> potential_params;
int seed = 42;
unique_ptr<hamiltonian_t> hami;

bool setup_params() {
    /*
     * check & setup parameters for the simulation 
     */
    // check
    misc::confirm<misc::ValueError>(Nsite > 1, "Nsite must > 1.");
    misc::confirm<misc::ValueError>(init_r.size() == init_p.size()  and init_r.size() == sigma_r.size() and init_p.size() == sigma_p.size(), 
                                    "argparse: init_r, init_p, sigma_r, sigma_p must have the same sizes.");
    misc::confirm<misc::ValueError>(init_s >= 0, "init_s must >= 0");
    misc::confirm<misc::ValueError>(mass > 0.0, "mass must > 0.");
    misc::confirm<misc::ValueError>(dt > 0.0, "dt must > 0.");
    misc::confirm<misc::ValueError>(Ntraj > 0, "Ntraj must > 0.");
    misc::confirm<misc::ValueError>(Nstep > 0, "Nstep must > 0.");
    misc::confirm<misc::ValueError>(Nstep > output_step, "Nstep must > output_step.");
    misc::confirm<misc::ValueError>(not potential_params.empty(), "potential_params cannot be empty.");
    // setup
    ndim = Nsite;
    edim = Nsite;
    hami = make_unique<hamiltonian_t>(edim);
    hami->set_params(potential_params);
    mass = hami->get_param("MASS");
    return true;
}

bool argparse(int argc, char** argv) 
{
    /**
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
    po::store(parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::notify(vm);    
    // help info
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return false;
    }
    // setup parameters
    return setup_params();
}


vector<trajectory_t> gen_swarm(int Ntraj) {
    /**
     * generate a swarm of trajectories
     */
    vector<complex<double>> init_c(edim, matrixop::ZEROZ);
    init_c.at(init_s) = matrixop::ONEZ;
    vector<trajectory_t> swarm;

    for (int itraj(0); itraj < Ntraj; ++itraj) {
        swarm.emplace_back(*hami);
        vector<double> r(ndim);
        vector<double> v(ndim);
        for (int i(0); i < ndim; ++i) {
            r.at(i) = randomer::normal(init_r.at(i), sigma_r.at(i));
            v.at(i) = randomer::normal(init_p.at(i), sigma_p.at(i)) / mass;
        }
        swarm.back().setup(mass, r, v, init_c, init_s);
    }
    return swarm;
}

bool check_end(const trajectory_t& traj) {
    /**
     * determine if a trajectory reaches the end of its simulation
     */
    return (traj.get_r().at(0) > 8.0 and traj.get_v().at(0) > 0.0) 
        or (traj.get_r().at(0) < -8.0 and traj.get_v().at(0) < 0.0);
}


void run() {
    /**
     * run simulation
     */

    // --- parameters --- //

    const int my_Ntraj = MPIer::assign_job(Ntraj).size();
    recorder_t recorder;

    // --- simulation --- //

    vector<trajectory_t> swarm = gen_swarm(my_Ntraj);
    for (int istep(0); istep < Nstep; ++istep) {
        // recording
        if (istep % output_step == 0) {
            recorder.stamp(swarm);
            if (all_of(swarm.begin(), swarm.end(), check_end)) {
                // simulation completes
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
    MPIer::barrier();

    // --- collect data --- // 

    const int Nrec = recorder.get_Nrec();
    vector<vector<int>> sarr_data(Nrec);
    vector<vector<double>> rarr_data(Nrec);
    vector<vector<double>> varr_data(Nrec);
    vector<vector<double>> KEarr_data(Nrec);
    vector<vector<double>> PEarr_data(Nrec);
    for (int irec(0); irec < Nrec; ++irec) {
        auto sarr = recorder.get_s_by_rec(irec);
        auto rarr = recorder.get_r_by_rec(irec);
        auto varr = recorder.get_v_by_rec(irec);
        auto KEarr = recorder.get_KE_by_rec(irec);
        auto PEarr = recorder.get_PE_by_rec(irec);

        if (MPIer::master) {
            sarr_data.at(irec) = move(sarr);
            rarr_data.at(irec) = move(rarr);
            varr_data.at(irec) = move(varr);
            KEarr_data.at(irec) = move(KEarr);
            PEarr_data.at(irec) = move(PEarr);
        }

        for (int r(1); r < MPIer::size; ++r) {
            if (MPIer::master) {
                MPIer::recv(r, sarr, rarr, varr, KEarr, PEarr);
                sarr_data.at(irec).insert(sarr_data.at(irec).end(), sarr.begin(), sarr.end());
                rarr_data.at(irec).insert(rarr_data.at(irec).end(), rarr.begin(), rarr.end());
                varr_data.at(irec).insert(varr_data.at(irec).end(), varr.begin(), varr.end());
                KEarr_data.at(irec).insert(KEarr_data.at(irec).end(), KEarr.begin(), KEarr.end());
                PEarr_data.at(irec).insert(PEarr_data.at(irec).end(), PEarr.begin(), PEarr.end());
            }
            else if (MPIer::rank == r) {
                MPIer::send(0, sarr, rarr, varr, KEarr, PEarr);
            }
            MPIer::barrier();
        }
    }
    MPIer::barrier(); 

    // --- post-procssing & output --- //

    if (MPIer::master) {
        // header output
        hami->output_params();
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
                    " sigma_p = ", sigma_p,
                    "");
        ioer::tabout("#", "t", 
                "n0T", "n0R", "n1T", "n1R",
                "px0T", "py0T", "px0R", "py0R",
                "px1T", "py1T", "px1R", "py1R",
                "Etot"
                );
        // dynamics output
        for (int irec(0); irec < Nrec; ++irec) {
            double t = irec * output_step * dt;
            const auto& sarr = sarr_data.at(irec);
            const auto& rarr = rarr_data.at(irec);
            const auto& varr = varr_data.at(irec);
            const auto& KEarr = KEarr_data.at(irec);
            const auto& PEarr = PEarr_data.at(irec);
            // population
            double n0T = 0.0, n0R = 0.0, n1T = 0.0, n1R = 0.0;
            // momentum 
            double px0T = 0.0, py0T = 0.0, px0R = 0.0, py0R = 0.0;
            double px1T = 0.0, py1T = 0.0, px1R = 0.0, py1R = 0.0;
            // enegry
            double KE = 0.0, PE = 0.0;
            for (int itraj(0); itraj < Ntraj; ++itraj) {
                if (sarr.at(itraj) == 0) {
                    if (rarr.at(0+itraj*ndim) < 0.0) {
                        n0R += 1.0;
                        px0R += varr.at(0+itraj*ndim);
                        py0R += varr.at(1+itraj*ndim);
                    }
                    else {
                        n0T += 1.0;
                        px0T += varr.at(0+itraj*ndim);
                        py0T += varr.at(1+itraj*ndim);
                    }
                }
                else if (sarr.at(itraj) == 1) {
                    if (rarr.at(0+itraj*ndim) < 0.0) {
                        n1R += 1.0;
                        px1R += varr.at(0+itraj*ndim);
                        py1R += varr.at(1+itraj*ndim);
                    }
                    else {
                        n1T += 1.0;
                        px1T += varr.at(0+itraj*ndim);
                        py1T += varr.at(1+itraj*ndim);
                    }
                }
                KE += KEarr.at(itraj);
                PE += PEarr.at(itraj);
            }

            px0R *= (n0R == 0 ? 0.0 : (mass / n0R));
            py0R *= (n0R == 0 ? 0.0 : (mass / n0R));
            px0T *= (n0T == 0 ? 0.0 : (mass / n0T));
            py0T *= (n0T == 0 ? 0.0 : (mass / n0T));
            px1R *= (n1R == 0 ? 0.0 : (mass / n1R));
            py1R *= (n1R == 0 ? 0.0 : (mass / n1R));
            px1T *= (n1T == 0 ? 0.0 : (mass / n1T));
            py1T *= (n1T == 0 ? 0.0 : (mass / n1T));

            n0R /= Ntraj;
            n0T /= Ntraj;
            n1R /= Ntraj;
            n1T /= Ntraj;

            KE /= Ntraj;
            PE /= Ntraj;

            ioer::tabout("#", t, 
                    n0T, n0R, n1T, n1R,
                    px0T, py0T, px0R, py0R,
                    px1T, py1T, px1R, py1R,
                    KE + PE
                    );
            // final output
            if (irec == Nrec - 1) {
                ioer::tabout(init_p.at(0), t, 
                        n0T, n0R, n1T, n1R,
                        px0T, py0T, px0R, py0R,
                        px1T, py1T, px1R, py1R,
                        KE + PE
                        );
            }
        }
    }
    MPIer::barrier();
}


int main(int argc, char** argv) {
    /**
     * entry
     */
    MPIer::setup();
    if (argparse(argc, argv) == false) {
        return -1;
    }
    randomer::seed(seed);
    if (MPIer::master) timer::tic();
    run();
    if (MPIer::master) ioer::info("# ", timer::toc());
    MPIer::barrier();
    MPIer::finalize();
    return 0;
}