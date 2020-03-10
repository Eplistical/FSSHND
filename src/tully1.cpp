#include <iostream>
#include <algorithm>
#include <complex>
#include <memory>
#include <vector>
#include "tully1_hamiltonian.hpp"
#include "fssh_trajectory.hpp"
#include "traj_recorder.hpp"
#include "misc/vector.hpp"
#include "misc/matrixop.hpp"
#include "misc/randomer.hpp"
#include "misc/ioer.hpp"
#include "misc/timer.hpp"
#include "misc/MPIer.hpp"
#include "misc/fmtstring.hpp"
#include "boost/program_options.hpp"

using namespace mqc;
using namespace std;
namespace po = boost::program_options;

using hamiltonian_t = Tully1_Hamiltonian;
using trajectory_t = FSSH_Trajectory<hamiltonian_t>;
using recorder_t = traj_recorder<trajectory_t>;

int ndim = 1;
int edim = 2;
double mass = 2000.0;
int init_s = 0;
vector<double> init_r { -6.0 };
vector<double> init_p { 20.0 };
vector<double> sigma_r { 0.0 };
vector<double> sigma_p { 0.0 };
int Ntraj = 2000;
int Nstep = 100000;
int output_step = 1000;
double dt = 0.1;
bool enable_hop = true;
bool enable_berry_force = true;
bool enable_log = true;
vector<double> potential_params;
int seed = 42;
unique_ptr<hamiltonian_t> hami;


void setup_params() {
    /**
     * check & setup parameters for the simulation 
     */
    // check
    misc::confirm<misc::ValueError>(init_r.size() == init_p.size()  and init_r.size() == sigma_r.size() and init_p.size() == sigma_p.size(), 
                                    "argparse: init_r, init_p, sigma_r, sigma_p must have the same sizes.");
    misc::confirm<misc::ValueError>(init_s >= 0, "init_s must >= 0");
    misc::confirm<misc::ValueError>(mass > 0.0, "mass must > 0.");
    misc::confirm<misc::ValueError>(dt > 0.0, "dt must > 0.");
    misc::confirm<misc::ValueError>(Ntraj > 0, "Ntraj must > 0.");
    misc::confirm<misc::ValueError>(Nstep > 0, "Nstep must > 0.");
    misc::confirm<misc::ValueError>(Nstep >= output_step, "Nstep must >= output_step.");
    // setup
    hami = make_unique<hamiltonian_t>();
    if (not potential_params.empty()) {
        hami->set_params(potential_params);
    }
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
        ("mass", po::value<decltype(mass)>(&mass), "mass")
        ("init_r", po::value<decltype(init_r)>(&init_r)->multitoken(), "init_r vector")
        ("init_p", po::value<decltype(init_p)>(&init_p)->multitoken(), "init_p vector")
        ("sigma_r", po::value<decltype(sigma_r)>(&sigma_r)->multitoken(), "sigma_r vector")
        ("sigma_p", po::value<decltype(sigma_p)>(&sigma_p)->multitoken(), "sigma_p vector")
        ("init_s", po::value<decltype(init_s)>(&init_s), "init_s")
        ("potential_params", po::value<decltype(potential_params)>(&potential_params)->multitoken(), "potential_params vector")
        ("enable_hop", po::value<decltype(enable_hop)>(&enable_hop), "enable hop")
        ("enable_berry_force", po::value<decltype(enable_berry_force)>(&enable_berry_force), "enable berry_force")
        ("enable_log", po::value<decltype(enable_log)>(&enable_log), "enable log")
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
    return true;
}

void logging(const string& msg, int rank = 0) {
    /**
     * print out log message
     */
    if (MPIer::rank == rank) {
        ioer::info(msg);
    }
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
        swarm.back().set_enable_hop(enable_hop);
        swarm.back().set_enable_berry_force(enable_berry_force);
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

    // --- setup --- //

    logging("# setting up simulation...");
    setup_params();
    const int my_Ntraj = MPIer::assign_job(Ntraj).size();
    recorder_t recorder;

    // --- simulation --- //

    logging("# initializing trajectories ...");
    vector<trajectory_t> swarm = gen_swarm(my_Ntraj);
    logging("# simulating ...");
    for (int istep(0); istep < Nstep; ++istep) {
        // recording
        if (istep % output_step == 0) {
            logging(misc::fmtstring("# step %d / %d", istep, Nstep));
            recorder.stamp(swarm);
            if (all_of(swarm.begin(), swarm.end(), check_end)) {
                // simulation completes
                for (int irec(istep / output_step); irec < Nstep / output_step; ++irec) {
                    recorder.stamp(swarm);
                }
                logging("# simulation completes.");
                break;
            }
        }
        // propagation
        for (trajectory_t& traj : swarm) {
            if (not check_end(traj)) {
                traj.integrator(dt);
            }
        }
    }
    MPIer::barrier();

    // --- collect data --- // 

    logging("# collecting data ...");
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

    logging("# postprocessing ...");
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
                    " enable_hop = ", enable_hop,
                    " enable_log = ", enable_log,
                    "");
        ioer::tabout("#", "t", 
                "n0T", "n0R", "n1T", "n1R",
                "p0T", "p0R", "p1T", "p1R",
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
            double p0T = 0.0, p0R = 0.0, p1T = 0.0, p1R = 0.0;
            // enegry
            double KE = 0.0, PE = 0.0;
            for (int itraj(0); itraj < Ntraj; ++itraj) {
                if (sarr.at(itraj) == 0) {
                    if (rarr.at(0+itraj*ndim) < 0.0) {
                        n0R += 1.0;
                        p0R += varr.at(0+itraj*ndim);
                    }
                    else {
                        n0T += 1.0;
                        p0T += varr.at(0+itraj*ndim);
                    }
                }
                else if (sarr.at(itraj) == 1) {
                    if (rarr.at(0+itraj*ndim) < 0.0) {
                        n1R += 1.0;
                        p1R += varr.at(0+itraj*ndim);
                    }
                    else {
                        n1T += 1.0;
                        p1T += varr.at(0+itraj*ndim);
                    }
                }
                KE += KEarr.at(itraj);
                PE += PEarr.at(itraj);
            }

            p0R *= (n0R == 0 ? 0.0 : (mass / n0R));
            p0T *= (n0T == 0 ? 0.0 : (mass / n0T));
            p1R *= (n1R == 0 ? 0.0 : (mass / n1R));
            p1T *= (n1T == 0 ? 0.0 : (mass / n1T));

            n0R /= Ntraj;
            n0T /= Ntraj;
            n1R /= Ntraj;
            n1T /= Ntraj;

            KE /= Ntraj;
            PE /= Ntraj;

            ioer::tabout("#", t, 
                    n0T, n0R, n1T, n1R,
                    p0T, p0R, p1T, p1R,
                    KE + PE
                    );
            // final output
            if (irec == Nrec - 1) {
                ioer::tabout(init_p.at(0),
                        n0T, n0R, n1T, n1R,
                        p0T, p0R, p1T, p1R,
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
    randomer::seed(MPIer::assign_random_seed(seed));
    timer::tic();
    try { 
        run();
    } catch(const std::exception& e) { 
        // exception caught
        ioer::info("# Thread ", MPIer::rank, " / ", MPIer::rank, ": Exception caught: ", e.what());
        ioer::info("# ", timer::toc());
        MPIer::abort();
    }
    if (MPIer::master) ioer::info("#", timer::toc());
    if (MPIer::master) ioer::info("#", timer::now());
    MPIer::barrier();
    MPIer::finalize();
    return 0;
}
