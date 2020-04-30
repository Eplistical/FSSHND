#include <iostream>
#include <algorithm>
#include <complex>
#include <memory>
#include <vector>
#include "fp2bath_hamiltonian.hpp"
#include "fp2bath_trajectory.hpp"
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

using hamiltonian_t = FP2Bath_Hamiltonian;
using trajectory_t = FP2Bath_Trajectory<hamiltonian_t>;
using recorder_t = traj_recorder<trajectory_t>;

int ndim = 2;
int edim = 2;
double mu_l = 0.0;
double mu_r = 0.0;
double gamma_l = 0.005;
double gamma_r = 0.005;
double mass = 333.33;
vector<double> init_r { 0.0, 0.0 };
vector<double> init_p { 0.0, 0.0 };
vector<double> sigma_r { 0.0, 0.0 };
vector<double> sigma_p { 0.0, 0.0 };
double kT = 0.05;
double init_kT = 0.05;
int Ntraj = 2000;
int Nstep = 100000;
int output_step = 1000;
double dt = 1.0;
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
    misc::confirm<misc::ValueError>(mass > 0.0, "mass must > 0.");
    misc::confirm<misc::ValueError>(dt > 0.0, "dt must > 0.");
    misc::confirm<misc::ValueError>(Ntraj > 0, "Ntraj must > 0.");
    misc::confirm<misc::ValueError>(Nstep > 0, "Nstep must > 0.");
    misc::confirm<misc::ValueError>(Nstep >= output_step, "Nstep must >= output_step.");
    // setup hamiltonian parameters
    hami = make_unique<hamiltonian_t>(mu_l, mu_r, gamma_l, gamma_r);
    if (not potential_params.empty()) {
        hami->set_params(potential_params);
    }
    // setup initial distribution parameters
    sigma_r.at(0) = (hami->get_param("OMEGA_X") > 0.0) ? 
        std::sqrt(init_kT / mass / std::pow(hami->get_param("OMEGA_X"), 2)) : 0.0;
    sigma_r.at(1) = (hami->get_param("OMEGA_Y") > 0.0) ? 
        std::sqrt(init_kT / mass / std::pow(hami->get_param("OMEGA_Y"), 2)) : 0.0;
    sigma_p.at(0) = std::sqrt(mass * init_kT);
    sigma_p.at(1) = std::sqrt(mass * init_kT);
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
        ("init_kT", po::value<decltype(init_kT)>(&init_kT), "kT for init distribution")
        ("potential_params", po::value<decltype(potential_params)>(&potential_params)->multitoken(), "potential_params vector")
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
    vector<trajectory_t> swarm;

    for (int itraj(0); itraj < Ntraj; ++itraj) {
        swarm.emplace_back(*hami);
        vector<double> r(ndim);
        vector<double> v(ndim);
        for (int i(0); i < ndim; ++i) {
            r.at(i) = randomer::normal(init_r.at(i), sigma_r.at(i));
            v.at(i) = randomer::normal(init_p.at(i), sigma_p.at(i)) / mass;
        }
        swarm.back().setup(mass, r, v, kT);
    }
    return swarm;
}

bool check_end(const trajectory_t& traj) {
    /**
     * determine if a trajectory reaches the end of its simulation
     */
    return false;
}


void run() {
    /**
     * run simulation
     */

    // --- setup --- //

    logging(misc::fmtstring("# MPI size = %d", MPIer::size));
    logging("# setting up simulation...");
    setup_params();
    const int my_Ntraj = MPIer::assign_job(Ntraj).size();
    recorder_t recorder;

    // --- simulation --- //

    logging("# initializing trajectories ...");
    vector<trajectory_t> swarm = gen_swarm(my_Ntraj);
    logging("# simulating ...");
    const int Nrec = Nstep / output_step + 1;
    for (int istep(0); istep < Nstep; ++istep) {
        // recording
        if (istep % output_step == 0) {
            logging(misc::fmtstring("# step %d / %d", istep, Nstep));
            if (all_of(swarm.begin(), swarm.end(), check_end)) {
                logging("# simulation completes.");
                break;
            }
            else {
                recorder.stamp(swarm);
            }
        }
        // propagation
        for (trajectory_t& traj : swarm) {
            if (not check_end(traj)) {
                traj.integrator(dt);
            }
        }
    }
    while (recorder.get_Nrec() < Nrec) {
        // appending records to Nrec length
        recorder.stamp(swarm);
    }
    misc::confirm<misc::ValueError>(Nrec == recorder.get_Nrec(), misc::fmtstring("Nrec = %d while recorder.get_Nrec() = %d.", Nrec, recorder.get_Nrec()));

    MPIer::barrier();

    // --- collect data --- // 

    logging("# collecting data ...");
    vector<vector<double>> rarr_data(Nrec);
    vector<vector<double>> varr_data(Nrec);
    vector<vector<double>> KEarr_data(Nrec);

    for (int irec(0); irec < Nrec; ++irec) {
        auto rarr = recorder.get_r_by_rec(irec);
        auto varr = recorder.get_v_by_rec(irec);
        auto KEarr = recorder.get_KE_by_rec(irec);

        if (MPIer::master) {
            rarr_data.at(irec) = move(rarr);
            varr_data.at(irec) = move(varr);
            KEarr_data.at(irec) = move(KEarr);
        }

        for (int r(1); r < MPIer::size; ++r) {
            if (MPIer::master) {
                MPIer::recv(r, rarr, varr, KEarr);
                rarr_data.at(irec).insert(rarr_data.at(irec).end(), rarr.begin(), rarr.end());
                varr_data.at(irec).insert(varr_data.at(irec).end(), varr.begin(), varr.end());
                KEarr_data.at(irec).insert(KEarr_data.at(irec).end(), KEarr.begin(), KEarr.end());
            }
            else if (MPIer::rank == r) {
                MPIer::send(0, rarr, varr, KEarr);
            }
            MPIer::barrier();
        }
    }
    MPIer::barrier(); 

    // --- post-procssing & output --- //

    auto cal_N = [](const vector<double>& r) {
        const auto H = hami->cal_H(r);
        const auto nablaH = hami->cal_nablaH(r);
        const double h = H.at(1+1*edim).real() - H.at(0+0*edim).real();
        const double gamma_l = hami->cal_gamma_l(r);
        const double gamma_r = hami->cal_gamma_r(r);
        const double gamma = gamma_l + gamma_r;
        const double mu_l = hami->cal_mu_l();
        const double mu_r = hami->cal_mu_r();
        const double f_l = misc::fermi((h - mu_l) / kT);
        const double f_r = misc::fermi((h - mu_r) / kT);
        const double f = (f_l * gamma_l + f_r * gamma_r) / gamma;
        return f;
    };

    logging("# postprocessing ...");
    if (MPIer::master) {
        // header output
        hami->output_params();
        ioer::info("# simulation paras: ",  
                    " ndim = ", ndim,
                    " edim = ", edim,
                    " mass = ", mass,
                    " mu_l = ", mu_l, 
                    " mu_r = ", mu_r,
                    " gamma_l = ", gamma_l, 
                    " gamma_r = ", gamma_r,
                    " Ntraj = ", Ntraj,
                    " Nstep = ", Nstep,
                    " output_step = ", output_step,
                    " dt = ", dt,
                    " kT = ", kT,
                    " init_r = ", init_r,
                    " init_p = ", init_p,
                    " sigma_r = ", sigma_r,
                    " sigma_p = ", sigma_p,
                    " init_kT = ", init_kT,
                    "");
        ioer::tabout("#", "t", 
                "n",
                "Ek"
                );
        // dynamics output
        for (int irec(0); irec < Nrec; ++irec) {
            double t = irec * output_step * dt;
            const auto& rarr = rarr_data.at(irec);
            const auto& varr = varr_data.at(irec);
            const auto& KEarr = KEarr_data.at(irec);

            // population
            double n0 = 0.0;
            // enegry
            double KE = 0.0;
            for (int itraj(0); itraj < Ntraj; ++itraj) {
                const vector<double> r{ rarr.at(0+itraj*ndim), rarr.at(1+itraj*ndim) };
                n0 += cal_N(r);
                KE += KEarr.at(itraj);
            }

            n0 /= Ntraj;
            KE /= Ntraj;

            ioer::tabout(" ", t, 
                        1.0 - n0, 
                        KE
                        );
            // final output
            if (irec == Nrec - 1) {
            }
        }

        // xv distribution
        ioer::info("## xv distribution:");
        ioer::info("## format: x1, y1,  ..., xN, yN, vx1, vy1, ..., vxN, vyN");
        for (int irec(0); irec < Nrec; ++irec) {
            ioer::tabout("## ", rarr_data.at(irec), varr_data.at(irec));
        }

    }

    MPIer::barrier();
}


void runtest() {
    /**
     * run test
     */

    // --- setup --- //

    logging("# setting up simulation...");
    setup_params();
    const int my_Ntraj = MPIer::assign_job(Ntraj).size();
    recorder_t recorder;


    if (MPIer::master) {
        // header output
        hami->output_params();
        ioer::info("# simulation paras: ",  
                    " ndim = ", ndim,
                    " edim = ", edim,
                    " mass = ", mass,
                    " mu_l = ", mu_l, 
                    " mu_r = ", mu_r,
                    " gamma_l = ", gamma_l, 
                    " gamma_r = ", gamma_r,
                    " Ntraj = ", Ntraj,
                    " Nstep = ", Nstep,
                    " output_step = ", output_step,
                    " dt = ", dt,
                    " kT = ", kT,
                    " init_r = ", init_r,
                    " init_p = ", init_p,
                    " sigma_r = ", sigma_r,
                    " sigma_p = ", sigma_p,
                    " init_kT = ", init_kT,
                    " enable_log = ", enable_log,
                    "");
    }
    
    // hamiltonian
    ioer::info("# Hamiltonian:");
    trajectory_t traj(*hami);
    int Nx = 101;
    vector<double> xarr = linspace(-20.0, 10.0, Nx);
    vector<double> yarr = linspace(-5.0, 5.0, Nx);

    for (int i(0); i < Nx; ++i) {
        double x = xarr.at(i);
        for (int j(0); j < Nx; ++j) {
            double y = yarr.at(j);
            vector<double> r{x, y};
            traj.setup(mass, r, init_p / mass, kT);
            auto H = hami->cal_H(r);
            auto nablaH = hami->cal_nablaH(r);
            ioer::tabout(x, y, 
                        H.at(0+0*edim).real(), H.at(1+1*edim).real(), 
                        traj.get_force(),
                        ""
                        );
        }
    }

    // initial distribution
    ioer::info("# initial dist (r & p):");
    auto swarm = gen_swarm(my_Ntraj);
    for (auto& traj : swarm) {
        ioer::info("## ", traj.get_r(), traj.get_v() * mass);
    }

    return ;
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
        if (string(argv[1]) == "test") {
            runtest();
        }
        else {
            run();
        }
    } catch(const std::exception& e) { 
        // exception caught
        ioer::info("# Thread ", MPIer::rank, " / ", MPIer::size, ": Exception caught: ", e.what());
        ioer::info("# ", timer::toc());
        MPIer::abort();
    }
    if (MPIer::master) ioer::info("#", timer::toc());
    if (MPIer::master) ioer::info("#", timer::now());
    MPIer::barrier();
    MPIer::finalize();
    return 0;
}
