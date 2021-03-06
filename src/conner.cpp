#include <iostream>
#include <algorithm>
#include <complex>
#include <memory>
#include <string>
#include <vector>
#include "conner_hamiltonian.hpp"
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

using hamiltonian_t = Conner_Hamiltonian;
using trajectory_t = FSSH_Trajectory<hamiltonian_t>;
using recorder_t = traj_recorder<trajectory_t>;

int Nsite = 2;
int ndim;
int edim;
double mass = 1000.0;
double kT = 1e-3;
double fric_gamma = -1.0;
int init_s = 0;
vector<double> init_r(Nsite, 0.0);
vector<double> init_p(Nsite, 0.0);
vector<double> sigma_r;
vector<double> sigma_p;
int Ntraj = 10000;
int Nstep = 1000000;
int output_step = 100;
double dt = 0.1;
bool enable_hop = true;
bool enable_log = true;
vector<double> potential_params;
int seed = 42;
unique_ptr<hamiltonian_t> hami;

void setup_params() {
    /*
     * check & setup parameters for the simulation 
     */
    // check
    misc::confirm<misc::ValueError>(Nsite > 1, "Nsite must > 1.");
    misc::confirm<misc::ValueError>(kT > 0.0, "kT must > 0.");
    misc::confirm<misc::ValueError>(init_s >= 0, "init_s must >= 0");
    misc::confirm<misc::ValueError>(dt > 0.0, "dt must > 0.");
    misc::confirm<misc::ValueError>(Ntraj > 0, "Ntraj must > 0.");
    misc::confirm<misc::ValueError>(Nstep > 0, "Nstep must > 0.");
    misc::confirm<misc::ValueError>(Nstep >= output_step, "Nstep must >= output_step.");
    // setup
    ndim = Nsite;
    edim = Nsite;
    hami = make_unique<hamiltonian_t>(edim);
    if (not potential_params.empty()) {
        hami->set_params(potential_params);
    }
    misc::confirm<misc::ValueError>(mass > 0.0, "mass must > 0.");

    sigma_r.assign(ndim, 0.5);
    sigma_p.assign(ndim, 1.0);

    /*
    sigma_r.assign(ndim, 0.0);
    sigma_p.assign(ndim, 0.0);
    */
    misc::confirm<misc::ValueError>(init_r.size() == Nsite and init_r.size() == init_p.size()  and init_r.size() == sigma_r.size() and init_p.size() == sigma_p.size(), 
                                    "argparse: init_r, init_p, sigma_r, sigma_p must have the same sizes (which should be Nsite).");
}

bool argparse(int argc, char** argv) 
{
    /**
     * parse cmd arguments
     */
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("Nsite", po::value<decltype(Nsite)>(&Nsite), "# site")
        ("Ntraj", po::value<decltype(Ntraj)>(&Ntraj), "# traj")
        ("Nstep", po::value<decltype(Nstep)>(&Nstep), "# step")
        ("output_step", po::value<decltype(output_step)>(&output_step), "# step for output")
        ("dt", po::value<decltype(dt)>(&dt), "single time step")
        ("mass", po::value<decltype(mass)>(&mass), "mass, will be overwritten if potential_params is set.")
        ("kT", po::value<decltype(kT)>(&kT), "temperature.")
        ("init_r", po::value<decltype(init_r)>(&init_r)->multitoken(), "init r")
        ("init_p", po::value<decltype(init_p)>(&init_p)->multitoken(), "init p")
        ("fric_gamma", po::value<decltype(fric_gamma)>(&fric_gamma), "friction gamma.")
        ("init_s", po::value<decltype(init_s)>(&init_s), "init_s")
        ("potential_params", po::value<decltype(potential_params)>(&potential_params)->multitoken(), "potential_params vector")
        ("enable_hop", po::value<decltype(enable_hop)>(&enable_hop), "enable_hop")
        ("enable_log", po::value<decltype(enable_log)>(&enable_log), "enable_log")
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
    if (enable_log and MPIer::rank == rank) {
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

    vector<double> rrandnum = randomer::vnormal(ndim * Ntraj, 0.0, 1.0);
    vector<double> prandnum = randomer::vnormal(ndim * Ntraj, 0.0, 1.0);

    for (int itraj(0); itraj < Ntraj; ++itraj) {
        swarm.emplace_back(*hami);
        vector<double> r(ndim);
        vector<double> v(ndim);
        for (int i(0); i < ndim; ++i) {
            r.at(i) = randomer::normal(init_r.at(i), sigma_r.at(i));
            v.at(i) = randomer::normal(init_p.at(i), sigma_p.at(i)) / mass;
        }
        swarm.back().setup(mass, r, v, init_c, init_s);
        swarm.back().set_gamma(fric_gamma);
        swarm.back().set_kT(kT);
        swarm.back().set_enable_hop(enable_hop);
    }
    return swarm;
}

bool check_end(const trajectory_t& traj) {
    /**
     * determine if a trajectory reaches the end of its simulation
     */
    const auto r = traj.get_r();
    if (r.at(0) < -hami->get_param("R") - 0.5 * hami->get_param("A")) return true;
    if (r.at(ndim-1) > hami->get_param("R") + 0.5 * hami->get_param("A")) return true;
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
    for (int istep(0); istep <= Nstep; ++istep) {
        // recording
        if (istep % output_step == 0) {
            logging(misc::fmtstring("# step %d / %d", istep, Nstep));
            if (all_of(swarm.begin(), swarm.end(), check_end)) {
                logging("# simulation completes.");
                for (int irec(recorder.get_Nrec()); irec < Nrec; ++irec) {
                    recorder.stamp(swarm);
                }
                break;
            }
            else {
                recorder.stamp(swarm);
            }
        }
        //ioer::info("thread ", MPIer::rank, " istep = ", istep);
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
    misc::confirm<misc::ValueError>(Nrec == recorder.get_Nrec(), misc::fmtstring("Nrec = %d while recorder.get_Nrec() = %d.", Nrec, recorder.get_Nrec()));
    vector<vector<double>> rarr_data(Nrec);
    vector<vector<double>> varr_data(Nrec);
    vector<vector<double>> KEarr_data(Nrec);
    vector<vector<double>> PEarr_data(Nrec);

    for (int irec(0); irec < Nrec; ++irec) {
        auto rarr = recorder.get_r_by_rec(irec);
        auto varr = recorder.get_v_by_rec(irec);
        auto KEarr = recorder.get_KE_by_rec(irec);
        auto PEarr = recorder.get_PE_by_rec(irec);

        if (MPIer::master) {
            rarr_data.at(irec) = move(rarr);
            varr_data.at(irec) = move(varr);
            KEarr_data.at(irec) = move(KEarr);
            PEarr_data.at(irec) = move(PEarr);
        }

        for (int r(1); r < MPIer::size; ++r) {

            if (MPIer::master) {
                MPIer::recv(r, rarr, varr, KEarr, PEarr);
                rarr_data.at(irec).insert(rarr_data.at(irec).end(), rarr.begin(), rarr.end());
                varr_data.at(irec).insert(varr_data.at(irec).end(), varr.begin(), varr.end());
                KEarr_data.at(irec).insert(KEarr_data.at(irec).end(), KEarr.begin(), KEarr.end());
                PEarr_data.at(irec).insert(PEarr_data.at(irec).end(), PEarr.begin(), PEarr.end());

            }
            else if (MPIer::rank == r) {
                MPIer::send(0, rarr, varr, KEarr, PEarr);
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
                    " Nsite = ", Nsite,
                    " ndim = ", ndim,
                    " edim = ", edim,
                    " mass = ", mass,
                    " kT = ", kT,
                    " Ntraj = ", Ntraj,
                    " Nstep = ", Nstep,
                    " output_step = ", output_step,
                    " dt = ", dt,
                    " init_r = ", init_r,
                    " init_p = ", init_p,
                    " sigma_r = ", sigma_r,
                    " sigma_p = ", sigma_p,
                    " fric_gamma = ", fric_gamma,
                    " seed = ", seed,
                    "");
        ioer::tabout("#", "t", 
                "npos",
                "nneg",
                "KEnegi",
                vector<string>(ndim-1, ""),
                "KEposi",
                vector<string>(ndim-1, ""),
                "KEneg",
                "KEpos",
                "Etot"
                );

        // dynamics output
        for (int irec(0); irec < Nrec; ++irec) {
            double t = irec * output_step * dt;
            const auto& rarr = rarr_data.at(irec);
            const auto& varr = varr_data.at(irec);
            const auto& KEarr = KEarr_data.at(irec);
            const auto& PEarr = PEarr_data.at(irec);
            // statistics
            vector<double> r, p;
            double npos = 0.0, nneg = 0.0;
            double KE = 0.0, PE = 0.0;
            double KEneg = 0.0, KEpos = 0.0;
            vector<double> KEposi(ndim, 0.0);
            vector<double> KEnegi(ndim, 0.0);
            for (int itraj(0); itraj < Ntraj; ++itraj) {
                if (rarr.at(ndim-1+itraj*ndim) > hami->get_param("R") + 0.5 * hami->get_param("A")) {
                    npos += 1.0;
                    // KE on different directions
                    for (int i(0); i < ndim; ++i) {
                        KEposi.at(i) += 0.5 * mass * std::pow(varr.at(i+itraj*ndim), 2);
                    }
                    // total KE 
                    KEpos += KEarr.at(itraj);
                }
                else if (rarr.at(0+itraj*ndim) < hami->get_param("R") - 0.5 * hami->get_param("A")) {
                    nneg += 1.0;
                    // KE on different directions
                    for (int i(0); i < ndim; ++i) {
                        KEnegi.at(i) += 0.5 * mass * std::pow(varr.at(i+itraj*ndim), 2);
                    }
                    // total KE 
                    KEneg += KEarr.at(itraj);
                }
                KE += KEarr.at(itraj);
                PE += PEarr.at(itraj);
            }

            if (npos > 0.0) {
                KEposi /= npos;
                KEpos /= npos;
            }
            if (nneg > 0.0) {
                KEnegi /= nneg;
                KEneg /= nneg;
            }
            npos /= Ntraj;
            nneg /= Ntraj;
            KE /= Ntraj;
            PE /= Ntraj;

            ioer::tabout("#", t, npos, nneg, KEnegi, KEposi, KEneg, KEpos, KE + PE);
            // final output
            if (irec == Nrec - 1) {
                ioer::tabout(hami->get_param("W"), npos, nneg, KEnegi, KEposi, KEneg, KEpos, KE + PE);
            }
        }
        /*
        // detailed traj
        ioer::info("## traj detail r:");
        for (int irec(0); irec < Nrec; ++irec) {
            ioer::info("## ", rarr_data.at(irec));
        }
        */

    }
    MPIer::barrier();
}


void runtest() {
    /**
     * run test
     */

    // --- DEBUG --- //


    logging("# setting up simulation...");
    setup_params();
    const int my_Ntraj = MPIer::assign_job(Ntraj).size();
    recorder_t recorder;

    if (MPIer::master) {
        // header output
        hami->output_params();
        ioer::info("# simulation paras: ",  
                    " Nsite = ", Nsite,
                    " ndim = ", ndim,
                    " edim = ", edim,
                    " mass = ", mass,
                    " kT = ", kT,
                    " Ntraj = ", Ntraj,
                    " Nstep = ", Nstep,
                    " output_step = ", output_step,
                    " dt = ", dt,
                    " init_r = ", init_r,
                    " init_p = ", init_p,
                    " sigma_r = ", sigma_r,
                    " sigma_p = ", sigma_p,
                    " fric_gamma et ", fric_gamma,
                    " seed = ", seed,
                    "");
    }

    // hamiltonian
    ioer::info("# hamiltonian: ");
    vector<complex<double>> init_c(edim, matrixop::ZEROZ);
    init_c.at(init_s) = matrixop::ONEZ;
    vector<double> init_v(ndim, 0.0);
    trajectory_t traj(*hami);
    for (double x(-20.0); x < 20.0; x += 0.1) {
        for (double y(-20.0); y < 20.0; y += 0.1) {
            vector<double> r{x, y};
            traj.setup(mass, r, init_v, init_c, init_s);
            auto H = hami->cal_H(r);
            auto eva = traj.get_eva();
            ioer::tabout(r, H.at(0+0*Nsite).real(), H.at(1+1*Nsite).real(), eva.at(0), eva.at(1));
        }
    }

    // init dist
    ioer::info("# initial dist:");
    auto swarm = gen_swarm(Ntraj);
    for (auto& traj : swarm) {
        ioer::info("## ", traj.get_r(), traj.get_v());
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
        run();
        //runtest();
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
