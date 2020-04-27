#include <iostream>
#include <algorithm>
#include <complex>
#include <memory>
#include <vector>
#include "three_state_1d_hamiltonian.hpp"
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
#include "boost/math/quadrature/trapezoidal.hpp"

using namespace mqc;
using namespace std;
namespace po = boost::program_options;

using hamiltonian_t = Three_State_1d_Hamiltonian;
using trajectory_t = FSSH_Trajectory<hamiltonian_t>;
using recorder_t = traj_recorder<trajectory_t>;
using boost::math::quadrature::trapezoidal;

int ndim = 1;
int edim = 3;
double mass = 1000.0;
double kT = 0.1;
double fric_gamma = -1.0;
int init_s = 0;
int Ntraj = 2000;
int Nstep = 100000;
int output_step = 100;
double Etot = 0.03;
double dt = 0.1;
bool enable_hop = true;
bool enable_deco = false;
bool enable_momrev = false;
bool enable_log = true;
vector<double> potential_params;
int seed = 42;
unique_ptr<hamiltonian_t> hami;

void setup_params() {
    /*
     * check & setup parameters for the simulation 
     */
    // check
    misc::confirm<misc::ValueError>(kT > 0.0, "kT must > 0.");
    misc::confirm<misc::ValueError>(init_s >= 0, "init_s must >= 0");
    misc::confirm<misc::ValueError>(dt > 0.0, "dt must > 0.");
    misc::confirm<misc::ValueError>(Ntraj > 0, "Ntraj must > 0.");
    misc::confirm<misc::ValueError>(Nstep > 0, "Nstep must > 0.");
    misc::confirm<misc::ValueError>(Nstep >= output_step, "Nstep must >= output_step.");
    // setup
    hami = make_unique<hamiltonian_t>();
    if (not potential_params.empty()) {
        // potential_params set, overwrite mass
        hami->set_params(potential_params);
    }
    misc::confirm<misc::ValueError>(mass > 0.0, "mass must > 0.");
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
        ("mass", po::value<decltype(mass)>(&mass), "mass, will be overwritten if potential_params is set.")
        ("kT", po::value<decltype(kT)>(&kT), "temperature.")
        ("fric_gamma", po::value<decltype(fric_gamma)>(&fric_gamma), "friction gamma.")
        ("init_s", po::value<decltype(init_s)>(&init_s), "init_s")
        ("Etot", po::value<decltype(Etot)>(&Etot), "Etot")
        ("potential_params", po::value<decltype(potential_params)>(&potential_params)->multitoken(), "potential_params vector")
        ("enable_hop", po::value<decltype(enable_hop)>(&enable_hop), "enable_hop")
        ("enable_deco", po::value<decltype(enable_deco)>(&enable_deco), "enable_deco")
        ("enable_momrev", po::value<decltype(enable_momrev)>(&enable_momrev), "enable_momrev")
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
    //ioer::info("# thread ", MPIer::rank, " : ", msg);
    if (enable_log and MPIer::rank == rank) {
        ioer::info(msg);
    }
}

vector<trajectory_t> gen_swarm(int my_Ntraj) {
    /**
     * generate a swarm of trajectories
     */
    vector<complex<double>> init_c(edim, matrixop::ZEROZ);
    init_c.at(init_s) = matrixop::ONEZ;
    vector<trajectory_t> swarm;
    // calculate vx
    const double vx = std::sqrt(2.0 * Etot / mass);
    // assign initial values
    for (int itraj(0); itraj < my_Ntraj; ++itraj) {
        swarm.emplace_back(*hami);
        vector<double> r(ndim);
        vector<double> v(ndim);
        r.at(0) = -4.0;
        v.at(0) = vx;

        swarm.back().setup(mass, r, v, init_c, init_s);
        swarm.back().set_gamma(fric_gamma);
        swarm.back().set_kT(kT);
        swarm.back().set_enable_hop(enable_hop);
        swarm.back().set_enable_momrev(enable_momrev);
    }
    return swarm;
}

bool check_end(const trajectory_t& traj) {
    /**
     * determine if a trajectory reaches the end of its simulation
     */
    const double x = traj.get_r().at(0);
    const double vx = traj.get_v().at(0);
    return (x < -4.0 and vx < 0.0) or (x > 4.0 and vx > 0.0);
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
        int itr = 0;
        for (trajectory_t& traj : swarm) {
            if (not check_end(traj)) {
                const auto last_v = traj.get_v();

                traj.integrator(dt);

                if (enable_deco) {
                    // simple decoherence
                    const auto v = traj.get_v();
                    if (last_v.at(0) * v.at(0) < 0.0) {
                        vector<complex<double>> newc(edim, 0.0); 
                        newc.at(traj.get_s()) = 1.0;
                        traj.set_c(newc);
                    }
                }
            }
            itr += 1;
        }
    }
    while (recorder.get_Nrec() < Nrec) {
        // appending records to Nrec length
        recorder.stamp(swarm);
    }
    misc::confirm<misc::ValueError>(Nrec == recorder.get_Nrec(), misc::fmtstring("Nrec = %d while recorder.get_Nrec() = %d.", Nrec, recorder.get_Nrec()));

    MPIer::barrier();

    // --- collect data --- // 


    // dynamic info
    logging("# collecting dynamics data ...");
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

    // hop info
    logging("# collecting hop data ...");
    int tot_Nhop_acc = 0;
    int tot_Nhop_fru = 0;
    int Nhop_acc = 0;
    int Nhop_fru = 0;
    for (const auto& traj : swarm) {
        Nhop_acc += traj.get_Nhop_accepted();
        Nhop_fru += traj.get_Nhop_frustrated();
    }
    if (MPIer::master) {
        tot_Nhop_acc = Nhop_acc;
        tot_Nhop_fru = Nhop_fru;
    }

    for (int r(1); r < MPIer::size; ++r) {
        if (MPIer::master) {
            MPIer::recv(r, Nhop_acc, Nhop_fru);
            tot_Nhop_acc += Nhop_acc;
            tot_Nhop_fru += Nhop_fru;
        }
        else if (MPIer::rank == r) {
            MPIer::send(0, Nhop_acc, Nhop_fru);
        }
        MPIer::barrier();
    }

    // --- post-procssing & output --- //

    logging("# postprocessing ...");
    if (MPIer::master) {
        // header output
        hami->output_params();
        ioer::info("# simulation paras: ",  
                    " ndim = ", ndim,
                    " edim = ", edim,
                    " mass = ", mass,
                    " kT = ", kT,
                    " Etot = ", Etot,
                    " Ntraj = ", Ntraj,
                    " Nstep = ", Nstep,
                    " output_step = ", output_step,
                    " dt = ", dt,
                    " fric_gamma = ", fric_gamma,
                    " enable_hop = ", enable_hop,
                    " enable_deco = ", enable_deco,
                    " enable_momrev = ", enable_momrev,
                    " seed = ", seed,
                    "");
        ioer::tabout("#", "t", 
                "nT", vector<string>(edim - 1, ""),
                "nR", vector<string>(edim - 1, ""),
                "KE",
                "PE",
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

            // statistics
            vector<double> r, p;
            double s = 0.0;
            double KE = 0.0, PE = 0.0;
            vector<double> nT(edim, 0.0);
            vector<double> nR(edim, 0.0);
            for (int itraj(0); itraj < Ntraj; ++itraj) {
                s += sarr.at(itraj);
                if (itraj == 0) {
                    r = vector<double>(rarr.begin() + itraj * ndim, rarr.begin() + (itraj + 1) * ndim);
                    p = mass * vector<double>(varr.begin() + itraj * ndim, varr.begin() + (itraj + 1) * ndim);
                }
                else {
                    r += vector<double>(rarr.begin() + itraj * ndim, rarr.begin() + (itraj + 1) * ndim);
                    p += mass * vector<double>(varr.begin() + itraj * ndim, varr.begin() + (itraj + 1) * ndim);
                }
                KE += KEarr.at(itraj);
                PE += PEarr.at(itraj);

                if (varr.at(0 + itraj * ndim) > 0.0) {
                    // nT
                    nT.at(sarr.at(itraj)) += 1.0;
                }
                else {
                    // nR
                    nR.at(sarr.at(itraj)) += 1.0;
                }
            }

            r /= Ntraj;
            p /= Ntraj;
            s /= Ntraj;

            KE /= Ntraj;
            PE /= Ntraj;

            nT /= Ntraj;
            nR /= Ntraj;

            ioer::tabout("#", t, nT, nR, KE, PE, KE + PE);
            // final output
            if (irec == Nrec - 1) {
                // dyn info
                ioer::tabout(Etot, nT, nR, KE, PE, KE + PE);
                // hop info
                ioer::info("## hop info: ", 
                            " Nhop_acc = ", tot_Nhop_acc, 
                            " Nhop_fru = ", tot_Nhop_fru, 
                            " acc ratio = ", static_cast<double>(tot_Nhop_acc) / (tot_Nhop_acc + tot_Nhop_fru),
                            ""
                            );
            }
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
                    " kT = ", kT,
                    " Etot = ", Etot,
                    " Ntraj = ", Ntraj,
                    " Nstep = ", Nstep,
                    " output_step = ", output_step,
                    " dt = ", dt,
                    " fric_gamma = ", fric_gamma,
                    " enable_hop = ", enable_hop,
                    " enable_deco = ", enable_deco,
                    " enable_momrev = ", enable_momrev,
                    " seed = ", seed,
                    "");
    }

    // hamiltonian
    ioer::info("# Hamiltonian:");
    vector<complex<double>> init_c(edim, matrixop::ZEROZ);
    init_c.at(init_s) = matrixop::ONEZ;
    vector<double> init_v(ndim, 0.0);
    trajectory_t traj(*hami);

    ioer::tabout("#0", arange(1, 100));
    ioer::tabout(
                "#x",
                "H00", "H11", "H22",
                "E0", "E1", "E2",
                "F0x", 
                "F1x", 
                "F2x", 
                "d01xR", "d01xI",
                "d02xR", "d02xI",
                "d12xR", "d12xI",
                "H01R", "H01I",
                "H02R", "H02I",
                "H12R", "H12I",
                ""
                );
    double x = -4.0;
    while (x < 4.0) {
        vector<double> r{x};
        traj.setup(mass, r, init_v, init_c, init_s);
        auto H = hami->cal_H(r);
        auto eva = traj.get_eva();
        auto dc = traj.get_dc();
        auto force = traj.get_force();

        ioer::tabout(x, 
                H.at(0+0*edim).real(), H.at(1+1*edim).real(), H.at(2+2*edim).real(),
                eva,

                real(force.at(0).at(0+0*edim)),
                real(force.at(0).at(1+1*edim)), 
                real(force.at(0).at(2+2*edim)), 

                real(dc.at(0).at(0+1*edim)), imag(dc.at(0).at(0+1*edim)),  // d01xR, d01xI
                real(dc.at(0).at(0+2*edim)), imag(dc.at(0).at(0+2*edim)),  // d02xR, d02xI
                real(dc.at(0).at(1+2*edim)), imag(dc.at(0).at(1+2*edim)),  // d12xR, d12xI

                H.at(0+1*edim).real(), H.at(0+1*edim).imag(), // H01R, H01I
                H.at(0+2*edim).real(), H.at(0+2*edim).imag(), // H02R, H02I
                H.at(1+2*edim).real(), H.at(1+2*edim).imag(), // H12R, H12I

                ""
                );
        x += 0.02;
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
        if (argc > 1 and string(argv[1]) == "test") {
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
