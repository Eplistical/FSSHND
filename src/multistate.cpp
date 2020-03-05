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
#include "misc/fmtstring.hpp"
#include "boost/program_options.hpp"

using namespace mqc;
using namespace std;
namespace po = boost::program_options;

using hamiltonian_t = Multistate_Hamiltonian;
using trajectory_t = FSSH_Trajectory<hamiltonian_t>;
using recorder_t = traj_recorder<trajectory_t>;

int Nsite = 2;
int ndim;
int edim;
double mass = 2000.0;
double kT = 0.1;
double fric_gamma = 0.0;
int init_s = 0;
vector<double> init_r(Nsite, 0.0);
vector<double> init_p(Nsite, 0.0);
vector<double> sigma_r;
vector<double> sigma_p;
int Ntraj = 2000;
int Nstep = 10000;
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
    misc::confirm<misc::ValueError>(fric_gamma >= 0.0, "fric_gamma must >= 0.");
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
        // potential_params set, overwrite mass
        hami->set_params(potential_params);
        mass = hami->get_param("MASS");
    }
    else {
        // potential_params not set, use mass to overwrite param["MASS"] in hamitonian
        hami->set_param("MASS", mass);
    }
    misc::confirm<misc::ValueError>(mass > 0.0, "mass must > 0.");
    sigma_r.assign(ndim, 0.0);
    sigma_p.assign(ndim, sqrt(mass * kT));
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
        ("init_r", po::value<decltype(init_r)>(&init_r)->multitoken(), "init_r vector")
        ("init_p", po::value<decltype(init_p)>(&init_p)->multitoken(), "init_p vector")
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

void logging(const string& msg) {
    /**
     * print out log message
     */
    if (enable_log and MPIer::master) {
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
    return false;
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
                logging("# simulation completes.");
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
    vector<vector<double>> ndarr_data(Nrec);
    for (int irec(0); irec < Nrec; ++irec) {
        auto sarr = recorder.get_s_by_rec(irec);
        auto rarr = recorder.get_r_by_rec(irec);
        auto varr = recorder.get_v_by_rec(irec);
        auto KEarr = recorder.get_KE_by_rec(irec);
        auto PEarr = recorder.get_PE_by_rec(irec);

        // calc diab pop
        vector<double> ndarr;
        auto trajarr = recorder.get_traj_by_rec(irec);
        for (auto& traj : trajarr) {
            auto nd = traj.cal_diab_pop();
            ndarr.insert(ndarr.end(), nd.begin(), nd.end());
        }

        if (MPIer::master) {
            sarr_data.at(irec) = move(sarr);
            rarr_data.at(irec) = move(rarr);
            varr_data.at(irec) = move(varr);
            KEarr_data.at(irec) = move(KEarr);
            PEarr_data.at(irec) = move(PEarr);
            ndarr_data.at(irec) = move(ndarr);
        }

        for (int r(1); r < MPIer::size; ++r) {
            if (MPIer::master) {
                MPIer::recv(r, sarr, rarr, varr, KEarr, PEarr, ndarr);
                sarr_data.at(irec).insert(sarr_data.at(irec).end(), sarr.begin(), sarr.end());
                rarr_data.at(irec).insert(rarr_data.at(irec).end(), rarr.begin(), rarr.end());
                varr_data.at(irec).insert(varr_data.at(irec).end(), varr.begin(), varr.end());
                KEarr_data.at(irec).insert(KEarr_data.at(irec).end(), KEarr.begin(), KEarr.end());
                PEarr_data.at(irec).insert(PEarr_data.at(irec).end(), PEarr.begin(), PEarr.end());
                ndarr_data.at(irec).insert(ndarr_data.at(irec).end(), ndarr.begin(), ndarr.end());
            }
            else if (MPIer::rank == r) {
                MPIer::send(0, sarr, rarr, varr, KEarr, PEarr, ndarr);
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
                "r", vector<string>(ndim - 1, ""),
                "p", vector<string>(ndim - 1, ""),
                "s", 
                "diab_pop", vector<string>(edim - 1, ""),
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
            const auto& ndarr = ndarr_data.at(irec);
            // statistics
            vector<double> r, p;
            vector<double> nd(edim, 0.0);
            double s = 0.0;
            double KE = 0.0, PE = 0.0;
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

                nd += vector<double> (ndarr.begin() + itraj * edim, ndarr.begin() + (itraj+1) * edim);
            }

            r /= Ntraj;
            p /= Ntraj;
            s /= Ntraj;
            KE /= Ntraj;
            PE /= Ntraj;

            nd /= Ntraj;

            ioer::tabout("#", t, r, p, s, nd, KE + PE);
            // final output
            if (irec == Nrec - 1) {
                ioer::tabout(r, p, s, nd, KE + PE);
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
    vector<double> r = init_r;
    ioer::info("# r = ", r);
    ioer::info("# H(r) = ", hami->cal_H(r));
    auto nablaH = hami->cal_nablaH(r);
    ioer::info("# nablaH(r):");
    int i = 0;
    for (auto Hx : nablaH) {
        ioer::info("# ", i, " : ");
        ioer::info("# ", Hx);
        i += 1;
    }

    vector<complex<double>> init_c(edim, matrixop::ZEROZ);
    init_c.at(init_s) = matrixop::ONEZ;
    vector<double> init_v(ndim, 0.0);
    trajectory_t traj(*hami);
    double x = -8.0;
    while (x < 8.0) {
        double y = -8.0;
        while (y < 8.0) {
            vector<double> r{x, y};
            traj.setup(mass, r, init_v, init_c, init_s);
            auto H = hami->cal_H(r);
            auto eva = traj.get_eva();
            ioer::tabout(r, H.at(0+0*Nsite).real(), H.at(1+1*Nsite).real(), eva.at(0), eva.at(1));
            y += 0.1;
        }
        x += 0.1;
    }


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
    randomer::seed(seed);
    timer::tic();
    try { 
        run();
    } catch(const std::exception& e) { 
        // exception caught
        ioer::info("# Thread ", MPIer::rank, " / ", MPIer::rank, ": Exception caught: ", e.what());
        ioer::info("# ", timer::toc());
        MPIer::abort();
    }
    if (MPIer::master) ioer::info("# ", timer::toc());
    MPIer::barrier();
    MPIer::finalize();
    return 0;
}
