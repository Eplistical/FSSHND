#include <iostream>
#include <algorithm>
#include <functional>
#include <complex>
#include <memory>
#include <vector>
#include "nuscatter_hamiltonian.hpp"
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

using hamiltonian_t = Nuscatter_Hamiltonian;
using trajectory_t = FSSH_Trajectory<hamiltonian_t>;
using recorder_t = traj_recorder<trajectory_t>;
using boost::math::quadrature::trapezoidal;

int ndim = 2;
int edim = 2;
double mass = 1000.0;
int init_s = 0;
string init_pos = "bottom";
double init_E = 0.3;
int Ntraj = 2000;
int Nstep = 100000;
int output_step = 1000;
double dt = 0.1;
bool enable_hop = true;
bool enable_log = true;
vector<double> potential_params;
int seed = 42;
unique_ptr<hamiltonian_t> hami;


void setup_params() {
    /**
     * check & setup parameters for the simulation 
     */
    // check
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
        ("init_E", po::value<decltype(init_E)>(&init_E), "initial kinetic energy")
        ("init_pos", po::value<decltype(init_pos)>(&init_pos), "initial position")
        ("init_s", po::value<decltype(init_s)>(&init_s), "inital surface")
        ("potential_params", po::value<decltype(potential_params)>(&potential_params)->multitoken(), "potential_params vector")
        ("enable_hop", po::value<decltype(enable_hop)>(&enable_hop), "enable hop")
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

void logging(const string& msg) {
    /**
     * print out log message
     */
    if (MPIer::master) {
        ioer::info(msg);
    }
}

double xdist(double x) {
    if (x < hami->get_param("r") - 0.5 * hami->get_param("R") or x > 0.5 * hami->get_param("R")) {
        return 0.0;
    }
    else {
        const double L = hami->get_param("R") - hami->get_param("r");
        return sqrt(2.0 / L) * sin(M_PI / L * (x + 0.5 * hami->get_param("R") - hami->get_param("r")));
    }
}

template<typename Callable> 
double wtrans(Callable func, double r, double p, double xmin = -10.0, double xmax = 10.0) {
    auto intr = [&func, &r, &p](double x) {
        return func(r-0.5*x) * func(r+0.5*x) * cos(x*p);
    };
    return trapezoidal(intr, xmin, xmax, 1e-6);
}


vector<vector<double>> rp_sample(size_t N, size_t Nstep_eql, size_t Nstep_collect, const vector<double>& rp0, const vector<double>& rpsigma) {
    // sample r&p from Wigner function
    vector<double> rpnow = rp0;
    double wnow = wtrans(xdist, rpnow[0], rpnow[1]);
    vector<double> rpnext(2);
    double wnext;
    // equilibrate
    for (size_t istep(0); istep < Nstep_eql; ++istep) {
        rpnext[0] = rpnow[0] + randomer::normal(0.0, rpsigma[0]);
        rpnext[1] = rpnow[1] + randomer::normal(0.0, rpsigma[1]);
        wnext = wtrans(xdist, rpnext[0], rpnext[1]);
        if (wnext > 0.0 and (wnext > wnow or randomer::rand() < wnext / wnow)) {
            rpnow = rpnext;
            wnow = wnext;
        }
    }
    // sampling
    vector<vector<double>> rst;
    rst.reserve(N);
    size_t Nstep_sample = Nstep_collect * N;
    for (size_t istep(0); istep < Nstep_sample; ++istep) {
        rpnext[0] = rpnow[0] + randomer::normal(0.0, rpsigma[0]);
        rpnext[1] = rpnow[1] + randomer::normal(0.0, rpsigma[1]);
        wnext = wtrans(xdist, rpnext[0], rpnext[1]);
        if (wnext > 0.0 and (wnext > wnow or randomer::rand() < wnext / wnow)) {
            rpnow = rpnext;
            wnow = wnext;
        }
        if (istep % Nstep_collect == 0) {
            rst.push_back(rpnow);
        }
    }
    return rst;
}

vector<trajectory_t> gen_swarm(int Ntraj) {
    /**
     * generate a swarm of trajectories
     */
    vector<complex<double>> init_c(edim, matrixop::ZEROZ);
    init_c.at(init_s) = matrixop::ONEZ;
    vector<trajectory_t> swarm;
    // sample r&p on the radial diretion
    const vector<double> rp0 { 0.5 * hami->get_param("r"), 0.0 }; 
    const vector<double> rpsigma { 0.1 * (hami->get_param("R") - hami->get_param("r")), 1.0 }; 
    const vector<vector<double>> rps = rp_sample(Ntraj, 10000, 40, rp0, rpsigma);
    // assign initial values
    for (int itraj(0); itraj < Ntraj; ++itraj) {
        swarm.emplace_back(*hami);
        vector<double> r(ndim);
        vector<double> v(ndim);
        if (init_pos == "bottom") {
            r.at(0) = rps.at(itraj).at(0);
            v.at(0) = rps.at(itraj).at(1) / mass;
            r.at(1) = -0.5 * hami->get_param("R");
            v.at(1) = 0.0;
        }
        else if (init_pos == "left") {
            r.at(0) = -0.5 * hami->get_param("R");
            v.at(0) = 0.0;
            r.at(1) = rps.at(itraj).at(0);
            v.at(1) = rps.at(itraj).at(1) / mass;
        }
        swarm.back().setup(mass, r, v, init_c, init_s);
        swarm.back().set_enable_hop(enable_hop);
    }
    return swarm;
}

bool check_end(const trajectory_t& traj) {
    /**
     * determine if a trajectory reaches the end of its simulation
     */
    const double m_theta = hami->cal_m_theta(traj.get_r());
    const vector<double> v = traj.get_v();
    return (m_theta < 0.0 and v.at(1) < 0.0) or (m_theta > M_PI/2 and v.at(0) < 0.0);
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
                    " enable_hop = ", enable_hop,
                    " enable_log = ", enable_log,
                    "");
        ioer::tabout("#", "t", 
                "n0L", "n0B", "n1L", "n1B",
                "px0L", "py0L", "px0B", "py0B",
                "px1L", "py1L", "px1B", "py1B",
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
            double n0L = 0.0, n0B = 0.0, n1L = 0.0, n1B = 0.0;
            // momentum 
            double px0L = 0.0, py0L = 0.0, px0B = 0.0, py0B = 0.0;
            double px1L = 0.0, py1L = 0.0, px1B = 0.0, py1B = 0.0;
            // enegry
            double KE = 0.0, PE = 0.0;
            for (int itraj(0); itraj < Ntraj; ++itraj) {
                if (sarr.at(itraj) == 0) {
                    if (rarr.at(0+itraj*ndim) > rarr.at(1+itraj*ndim)) {
                        n0B += 1.0;
                        px0B += varr.at(0+itraj*ndim);
                        py0B += varr.at(1+itraj*ndim);
                    }
                    else {
                        n0L += 1.0;
                        px0L += varr.at(0+itraj*ndim);
                        py0L += varr.at(1+itraj*ndim);
                    }
                }
                else if (sarr.at(itraj) == 1) {
                    if (rarr.at(0+itraj*ndim) > rarr.at(1+itraj*ndim)) {
                        n1B += 1.0;
                        px1B += varr.at(0+itraj*ndim);
                        py1B += varr.at(1+itraj*ndim);
                    }
                    else {
                        n1L += 1.0;
                        px1L += varr.at(0+itraj*ndim);
                        py1L += varr.at(1+itraj*ndim);
                    }
                }
                KE += KEarr.at(itraj);
                PE += PEarr.at(itraj);
            }

            px0B *= (n0B == 0 ? 0.0 : (mass / n0B));
            py0B *= (n0B == 0 ? 0.0 : (mass / n0B));
            px0L *= (n0L == 0 ? 0.0 : (mass / n0L));
            py0L *= (n0L == 0 ? 0.0 : (mass / n0L));
            px1B *= (n1B == 0 ? 0.0 : (mass / n1B));
            py1B *= (n1B == 0 ? 0.0 : (mass / n1B));
            px1L *= (n1L == 0 ? 0.0 : (mass / n1L));
            py1L *= (n1L == 0 ? 0.0 : (mass / n1L));

            n0B /= Ntraj;
            n0L /= Ntraj;
            n1B /= Ntraj;
            n1L /= Ntraj;

            KE /= Ntraj;
            PE /= Ntraj;

            ioer::tabout("#", t, 
                    n0L, n0B, n1L, n1B,
                    px0L, py0L, px0B, py0B,
                    px1L, py1L, px1B, py1B,
                    KE + PE
                    );
            // final output
            if (irec == Nrec - 1) {
                ioer::tabout(init_E,
                        n0L, n0B, n1L, n1B,
                        px0L, py0L, px0B, py0B,
                        px1L, py1L, px1B, py1B,
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