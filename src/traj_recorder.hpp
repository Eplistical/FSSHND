#ifndef _TRAJ_RECORDER_HPP
#define _TRAJ_RECORDER_HPP

/*
 * class to record history of trajectories
 */

#include <complex>
#include <iterator>
#include <vector>
#include "misc/exceptions.hpp"
#include "misc/crasher.hpp"


namespace mqc {


    // --- DECLARATION --- //


    template <typename TrajectoryType>
    class traj_recorder {
        public:
            using trajectory_t = TrajectoryType;

        public:
            // --- ctor/dtor --- //
            traj_recorder();
            ~traj_recorder() = default;

        public:
            // --- modifier --- //
            void clear();
            void stamp(const std::vector<trajectory_t>& swarm);

        private:
            // --- extractor-util --- //
            inline std::vector<double> get_r_ele(int irec, int itraj) const;
            inline std::vector<double> get_v_ele(int irec, int itraj) const;
            inline std::vector<std::complex<double>> get_c_ele(int irec, int itraj) const;
            inline int get_s_ele(int irec, int itraj) const;
            inline double get_t_ele(int irec, int itraj) const;
            inline double get_KE_ele(int irec, int itraj) const;
            inline double get_PE_ele(int irec, int itraj) const;
        public:
            // --- extractor --- //
            std::vector<double> get_r_by_traj(int itraj = -1) const;
            std::vector<double> get_r_by_rec(int irec = -1) const;
            std::vector<double> get_v_by_traj(int itraj = -1) const;
            std::vector<double> get_v_by_rec(int irec = -1) const;
            std::vector<std::complex<double>> get_c_by_traj(int itraj = -1) const;
            std::vector<std::complex<double>> get_c_by_rec(int irec = -1) const;
            std::vector<int> get_s_by_traj(int itraj = -1) const;
            std::vector<int> get_s_by_rec(int irec = -1) const;
            std::vector<double> get_t_by_traj(int itraj = -1) const;
            std::vector<double> get_t_by_rec(int irec = -1) const;
            std::vector<double> get_KE_by_traj(int itraj = -1) const;
            std::vector<double> get_KE_by_rec(int irec = -1) const;
            std::vector<double> get_PE_by_traj(int itraj = -1) const;
            std::vector<double> get_PE_by_rec(int irec = -1) const;
        public:
            // --- getter/setter --- //
            int get_ndim() const noexcept { return m_ndim; }
            int get_edim() const noexcept { return m_edim; }
            int get_Ntraj() const noexcept { return m_Ntraj; }
            int get_Nrec() const noexcept { return m_Nrec; }

        private:
            int m_ndim;
            int m_edim;
            int m_Ntraj;
            int m_Nrec;
            // history format: : traj1(t1), traj2(t1), ..., trajN(t1), traj1(t2), ..., trajN(t2), ...
            std::vector<trajectory_t> m_history; 
    };


    // --- IMPLEMENTATION --- //


    template <typename TrajectoryType> 
        traj_recorder<TrajectoryType>::traj_recorder() 
        : m_ndim(0), m_edim(0), m_Ntraj(0), m_Nrec(0), m_history()
        { }
    

    template <typename TrajectoryType> 
        void traj_recorder<TrajectoryType>::clear() {
            /*
             * clear history
             */
            m_ndim = 0;
            m_edim = 0;
            m_Ntraj = 0;
            m_Nrec = 0;
            std::vector<trajectory_t>().swap(m_history);
        }

    template <typename TrajectoryType> 
        void traj_recorder<TrajectoryType>::stamp(const std::vector<trajectory_t>& swarm) {
            /*
             * record info from a swarm
             */
            // check
            if (m_history.empty()) {
                misc::confirm<misc::ValueError>(swarm.size() > 0, "stamp: empty swarm.");
                m_Ntraj = swarm.size();
                m_ndim = swarm.at(0).get_ndim();
                m_edim = swarm.at(0).get_edim();
            }
            else {
                misc::confirm<misc::ValueError>(m_Ntraj == swarm.size(), "stamp: inconsistent Ntraj.");
                misc::confirm<misc::ValueError>(m_ndim == swarm.at(0).get_ndim(), "stamp: inconsistent ndim.");
                misc::confirm<misc::ValueError>(m_edim == swarm.at(0).get_edim(), "stamp: inconsistent edim.");
            }
            // copy, shrink memory, and add to history
            std::vector<trajectory_t> tmpswarm(swarm);
            for (trajectory_t& traj : tmpswarm) {
                traj.die();
            }
            m_history.insert(m_history.end(), 
                std::make_move_iterator(tmpswarm.begin()), 
                std::make_move_iterator(tmpswarm.end())
            );
            m_Nrec += 1;
        }


    // --- extractor-util --- //


    template <typename TrajectoryType> 
        inline std::vector<double> traj_recorder<TrajectoryType>::get_r_ele(int irec, int itraj) const {
            return m_history.at(itraj+irec*m_Ntraj).get_r();
        }

    template <typename TrajectoryType> 
        inline std::vector<double> traj_recorder<TrajectoryType>::get_v_ele(int irec, int itraj) const {
            return m_history.at(itraj+irec*m_Ntraj).get_v();
        }

    template <typename TrajectoryType> 
        inline std::vector<std::complex<double>> traj_recorder<TrajectoryType>::get_c_ele(int irec, int itraj) const {
            return m_history.at(itraj+irec*m_Ntraj).get_c();
        }

    template <typename TrajectoryType> 
        inline int traj_recorder<TrajectoryType>::get_s_ele(int irec, int itraj) const {
            return m_history.at(itraj+irec*m_Ntraj).get_s();
        }

    template <typename TrajectoryType> 
        inline double traj_recorder<TrajectoryType>::get_t_ele(int irec, int itraj) const {
            return m_history.at(itraj+irec*m_Ntraj).get_t();
        }

    template <typename TrajectoryType> 
        inline double traj_recorder<TrajectoryType>::get_KE_ele(int irec, int itraj) const {
            return m_history.at(itraj+irec*m_Ntraj).cal_KE();
        }

    template <typename TrajectoryType> 
        inline double traj_recorder<TrajectoryType>::get_PE_ele(int irec, int itraj) const {
            return m_history.at(itraj+irec*m_Ntraj).cal_PE();
        }


    // --- extractor --- //


    template <typename TrajectoryType> 
        std::vector<double> traj_recorder<TrajectoryType>::get_r_by_traj(int itraj) const {
            /*
             * extract r for the itraj^th trajectory
             *  if itraj == -1, extract for all trajectories
             *  return format will be : traj1(t1), ..., traj1(tM), ..., trajN(t1), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(itraj < m_Ntraj and itraj >= -1, "traj_recorder: invalid itraj.");
            // extract
            std::vector<double> rst;
            rst.reserve(m_ndim * m_Nrec * (itraj == -1 ? m_Ntraj : 1));
            if (itraj == -1) {
                for (int ktraj(0); ktraj < m_Ntraj; ++ktraj) {
                    auto tmp = get_r_by_traj(ktraj);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int irec(0); irec < m_Nrec; ++irec) {
                    auto tmp = get_r_ele(irec, itraj);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            return rst;
        }

    template <typename TrajectoryType> 
        std::vector<double> traj_recorder<TrajectoryType>::get_r_by_rec(int irec) const {
            /*
             * extract r for the irec^th record
             *  if irec == -1, extract for all records
             *  return format will be : traj1(t1), ..., trajN(t1), ..., traj1(tM), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(irec < m_Nrec and irec >= -1, "rec_recorder: invalid irec.");
            // extract
            std::vector<double> rst;
            rst.reserve(m_ndim * m_Ntraj * (irec == -1 ? m_Nrec : 1));
            if (irec == -1) {
                for (int krec(0); krec < m_Nrec; ++krec) {
                    auto tmp = get_r_by_rec(krec);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int itraj(0); itraj < m_Ntraj; ++itraj) {
                    auto tmp = get_r_ele(irec, itraj);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            return rst;
        }

    template <typename TrajectoryType> 
        std::vector<double> traj_recorder<TrajectoryType>::get_v_by_traj(int itraj) const {
            /*
             * extract v for the itraj^th trajectory
             *  if itraj == -1, extract for all trajectories
             *  return format will be : traj1(t1), ..., traj1(tM), ..., trajN(t1), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(itraj < m_Ntraj and itraj >= -1, "traj_recorder: invalid itraj.");
            // extract
            std::vector<double> rst;
            rst.reserve(m_ndim * m_Nrec * (itraj == -1 ? m_Ntraj : 1));
            if (itraj == -1) {
                for (int ktraj(0); ktraj < m_Ntraj; ++ktraj) {
                    auto tmp = get_v_by_traj(ktraj);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int irec(0); irec < m_Nrec; ++irec) {
                    auto tmp = get_v_ele(irec, itraj);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            return rst;
        }

    template <typename TrajectoryType> 
        std::vector<double> traj_recorder<TrajectoryType>::get_v_by_rec(int irec) const {
            /*
             * extract v for the irec^th record
             *  if irec == -1, extract for all records
             *  return format will be : traj1(t1), ..., trajN(t1), ..., traj1(tM), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(irec < m_Nrec and irec >= -1, "rec_recorder: invalid irec.");
            // extract
            std::vector<double> rst;
            rst.reserve(m_ndim * m_Ntraj * (irec == -1 ? m_Nrec : 1));
            if (irec == -1) {
                for (int krec(0); krec < m_Nrec; ++krec) {
                    auto tmp = get_v_by_rec(krec);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int itraj(0); itraj < m_Ntraj; ++itraj) {
                    auto tmp = get_v_ele(irec, itraj);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            return rst;
        }

    template <typename TrajectoryType> 
        std::vector<std::complex<double>> traj_recorder<TrajectoryType>::get_c_by_traj(int itraj) const {
            /*
             * extract c for the itraj^th trajectory
             *  if itraj == -1, extract for all trajectories
             *  return format will be : traj1(t1), ..., traj1(tM), ..., trajN(t1), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(itraj < m_Ntraj and itraj >= -1, "traj_recorder: invalid itraj.");
            // extract
            std::vector<std::complex<double>> rst;
            rst.reserve(m_edim * m_Nrec * (itraj == -1 ? m_Ntraj : 1));
            if (itraj == -1) {
                for (int ktraj(0); ktraj < m_Ntraj; ++ktraj) {
                    auto tmp = get_c_by_traj(ktraj);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int irec(0); irec < m_Nrec; ++irec) {
                    auto tmp = get_c_ele(irec, itraj);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            return rst;
        }

    template <typename TrajectoryType> 
        std::vector<std::complex<double>> traj_recorder<TrajectoryType>::get_c_by_rec(int irec) const {
            /*
             * extract c for the irec^th record
             *  if irec == -1, extract for all records
             *  return format will be : traj1(t1), ..., trajN(t1), ..., traj1(tM), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(irec < m_Nrec and irec >= -1, "rec_recorder: invalid irec.");
            // extract
            std::vector<std::complex<double>> rst;
            rst.reserve(m_edim * m_Ntraj * (irec == -1 ? m_Nrec : 1));
            if (irec == -1) {
                for (int krec(0); krec < m_Nrec; ++krec) {
                    auto tmp = get_c_by_rec(krec);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int itraj(0); itraj < m_Ntraj; ++itraj) {
                    auto tmp = get_c_ele(irec, itraj);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            return rst;
        }

    template <typename TrajectoryType> 
        std::vector<int> traj_recorder<TrajectoryType>::get_s_by_traj(int itraj) const {
            /*
             * extract s for the itraj^th trajectory
             *  if itraj == -1, extract for all trajectories
             *  return format will be : traj1(t1), ..., traj1(tM), ..., trajN(t1), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(itraj < m_Ntraj and itraj >= -1, "traj_recorder: invalid itraj.");
            // extract
            std::vector<int> rst;
            rst.reserve(m_Nrec * (itraj == -1 ? m_Ntraj : 1));
            if (itraj == -1) {
                for (int ktraj(0); ktraj < m_Ntraj; ++ktraj) {
                    auto tmp = get_s_by_traj(ktraj);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int irec(0); irec < m_Nrec; ++irec) {
                    rst.push_back(get_s_ele(irec, itraj));
                }
            }
            return rst;
        }

    template <typename TrajectoryType> 
        std::vector<int> traj_recorder<TrajectoryType>::get_s_by_rec(int irec) const {
            /*
             * extract s for the irec^th record
             *  if irec == -1, extract for all records
             *  return format will be : traj1(t1), ..., trajN(t1), ..., traj1(tM), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(irec < m_Nrec and irec >= -1, "rec_recorder: invalid irec.");
            // extract
            std::vector<int> rst;
            rst.reserve(m_Ntraj * (irec == -1 ? m_Nrec : 1));
            if (irec == -1) {
                for (int krec(0); krec < m_Nrec; ++krec) {
                    auto tmp = get_s_by_rec(krec);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int itraj(0); itraj < m_Ntraj; ++itraj) {
                    rst.push_back(get_s_ele(irec, itraj));
                }
            }
            return rst;
        }

    template <typename TrajectoryType> 
        std::vector<double> traj_recorder<TrajectoryType>::get_t_by_traj(int itraj) const {
            /*
             * extract t for the itraj^th trajectory
             *  if itraj == -1, extract for all trajectories
             *  return format will be : traj1(t1), ..., traj1(tM), ..., trajN(t1), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(itraj < m_Ntraj and itraj >= -1, "traj_recorder: invalid itraj.");
            // extract
            std::vector<double> rst;
            rst.reserve(m_Nrec * (itraj == -1 ? m_Ntraj : 1));
            if (itraj == -1) {
                for (int ktraj(0); ktraj < m_Ntraj; ++ktraj) {
                    auto tmp = get_t_by_traj(ktraj);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int irec(0); irec < m_Nrec; ++irec) {
                    rst.push_back(get_t_ele(irec, itraj));
                }
            }
            return rst;
        }

    template <typename TrajectoryType> 
        std::vector<double> traj_recorder<TrajectoryType>::get_t_by_rec(int irec) const {
            /*
             * extract t for the irec^th record
             *  if irec == -1, extract for all records
             *  return format will be : traj1(t1), ..., trajN(t1), ..., traj1(tM), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(irec < m_Nrec and irec >= -1, "rec_recorder: invalid irec.");
            // extract
            std::vector<double> rst;
            rst.reserve(m_Ntraj * (irec == -1 ? m_Nrec : 1));
            if (irec == -1) {
                for (int krec(0); krec < m_Nrec; ++krec) {
                    auto tmp = get_t_by_rec(krec);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int itraj(0); itraj < m_Ntraj; ++itraj) {
                    rst.push_back(get_t_ele(irec, itraj));
                }
            }
            return rst;
        }

    template <typename TrajectoryType> 
        std::vector<double> traj_recorder<TrajectoryType>::get_KE_by_traj(int itraj) const {
            /*
             * extract KE for the itraj^th trajectory
             *  if itraj == -1, extract for all trajectories
             *  return format will be : traj1(t1), ..., traj1(tM), ..., trajN(t1), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(itraj < m_Ntraj and itraj >= -1, "traj_recorder: invalid itraj.");
            // extract
            std::vector<double> rst;
            rst.reserve(m_Nrec * (itraj == -1 ? m_Ntraj : 1));
            if (itraj == -1) {
                for (int ktraj(0); ktraj < m_Ntraj; ++ktraj) {
                    auto tmp = get_KE_by_traj(ktraj);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int irec(0); irec < m_Nrec; ++irec) {
                    rst.push_back(get_KE_ele(irec, itraj));
                }
            }
            return rst;
        }

    template <typename TrajectoryType> 
        std::vector<double> traj_recorder<TrajectoryType>::get_KE_by_rec(int irec) const {
            /*
             * extract KE for the irec^th record
             *  if irec == -1, extract for all records
             *  return format will be : traj1(t1), ..., trajN(t1), ..., traj1(tM), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(irec < m_Nrec and irec >= -1, "rec_recorder: invalid irec.");
            // extract
            std::vector<double> rst;
            rst.reserve(m_Ntraj * (irec == -1 ? m_Nrec : 1));
            if (irec == -1) {
                for (int krec(0); krec < m_Nrec; ++krec) {
                    auto tmp = get_KE_by_rec(krec);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int itraj(0); itraj < m_Ntraj; ++itraj) {
                    rst.push_back(get_KE_ele(irec, itraj));
                }
            }
            return rst;
        }

    template <typename TrajectoryType> 
        std::vector<double> traj_recorder<TrajectoryType>::get_PE_by_traj(int itraj) const {
            /*
             * extract PE for the itraj^th trajectory
             *  if itraj == -1, extract for all trajectories
             *  return format will be : traj1(t1), ..., traj1(tM), ..., trajN(t1), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(itraj < m_Ntraj and itraj >= -1, "traj_recorder: invalid itraj.");
            // extract
            std::vector<double> rst;
            rst.reserve(m_Nrec * (itraj == -1 ? m_Ntraj : 1));
            if (itraj == -1) {
                for (int ktraj(0); ktraj < m_Ntraj; ++ktraj) {
                    auto tmp = get_PE_by_traj(ktraj);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int irec(0); irec < m_Nrec; ++irec) {
                    rst.push_back(get_PE_ele(irec, itraj));
                }
            }
            return rst;
        }

    template <typename TrajectoryType> 
        std::vector<double> traj_recorder<TrajectoryType>::get_PE_by_rec(int irec) const {
            /*
             * extract PE for the irec^th record
             *  if irec == -1, extract for all records
             *  return format will be : traj1(t1), ..., trajN(t1), ..., traj1(tM), ..., trajN(tM)
             */
            // check
            misc::confirm<misc::IndexError>(irec < m_Nrec and irec >= -1, "rec_recorder: invalid irec.");
            // extract
            std::vector<double> rst;
            rst.reserve(m_Ntraj * (irec == -1 ? m_Nrec : 1));
            if (irec == -1) {
                for (int krec(0); krec < m_Nrec; ++krec) {
                    auto tmp = get_PE_by_rec(krec);
                    rst.insert(rst.end(), std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
            }
            else {
                for (int itraj(0); itraj < m_Ntraj; ++itraj) {
                    rst.push_back(get_PE_ele(irec, itraj));
                }
            }
            return rst;
        }
} // namespace mqc



#endif // _TRAJ_RECORDER_HPP
