#ifndef _RANDOMER_HPP
#define _RANDOMER_HPP
/* module for random
 * require C++11
 * written based on <random> and <vector>
 *
 * Gaohan
 */
#include <cmath>
#include <random>
#include <functional>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <bitset>
#include "types.hpp"
#include "crasher.hpp"

namespace randomer {
	static std::random_device ran_dev;
	static std::mt19937 rng(ran_dev());

	inline VOID_T seed(std::mt19937::result_type val = rng.default_seed) noexcept {
		rng.seed(val);
	}

	template <typename IntegralType>
        inline typename std::enable_if<std::is_integral<IntegralType>::value, IntegralType>::type
        randint(IntegralType lb = IntegralType(0), IntegralType ub = IntegralType(100)) {
            return std::uniform_int_distribution<IntegralType>(lb, ub)(rng);
        }

	template <typename SizeType, typename IntegralType>
        inline typename std::enable_if<std::is_integral<IntegralType>::value, std::vector<IntegralType>>::type
        vrandint(SizeType N, IntegralType lb = IntegralType(0), IntegralType ub = IntegralType(100)) {
            std::uniform_int_distribution<IntegralType> dis(lb, ub);
            std::vector<IntegralType> rst(N);
            for(SizeType j(0); j < N; ++j) {
                rst[j] = dis(rng);
            }
            return rst;
        }

	template <typename FloatType = DOUBLE_T>
        inline typename std::enable_if<std::is_floating_point<FloatType>::value, FloatType>::type
        rand(FloatType lb = FloatType(0.0), FloatType ub = FloatType(1.0)) {
            return std::uniform_real_distribution<FloatType>(lb, ub)(rng);
        }

	template <typename SizeType, typename FloatType = DOUBLE_T>
        inline typename std::enable_if<std::is_floating_point<FloatType>::value, std::vector<FloatType>>::type
        vrand(SizeType N, FloatType lb = FloatType(0.0), FloatType ub = FloatType(1.0)) {
            std::uniform_real_distribution<FloatType> dis(lb, ub);
            std::vector<FloatType> rst(N);
            for(SizeType j(0); j < N; ++j) {
                rst[j] = dis(rng);
            }
            return rst;
        }

	template <typename FloatType>
        inline typename std::enable_if<std::is_floating_point<FloatType>::value, FloatType>::type
        normal(FloatType mu = FloatType(0.0), FloatType sigma = FloatType(1.0)) {
            return std::normal_distribution<FloatType>(mu, sigma)(rng);
        }

	template <typename SizeType, typename FloatType>
        inline typename std::enable_if<std::is_floating_point<FloatType>::value, std::vector<FloatType>>::type
        vnormal(SizeType N, FloatType mu = FloatType(0.0), FloatType sigma = FloatType(1.0)) {
            std::normal_distribution<FloatType> dis(mu, sigma);
            std::vector<FloatType> rst(N);
            for(SizeType j = 0; j < N; ++j) {
                rst[j] = dis(rng);
            }
            return rst;
	}

	template <typename InputIterator, typename IntegralType = INT_T>
        inline typename std::enable_if<std::is_integral<IntegralType>::value, IntegralType>::type
        discrete(const InputIterator& begin, const InputIterator& end) {
            return std::discrete_distribution<IntegralType>(begin, end)(rng);
        }

	template <typename IntegralType>
        inline typename std::enable_if<std::is_integral<IntegralType>::value, IntegralType>::type
        choice(IntegralType N) noexcept {
            return std::uniform_int_distribution<IntegralType>(0, N-1)(rng);
    }

    template <typename ElementType>
        inline ElementType choice(const std::vector<ElementType>& v) {
            return v.at(choice(v.size()));
        }

	template <typename IntegralType>
        typename std::enable_if<std::is_integral<IntegralType>::value, std::vector<IntegralType>>::type
        __choice_hash(IntegralType N, IntegralType m) {
            // randomly pick m numbers from 0 to N - 1, Hashset way
            std::uniform_int_distribution<> dis(0, N - 1);
            std::unordered_set<IntegralType> record;
            if (2 * m <= N) {
                while (record.size() < m) {
                    record.insert(dis(rng));
                }
                std::vector<IntegralType> rst(record.begin(), record.end());
                return rst;
            }
            else {
                m = N - m;
                while (record.size() < m) {
                    record.insert(dis(rng));
                }
                m = N - m;
                std::vector<IntegralType> rst(m);
                IntegralType j(0);
                for (IntegralType i(0); i < N; ++i) {
                    if (record.find(i) == record.end()) {
                        rst[j] = i;
                        ++j;
                    }
                }
                std::shuffle(rst.begin(), rst.end(), rng);
                return rst;
            }
	}

	template <typename IntegralType>
        typename std::enable_if<std::is_integral<IntegralType>::value, std::vector<IntegralType>>::type
        __choice_bit(IntegralType N, IntegralType m) {
            // randomly pick m numbers from 0 to N - 1, Bitset way
            std::uniform_int_distribution<> dis(0, N - 1);
            constexpr SIZE_T BitSetMax = 1000000;
            if (N <= BitSetMax) {
                std::bitset<BitSetMax> record;
                std::vector<IntegralType> rst(m);

                IntegralType Npick = 2 * m <= N ? m : N - m ;

                IntegralType j(0);
                while (j < Npick) {
                    IntegralType tmp(dis(rng));
                    if (!record.test(tmp)) {
                        record.set(tmp);
                        ++j;
                    }
                }

                j = 0;
                BOOL_T target_flag(2 * m <= N);
                for (IntegralType i(0); i < N; ++i) {
                    if (record.test(i) == target_flag) {
                        rst[j] = i;
                        ++j;
                    }
                }
                std::shuffle(rst.begin(), rst.end(), rng);
                return rst;
            }
        }

	template <typename IntegralType>
        typename std::enable_if<std::is_integral<IntegralType>::value, std::vector<IntegralType>>::type
        choice(IntegralType N, IntegralType m) {
            misc::crash(N >= m, "randomer::choice : N must >= m");
            if (N > 1000000 or static_cast<DOUBLE_T>(N) / m > 1000.0) {
                return __choice_hash(N, m);
            }
            else {
                return __choice_bit(N, m);
            }
        }

    template <typename ElementType, typename IntegralType>
        typename std::enable_if<std::is_integral<IntegralType>::value, std::vector<ElementType>>::type
        choice(const std::vector<ElementType>& v, IntegralType m) {
			std::vector<IntegralType> idx(choice(v.size(), m));
			std::vector<ElementType> rst(m);
			for (IntegralType i(0); i< m; ++i) {
				rst[i] = v.at(idx.at(i));
			}
			return rst;
        }

	template<size_t dim = 3>
		inline std::vector<DOUBLE_T> maxwell_dist(DOUBLE_T mass, DOUBLE_T kT, SIZE_T N = 1) {
			return vnormal(N * dim, 0.0, std::sqrt(kT / mass));
		}


	template <typename FloatType>
        typename std::enable_if<std::is_floating_point<FloatType>::value, std::vector<FloatType>>::type
        MHsample(FloatType (*prob_func)(FloatType), 
            SIZE_T N, 
            SIZE_T Nstep_eql = 10000, 
            SIZE_T Nstep_collect = 100, 
            FloatType x0 = 0.0, 
            FloatType sigma = 1.0
            ) {
        // Sample 1D function with the Metropolis–Hastings algorithm
        FloatType xnow = x0;
        FloatType fxnow = prob_func(xnow);
        FloatType xnext, fxnext;
        // equilibrate
        for (SIZE_T istep(0); istep < Nstep_eql; ++istep) {
            xnext = xnow + randomer::normal(0.0, sigma);
            fxnext = prob_func(xnext);
            if (randomer::rand() < fxnext / fxnow) {
                xnow = xnext;
                fxnow = fxnext;
            }
        }
        // sampling
        std::vector<FloatType> rst;
        rst.reserve(N);
        SIZE_T Nstep_sample = Nstep_collect * N;
        for (SIZE_T istep(0); istep < Nstep_sample; ++istep) {
            xnext = xnow + randomer::normal(0.0, sigma);
            fxnext = prob_func(xnext);
            if (randomer::rand() < fxnext / fxnow) {
                xnow = xnext;
                fxnow = fxnext;
            }
            if (istep % Nstep_collect == 0) {
                rst.push_back(xnow);
            }
        }
        return rst;
    }


	template <typename FloatType, typename StateType>
        typename std::enable_if<std::is_floating_point<FloatType>::value, std::vector<StateType>>::type
        vMHsample(const std::function<FloatType(const StateType&)>& prob_func,
            SIZE_T N, 
            SIZE_T Nstep_eql, 
            SIZE_T Nstep_collect, 
            const StateType& x0, 
            const StateType& sigma) {
        // Sample general function with the Metropolis–Hastings algorithm
        StateType xnow = x0;
        StateType xnext(x0.size());
        FloatType fxnow = prob_func(xnow);
        FloatType fxnext;
        // equilibrate
        for (SIZE_T istep(0); istep < Nstep_eql; ++istep) {
            for (int i(0); i < xnow.size(); ++i) {
                xnext.at(i) = xnow.at(i) + randomer::normal(0.0, sigma.at(i));
            }
            fxnext = prob_func(xnext);
            if (fxnext >= fxnow or randomer::rand() < fxnext / fxnow) {
                xnow = xnext;
                fxnow = fxnext;
            }
        }
        // sampling
        std::vector<StateType> rst;
        rst.reserve(N);
        SIZE_T Nstep_sample = Nstep_collect * N;
        for (SIZE_T istep(0); istep < Nstep_sample; ++istep) {
            for (int i(0); i < xnow.size(); ++i) {
                xnext.at(i) = xnow.at(i) + randomer::normal(0.0, sigma.at(i));
            }
            fxnext = prob_func(xnext);
            if (fxnext >= fxnow or randomer::rand() < fxnext / fxnow) {
                xnow = xnext;
                fxnow = fxnext;
            }
            if (istep % Nstep_collect == 0) {
                rst.push_back(xnow);
            }
        }
        return rst;
    }
}

#endif // _RANDOMER_HPP
