#ifndef _FERMI_HPP
#define _FERMI_HPP
/* fermi function
 *
 * f(x) = 1 / (1 + e**x) 
 */

#include <cmath>
#include <vector>
#include "types.hpp"

namespace misc {

    inline DOUBLE_T fermi(DOUBLE_T x) noexcept {
        return 1.0 / (1.0 + std::exp(x));
    }

    inline COMPLEX_T fermi(const COMPLEX_T& z) noexcept {
        return 1.0 / (1.0 + std::exp(z));
    }

    template <typename Type>
        inline std::vector<Type> fermi(const std::vector<Type>& x) noexcept {
            const SIZE_T N(x.size());
            std::vector<Type> rst(N);
            for (SIZE_T i(0); i < N; ++i) {
                rst[i] = fermi(x[i]);
            }
            return rst;
        }

};

#endif
