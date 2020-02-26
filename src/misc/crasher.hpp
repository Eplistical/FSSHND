#ifndef _CRASHER_HPP
#define _CRASHER_HPP
// module for exception handling

#include <string>
#include <stdexcept>
#include "types.hpp"
#include "exceptions.hpp"

namespace misc {
    
    template <typename ExceptionType = Error>
        inline void confirm(BOOL_T condition, const STRING_T& msg) {
            if (not condition) {
                throw(ExceptionType(msg));
            }
        }


    template <typename ExceptionType = Error>
        inline void crash(BOOL_T condition, const STRING_T& msg) {
            confirm(not condition, msg);
        }

};
#endif
