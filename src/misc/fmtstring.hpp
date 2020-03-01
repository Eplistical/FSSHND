#ifndef _FMTSTRING_HPP
#define _FMTSTRING_HPP
/* fmtstring function
 *
 * generate a formatted string with sprintf, and return the string
 *
 */

#include <cstdio>
#include <string>
#include "types.hpp"

namespace misc {
    
   template <typename ... Types>
       STRING_T fmtstring(const CHAR_T* fmt, Types&&... x)
       {
           static CHAR_T buf[512];
           std::sprintf(buf, fmt, std::forward<Types>(x) ...);
           return STRING_T(buf);
       }
}

#endif
