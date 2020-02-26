#ifndef _TIMER_HPP
#define _TIMER_HPP
/* module for timer 
 * require C++11
 * written based on <chrono> and <ctime>
 *
 * Gaohan
 */
#include <iostream>
#include <chrono>
#include <ctime>
#include <sstream>
#include <unordered_map>
#include "types.hpp"

namespace timer {
    using namespace std;
	using namespace std::chrono;
	using TicTimeType = unordered_map<INT_T, time_point<system_clock> >;

	static TicTimeType tic_time;

    inline STRING_T s2dhms(DOUBLE_T sec) 
    {
        // convert second to d:h:m:s
        int d, h, m, s, ms;
        d = static_cast<int>(sec / 86400);
        sec -= d * 86400;
        h = static_cast<int>(sec / 3600);
        sec -= h * 3600;
        m = static_cast<int>(sec / 60);
        sec -= m * 60;
        s = static_cast<int>(sec);
        sec -= s;
        ms = static_cast<int>(sec * 1000 + 0.5);
		ostringstream buf;
        if (d > 0) buf << d << " d ";
        if (h > 0) buf << h << " h ";
        if (m > 0) buf << m << " m ";
        if (s > 0) buf << s << " s ";
        buf << ms << " ms";
        return buf.str();
    }

    inline DOUBLE_T elapsed(INT_T index = 0) noexcept 
    {
		duration<DOUBLE_T> rst = system_clock::now() - tic_time[index];
        return rst.count();
    }

	inline VOID_T tic(INT_T index = 0) noexcept
    {
		tic_time[index] = system_clock::now();
	}

	inline STRING_T toc(INT_T index = 0) noexcept
    {
		if (tic_time.find(index) == tic_time.end()) {
			return STRING_T("timer::toc: error!, index not found");
		}
		//duration<DOUBLE_T> elapsed = system_clock::now() - tic_time[index];
		ostringstream buf;
		//buf << "Time elapsed: " << elapsed.count() << " s" ;
		buf << "Time elapsed: " << s2dhms(elapsed(index));
		return buf.str();
	}

	inline STRING_T now() noexcept
    {
		const time_t _now = system_clock::to_time_t(system_clock::now());
		ostringstream buf;
		buf << "now: " << ctime(&_now);
		return buf.str();
	}

};
#endif
