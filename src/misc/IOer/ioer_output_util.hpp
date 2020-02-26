#ifndef _IOER_OUTPUT_T_HPP
#define _IOER_OUTPUT_T_HPP
/* module for output utilities in ioer
 *
 * Gaohan
 */
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <type_traits>
#include "../type_traiter.hpp"
#include "ioer_macros.hpp"
#include "stream_io_mgr.hpp"
#include "ioer_exceptions.hpp"
#include "ioer_output_util_aux.hpp"

namespace ioer 
{
    using std::string;
    using std::size_t;
    using std::ios;
    using std::ios_base;
    using std::iostream;
    using std::ostream;
    using std::streamsize;
    using std::setw;
    using std::setprecision;
    using std::is_arithmetic;
	using std::is_same;
    using type_traiter::is_direct_outputable;
    using type_traiter::is_sequence_container;
    using type_traiter::is_complex;
    using type_traiter::is_fundamental;
    using type_traiter::is_vector;
    using type_traiter::is_deque;
    using type_traiter::is_list;
    using type_traiter::is_forward_list;
    using type_traiter::is_c_string;
    using type_traiter::is_string;
    using type_traiter::is_array_container;

    class output_t 
    {
		private:
			stream_io_mgr& stream_io_mgr_sing = stream_io_mgr::getInstance();

        private:
            string _path;

            size_t _width = 16;
            size_t _precision = 8;
            string _dlm = " ";
            bool _flush = false;

            // functors
            _pair_output_functor_t _pair_output_functor_obj;

        public:
            output_t(const string& path = STDIO_PATH, ios::openmode mode = ios::out) noexcept
                : _path(path), _pair_output_functor_obj(STDIO_PATH, _dlm, _width, _precision, true, false, false)
                { 
                    open(_path, mode | ios::out);
                }

            // no copy 
            output_t(const output_t&) = delete;
            output_t& operator=(const output_t&) = delete;

            // getter
            size_t width(void) const noexcept
            {
                return _width; 
            }

            size_t precision(void) const noexcept
            {
                return _precision; 
            }

            string dlm(void) const noexcept
            {
                return _dlm; 
            }

            bool is_flush(void) const noexcept
            {
                return _flush; 
            }

            // setter
            void set_width(size_t width) noexcept
            {
                _width = width; 
            }

            void set_precision(size_t precision)
            {
                _precision = precision; 
            }

            void set_dlm(const string& dlm) noexcept
            {
                _dlm = dlm; 
            }

            void set_flush(bool flush) noexcept
            {
                _flush = flush; 
            }

            // open/close
            void open(const string& path, ios::openmode mode = ios::out) 
            {
                stream_io_mgr_sing.open(path, mode);
                _path = path;
            }

            void close() noexcept
            {
                stream_io_mgr_sing.close(_path);
                _path = ioer::STDIO_PATH;
            }

            // utilities

            // -- basic -- //
            void newline(void)
            {
                stream_io_mgr_sing.at(_path) << "\n"; 
                if (_flush) stream_io_mgr_sing.at(_path) << flush; 
            }

            void drawline(char c, size_t len = 32) 
            {
                stream_io_mgr_sing.at(_path) << string(len, c) << "\n"; 
                if (_flush) stream_io_mgr_sing.at(_path) << flush; 
            }

			// -- operator<< -- //
			template <typename ParamType>
				typename enable_if<is_direct_outputable<ParamType>::value, output_t&>::type
				operator<<(const ParamType& x) 
				{
                    stream_io_mgr_sing.at(_path) << x;
					return *this;
				}

			template <typename ParamType>
				typename enable_if< is_same<ParamType, ostream>::value || 
									is_same<ParamType, ios>::value || 
									is_same<ParamType, ios_base>::value, 
									output_t&>::type
				operator<<(ParamType& (*x)(ParamType&)) 
				{
                    stream_io_mgr_sing.at(_path) << x;
					return *this;
				}

			output_t& operator<<(streambuf* x) 
			{
				stream_io_mgr_sing.at(_path) << x;
				return *this;
			}

            // -- info -- //
            template <typename ParamType, typename enable_if< is_sequence_container<ParamType>::value, int>::type = 0>
                void _info_helper(const ParamType& x, int) 
                {
                    for (const auto& xi : x) {
                        stream_io_mgr_sing.at(_path) << xi << " ";
                    }
                }

            template <typename ParamType>
                void _info_helper(const ParamType& x, ...) 
                {
                    stream_io_mgr_sing.at(_path) << x;
                }

            template <typename ParamType>
                void _info(const ParamType& x) 
                {
                    _info_helper(x, 0);
                    if (_flush) stream_io_mgr_sing.at(_path) << flush; 
                }

            template <typename ParamType, typename ... Types>
                void _info(const ParamType& x, const Types& ... otherx) 
                {
					_info(x);
                    _info(otherx ...);
                }

            template <typename ... Types>
                void info_nonewline(const Types& ... x) 
                {
                    _info(x ...);
                }

            template <typename ... Types>
                void info(const Types& ... x) 
                {
                    info_nonewline(x ...);
                    newline();
                }

            // -- tabout -- //
            template <typename ParamType, typename enable_if< is_sequence_container<ParamType>::value, int>::type = 0>
                void _tabout_helper(const ParamType& x, int) 
                {
                    for (const auto& xi : x) {
                        stream_io_mgr_sing.at(_path) << setw(_width) << setprecision(_precision) << xi;
                    }
                }

            template <typename ParamType>
                void _tabout_helper(const ParamType& x, ...) 
                {
                    stream_io_mgr_sing.at(_path) << setw(_width) << setprecision(_precision) << x;
                }

            template <typename ParamType>
                void _tabout(const ParamType& x) 
                {
                    _tabout_helper(x, 0);
                    if (_flush) stream_io_mgr_sing.at(_path) << flush; 
                }

            template <typename ParamType, typename ... Types>
                void _tabout(const ParamType& x, const Types& ... otherx) 
                {
                    _tabout(x);
                    _tabout(otherx ...);
                }

            template <typename ... Types>
                void tabout_nonewline(const Types& ... x) 
                {
                    _tabout(x ...);
                }

            template <typename ... Types>
                void tabout(const Types& ... x) 
                {
                    tabout_nonewline(x ...);
                    newline();
                }

            // -- pair output -- //
            const _pair_output_functor_t&
                keyval(bool _keyfirst = true, bool _nonewline = false)
                {
                    _pair_output_functor_obj.config(_path, _dlm, _width, _precision, _keyfirst, _nonewline, _flush); 
                    return _pair_output_functor_obj;
                }

            // -- binary output -- //
            template <typename ParamType, typename enable_if<is_complex<ParamType>::value, int>::type = 0>
                void _write_helper(const ParamType& x, int)
                {
                    using ValType = typename ParamType::value_type;
                    ValType tmp[] = {x.real(), x.imag()};
                    stream_io_mgr_sing.at(_path).write
                        (reinterpret_cast<const char*>(tmp), static_cast<streamsize>(2 * sizeof(ValType)));
                }

            template <typename ParamType, typename enable_if<is_c_string<ParamType>::value, int>::type = 0>
                void _write_helper(const ParamType& x, int)
                {
                    size_t N = 0;
                    while(x[N] != '\0') ++N;
                    stream_io_mgr_sing.at(_path).write
                        (reinterpret_cast<const char*>(&x), static_cast<streamsize>(N * sizeof(char)));
                }

            template <typename ParamType, 
				typename enable_if< is_array_container<ParamType>::value || is_vector<ParamType>::value || is_string<ParamType>::value, int>::type = 0>
                    void _write_helper(const ParamType& x, int)
                    {
                        stream_io_mgr_sing.at(_path).write
                            (reinterpret_cast<const char*>(x.data()), static_cast<streamsize>(x.size() * sizeof(typename ParamType::value_type)));
                    }

            template <typename ParamType,
                typename enable_if<	is_deque<ParamType>::value || is_forward_list<ParamType>::value || is_list<ParamType>::value, int>::type = 0>
					void _write_helper(const ParamType& x, int)
                {
                    for(const auto& it : x) {
                        stream_io_mgr_sing.at(_path).write
                            (reinterpret_cast<const char*>(&x), static_cast<streamsize>(sizeof(typename ParamType::value_type)));
                    }
                }

            template <typename ParamType>
                void _write_helper(const ParamType& x, ...)
                {
                    stream_io_mgr_sing.at(_path).write
                        (reinterpret_cast<const char*>(&x), static_cast<streamsize>(sizeof(x)));
                }


            template <typename ParamType>
                void _write(const ParamType& x) 
                {
                    _write_helper(x, 0);
                }

            template <typename ParamType, typename ... Types>
                void _write(const ParamType& x, const Types& ... otherx) 
                {
                    _write(x);
                    _write(otherx ...);
                }

            template <typename ... Types>
                void write(const Types& ... x) 
                {
                    _write(x ...);
                }
    }; // class output_t

    // -- non-class functions for standard output -- //
    //extern output_t STDOUT;
	static output_t STDOUT;

    inline void newline(void) 
    {
        STDOUT.newline();
    }

    inline void drawline(char c, size_t len = 32)
    {
        STDOUT.drawline(c, len);
    }

    template <typename ... Types>
        inline void info_nonewline(const Types& ... x) 
        {
            STDOUT.info_nonewline(x ...);
        }

    template <typename ... Types>
        inline void info(const Types& ... x) 
        {
            STDOUT.info(x ...);
        }

    template <typename ... Types>
        inline void tabout(const Types& ... x) 
        {
            STDOUT.tabout(x ...);
        }

    template <typename ... Types>
        inline void tabout_nonewline(const Types& ... x) 
        {
            STDOUT.tabout_nonewline(x ...);
        }

    template <typename ... Types>
        inline void write(const Types& ... x) 
        {
            STDOUT.write(x ...);
        }

    inline const _pair_output_functor_t&
        keyval(bool _keyfirst = true, bool _nonewline = false)
        {
            return STDOUT.keyval(_keyfirst, _nonewline);
        }
};

#endif // _IOER_OUTPUT_T_HPP
