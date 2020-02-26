#ifndef _IOER_INPUT_T_HPP
#define _IOER_INPUT_T_HPP
/* module for input utilities in ioer
 *
 * Gaohan
 */
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include "../type_traiter.hpp"
#include "ioer_macros.hpp"
#include "stream_io_mgr.hpp"
#include "ioer_exceptions.hpp"

namespace ioer
{
	using namespace std;
	using namespace type_traiter;

	class input_t 
	{
		private:
			stream_io_mgr& stream_io_mgr_sing = stream_io_mgr::getInstance();

		private:
			string _path;
            // for formatted input
            string _raw_text;
            map<string, string> _paradict;


		public:
			input_t(const string& path = STDIO_PATH, ios::openmode mode = ios::in) noexcept
				: 	_path(path)
			{ 
				stream_io_mgr_sing.open(path, mode | ios::in);
			}

			// no copy 
			input_t(const input_t&) = delete;
			input_t& operator=(const input_t&) = delete;
			
			// open/close
			void open(const string& path, ios::openmode mode = ios::in) 
			{
				stream_io_mgr_sing.open(path, mode);
				_path = path;
			}

			void close(void) 
			{
				stream_io_mgr_sing.close(_path);
				_path = STDIO_PATH;
			}

			// utilities

			// -- binary input -- //
			template<typename ParamType, typename enable_if<is_complex<ParamType>::value, int>::type = 0>
				void _read_helper(ParamType& x, int)
				{
					using ValType = typename ParamType::value_type;
					ValType tmp[2];
					stream_io_mgr_sing.at(_path).read
						(reinterpret_cast<char*>(tmp), static_cast<std::streamsize>(2 * sizeof(ValType)));
					x.real(tmp[0]);
					x.imag(tmp[1]);
				}

			template<typename ParamType,
				typename enable_if<	is_array_container<ParamType>::value || is_vector<ParamType>::value || is_string<ParamType>::value, int>::type = 0>
					void _read_helper(ParamType& x, int)
					{
						stream_io_mgr_sing.at(_path).read
							(reinterpret_cast<char*>(&x[0]), static_cast<std::streamsize>(x.size() * sizeof(typename ParamType::value_type)));
					}

			template<typename ParamType,
				typename enable_if<	is_deque<ParamType>::value || is_forward_list<ParamType>::value || is_list<ParamType>::value, int>::type = 0>
					void _read_helper(ParamType& x, int)
				{
					for(auto& it : x) {
						stream_io_mgr_sing.at(_path).read
							(reinterpret_cast<char*>(&x), static_cast<std::streamsize>(sizeof(typename ParamType::value_type)));
					}
				}

			template<typename ParamType>
				void _read_helper(ParamType& x, ...)
				{
					stream_io_mgr_sing.at(_path).read
						(reinterpret_cast<char*>(&x), static_cast<std::streamsize>(sizeof(x)));
				}


			template<typename ParamType>
				void _read(ParamType& x) 
				{
					_read_helper(x, 0);
				}

			template<typename ParamType, typename ... Types>
				void _read(ParamType& x, Types& ... otherx) 
				{
					_read(x);
					_read(otherx ...);
				}

			template<typename ... Types>
				void read(Types& ... x) 
				{
					_read(x ...);
				}

			// -- formatted input -- //

            void read_text(char commentsign = '#')
            {
                // only applicable for format as following:
                // 
                // value        # key
                //
                // 200          # mass
                // # this is a comment
                //
                // 2.3          # velocity
                // 1.out        # outfile

                string line, key, val;
                bool firstchar(true), val_done(false);
                _raw_text.clear();

                while (!stream_io_mgr_sing.at(_path).eof()) {
                    getline(stream_io_mgr_sing.at(_path), line);
                    _raw_text += line;
                    _raw_text += "\n";

                    // a blank line
                    if (!line.empty()) {
                        key.clear();
                        val.clear();
                        firstchar = true;
                        val_done = false;

                        for (char c : line) {
                            if (c != ' ' and c != '\t') {
                                // a comment line
                                if (firstchar and c == commentsign) {
                                    break;
                                }
                                // line w/ key - value pair
                                else {
                                    if (!val_done and c == commentsign) {
                                        val_done = true;
                                    }
                                    else if (val_done) {
                                        key.push_back(c);
                                    }
                                    else {
                                        val.push_back(c);
                                    }
                                }
                                firstchar = false;
                            }
                        }
                        // add record to para dict
                        if (!key.empty()) {
                            _paradict[key] = val;
                        }
                    }
                }
            }

            string get_raw_text() 
            {
                return _raw_text;
            }

			template<typename ParamType>
                typename enable_if<is_integral<ParamType>::value, void>::type
                extract_para(const string& key, ParamType& val) 
                {
                    val = static_cast<ParamType>(stoll(_paradict.at(key)));
                }

			template<typename ParamType>
                typename enable_if<is_floating_point<ParamType>::value, void>::type
                extract_para(const string& key, ParamType& val) 
                {
                    val = static_cast<ParamType>(stod(_paradict.at(key)));
                }

			template<typename ParamType>
                typename enable_if<is_same<ParamType, string>::value, void>::type
                extract_para(const string& key, ParamType& val) 
                {
                    val = _paradict.at(key);
                }
            
	}; // class input_t
};
#endif
