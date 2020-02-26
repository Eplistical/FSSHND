#ifndef _IOER_STREAM_IO_MGR_HPP
#define _IOER_STREAM_IO_MGR_HPP
// module for ioer::stream_io_mgr
// stream I/O manager class
// THIS IS A SINGLETON

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <memory>
#include "ioer_macros.hpp"
#include "ioer_exceptions.hpp"
#include "io_pathtracker.hpp"

namespace ioer 
{
	using namespace std;
	const string STDIO_PATH("");

	class stream_io_mgr 
	{
		public:
			using path_type = string;
			using file_type = unique_ptr<fstream>;

		private:
			stream_io_mgr() = default;

		public:
			~stream_io_mgr() = default;
			static stream_io_mgr& getInstance() 
			{
				static stream_io_mgr stream_io_mgr_singleton;
				return stream_io_mgr_singleton;
			}

			// no copy
			stream_io_mgr(const stream_io_mgr& other) = delete;
			stream_io_mgr& operator=(const stream_io_mgr& other) = delete;

		public:
			void open(const path_type& path, ios::openmode mode)
			{
				if (not is_stdio(path)) {
					if (path_exist(path))
						return;
					_io_pathtracker_sing.insert(path, file_specification::STREAM);
					_fpdict[path] = unique_ptr<fstream>(new fstream(path, mode));
					if (not _fpdict.at(path)->is_open()){
						ostringstream errmsg;
						errmsg << "ioer::stream_io_mgr::open : Unable to open file "
							<< path << ".";
						throw FileNotOpenedError(errmsg.str());
					}
				}
			}

			void close(const path_type& path) noexcept
			{
				if (path_exist(path)) {
					_fpdict.at(path)->close();
					_fpdict.erase(path);
					_io_pathtracker_sing.erase(path);
				}
			}

		public:
			size_t size() const noexcept
			{
				return _fpdict.size(); 
			}

			iostream& operator[](const path_type& path) 
			{
				return static_cast<iostream&>
					( (is_stdio(path) ? cout : *_fpdict[path]) );
			}

			iostream& at(const path_type& path)
			{
				try {
					return static_cast<iostream&>
						( (is_stdio(path) ? cout : *_fpdict.at(path)) );
				} catch(const out_of_range& e) {
					ostringstream errmsg;
					errmsg << "ioer::stream_io_mgr::at : "
						<< path << " not opened.";
					throw FileNotOpenedError(errmsg.str());
				}
			}


		private:
			bool path_exist(const path_type& path) const noexcept
			{ 
				return (_fpdict.find(path) != _fpdict.end()); 
			}

			bool is_stdio(const path_type& path) const noexcept
			{ 
				return (path == ioer::STDIO_PATH); 
			}


		private:
			unordered_map<path_type, file_type> _fpdict;
			io_pathtracker& _io_pathtracker_sing = io_pathtracker::getInstance();

	}; // class stream_io_mgr

}; // namespace ioer

#endif // _IOER_STREAM_IO_MGR_HPP
