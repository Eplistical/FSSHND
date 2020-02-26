#ifndef _IOER_HDF5_IO_MGR_HPP
#define _IOER_HDF5_IO_MGR_HPP
// module for ioer::hdf5_io_mgr
// stream I/O manager class
// THIS IS A SINGLETON

#ifdef IOER_WITH_HDF5

#include <string>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include "ioer_macros.hpp"
#include "ioer_exceptions.hpp"
#include "io_pathtracker.hpp"
#include "H5Cpp.h"

namespace ioer 
{
	using namespace std;
	using namespace H5;

	class hdf5_io_mgr 
	{
		public:
			using file_type = H5::H5File;
			using path_type = std::string;

		private:
			hdf5_io_mgr() = default;

		public:
			~hdf5_io_mgr() = default;
			static hdf5_io_mgr& getInstance() 
			{
				static hdf5_io_mgr hdf5_io_mgr_singleton;
				return hdf5_io_mgr_singleton;
			}

			// no copy
			hdf5_io_mgr(const hdf5_io_mgr& other) = delete;
			hdf5_io_mgr& operator=(const hdf5_io_mgr& other) = delete;

		public:
			void open(const path_type& path, ios::openmode mode)
			{
				if (not path_exist(path)) {
					_io_pathtracker_sing.insert(path, file_specification::HDF5);

					if (mode & ios::in and mode & ios::out) {
						_fpdict[path] = H5File(path, H5F_ACC_RDWR);
					}
					else if (mode & ios::in) {
						_fpdict[path] = H5File(path, H5F_ACC_RDONLY);
					}
					else if (mode & ios::out) {
						_fpdict[path] = H5File(path, H5F_ACC_TRUNC);
					}
				}
			}

			void close(const path_type& path) noexcept
			{
				if (path_exist(path)) {
					_fpdict.at(path).close();
					_fpdict.erase(path);
					_io_pathtracker_sing.erase(path);
				}
			}

		public:
			size_t size() const noexcept
			{
				return _fpdict.size(); 
			}

			file_type& operator[](const path_type& path) 
			{
				return _fpdict[path];
			}

			file_type& at(const path_type& path)
			{
				try {
					return _fpdict.at(path);
				} catch (const std::out_of_range& e) {
					ostringstream errmsg;
					errmsg << "ioer::hdf5_io_mgr::at : "
						<< path << " not opened.";
					throw FileNotOpenedError(errmsg.str());
				}
			}


		protected:
			bool path_exist(const path_type& path) const noexcept
			{ 
				return (_fpdict.find(path) != _fpdict.end()); 
			}


		private:
			unordered_map<path_type, file_type> _fpdict;
			io_pathtracker& _io_pathtracker_sing = io_pathtracker::getInstance();

	}; // class hdf5_io_mgr

}; // namespace ioer

#endif // ifdef IOER_WITH_HDF5
#endif // ifndef _IOER_HDF5_IO_MGR_HPP
