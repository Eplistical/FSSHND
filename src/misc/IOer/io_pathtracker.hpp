#ifndef _IO_PATHTRACKER_HPP
#define _IO_PATHTRACKER_HPP
// module for ioer::io_pathracker
// Track all opened files
// THIS IS A SINGLETON.

#include <string>
#include <sstream>
#include <utility>
#include <unordered_map>
#include "ioer_macros.hpp"
#include "ioer_exceptions.hpp"

namespace ioer 
{

	class io_pathtracker 
	{
		private:
			io_pathtracker() = default;

		public:
			~io_pathtracker() = default;

		public:
			inline static io_pathtracker& getInstance()
			{
				static io_pathtracker io_pathtracker_singleton;
				return io_pathtracker_singleton;
			}

		public:
			size_t size() const noexcept
			{
				return _pathdict.size();
			}

			bool empty() const noexcept
			{
				return _pathdict.empty();
			}

			bool contains(const std::string& path) const noexcept
			{
				return (_pathdict.find(path) != _pathdict.end());
			}

		public:
			void insert(const std::string& path, const file_specification& ftype) 
			{
				if (not contains(path)) {
					_pathdict.insert(std::make_pair(path, ftype));
				}
				else {
					std::ostringstream errmsg;
					errmsg << "ioer::io_pathtracker::insert : "
						<< path << " already opened.";
					throw FileOpenedError(errmsg.str());
				}
			}

			void erase(const std::string& path)
			{
				if (contains(path)) {
					_pathdict.erase(path);
				}
				else {
					std::ostringstream errmsg;
					errmsg << "ioer::io_pathtracker::erase : "
						<< path << " not opened.";
					throw FileNotOpenedError(errmsg.str());
				}
			}

			file_specification ftype(const std::string& path) 
			{
				try {
					return _pathdict.at(path);
				} catch (const std::out_of_range& err) {
					std::ostringstream errmsg;
					errmsg << "ioer::io_pathtracker::ftype : "
						<< path << " not opened.";
					throw FileNotOpenedError(errmsg.str());
				}
			}


		private:
			std::unordered_map<std::string, file_specification> _pathdict;
	};

}; // ioer

#endif
