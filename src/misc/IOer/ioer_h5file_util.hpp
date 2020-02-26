#ifndef _IOER_H5FILE_UTIL_HPP
#define _IOER_H5FILE_UTIL_HPP
/* module for hdf5 file utilities in ioer
 *
 * Gaohan
 */

#ifdef IOER_WITH_HDF5

#include <iostream>
#include <vector>
#include <cassert>
#include <map>
#include <utility>
#include "../type_traiter.hpp"
#include "../types.hpp"
#include "ioer_macros.hpp"
#include "ioer_exceptions.hpp"
#include "hdf5_io_mgr.hpp"

namespace ioer
{
	using namespace std;
	using namespace type_traiter;

	static const map<type_index, H5::DataType> h5_dtype_dict
		= {
			{typeid(CHAR_T), H5::PredType::NATIVE_CHAR},
			{typeid(INT_T), H5::PredType::NATIVE_INT},
			{typeid(UNSIGNED_T), H5::PredType::NATIVE_UINT},
			{typeid(INT32_T), H5::PredType::NATIVE_INT32},
			{typeid(INT64_T), H5::PredType::NATIVE_INT64},
			{typeid(UINT32_T), H5::PredType::NATIVE_UINT32},
			{typeid(UINT64_T), H5::PredType::NATIVE_UINT64},
			{typeid(LONG_T), H5::PredType::NATIVE_LONG},
			{typeid(LONG_LONG_T), H5::PredType::NATIVE_LLONG},
			{typeid(UNSIGNED_LONG_T), H5::PredType::NATIVE_ULONG},
			{typeid(UNSIGNED_LONG_LONG_T), H5::PredType::NATIVE_ULLONG},
			{typeid(DOUBLE_T), H5::PredType::NATIVE_DOUBLE},
			{typeid(LONG_DOUBLE_T), H5::PredType::NATIVE_LDOUBLE}
		};

	class h5file_t 
	{
		public:
			using path_type = hdf5_io_mgr::path_type;
			using file_type = hdf5_io_mgr::file_type;
			using key_type = string;

		public:
			h5file_t(const string& path, ios::openmode mode)
			{ 
				open(path, mode);
			}

			// no copy 
			h5file_t(const h5file_t&) = delete;
			h5file_t& operator=(const h5file_t&) = delete;

		public:
			// open/close
			void open(const path_type& path, ios::openmode mode) 
			{
				_path = path;
				_mode = mode;
				_hdf5_io_mgr_sing.open(path, mode);
			}

			void close(void) 
			{
				_hdf5_io_mgr_sing.close(_path);
			}

		public:
			// -- create_dataset -- //
			template <typename DType>
				void _create_dataset(const key_type& ds_key, const vector<DType>& buf) {
					if (not buf.empty()) {
						hsize_t dim(buf.size());
						H5::DataType dtype(getH5PredType(buf.at(0)));
						H5::DataSpace dspace(1, &dim);
						H5::DataSet ds(_hdf5_io_mgr_sing.at(_path).createDataSet(ds_key, dtype, dspace));
						ds.write(buf.data(), dtype);
					}
				}

			template <typename DType, typename ... Types>
				void _create_dataset(const key_type& ds_key, const vector<DType>& buf, const Types& ... others) {
					_create_dataset(ds_key, buf);
					_create_dataset(others ...);
				}

			template <typename ... Types>
				void create_dataset(const Types& ... x) {
					_create_dataset(x ...);
				}


			// -- create_attr -- //
			template <typename AType>
				typename enable_if<is_arithmetic<AType>::value and
										not is_c_string<AType>::value, void>::type
				_create_attr(const key_type& ds_key, const key_type& attr_key, const AType& val) {
					H5::DataType dtype(getH5PredType(val));
					H5::DataSpace dspace(H5S_SCALAR);
					_hdf5_io_mgr_sing.at(_path).openDataSet(ds_key).createAttribute(attr_key, dtype, dspace).write(dtype, &val);
				}

			template <typename AType>
				typename enable_if<is_string<AType>::value or
										is_c_string<AType>::value, void>::type
				_create_attr(const key_type& ds_key, const key_type& attr_key, const AType& val) {
					string buf(val);
					H5::StrType dtype(PredType::C_S1, buf.size());
					H5::DataSpace dspace(H5S_SCALAR);
					_hdf5_io_mgr_sing.at(_path).openDataSet(ds_key).createAttribute(attr_key, dtype, dspace).write(dtype, buf);
				}

			template <typename AType, typename ... Types>
				void _create_attr(const key_type& ds_key, const key_type& attr_key, const AType& val, const Types& ... others) {
					_create_attr(ds_key, attr_key, val);
					_create_attr(ds_key, others ...);
				}

			template <typename ... Types>
				void create_attr(const key_type& ds_key, const Types& ... x) {
					_create_attr(ds_key, x ...);
				}

		public:
			// -- read_dataset -- //
			template <typename DType>
				void _read_dataset(const key_type& ds_key, vector<DType>& buf) {
					static vector<hsize_t> dims;
					H5::DataSet ds = _hdf5_io_mgr_sing.at(_path).openDataSet(ds_key);
					H5::DataSpace dataspace = ds.getSpace();
					dims.resize(dataspace.getSimpleExtentNdims());
					dataspace.getSimpleExtentDims(&dims[0]);
					// read data, only support 1D
					assert(dims.size() == 1);
					buf.resize(dims[0]);
					ds.read(buf.data(), getH5PredType(buf[0]));
				}

			template <typename DType, typename ... Types>
				void _read_dataset(const key_type& ds_key, vector<DType>& buf, Types& ... others) {
					_read_dataset(ds_key, buf);
					_read_dataset(others ...);
				}

			template <typename ... Types>
				void read_dataset(Types& ... x) {
					_read_dataset(x ...);
				}

			// -- read_attr -- //
			void _read_attr(const key_type& ds_key) {}

			template <typename AType>
				typename enable_if<is_arithmetic<AType>::value, void>::type
				_read_attr(const key_type& ds_key, const key_type& attr_key, AType& val) {
					H5::DataType dtype(getH5PredType(val));
					_hdf5_io_mgr_sing.at(_path).openDataSet(ds_key).openAttribute(attr_key).read(dtype, &val);
				}

			template <typename AType>
				typename enable_if<is_string<AType>::value, void>::type
				_read_attr(const key_type& ds_key, const key_type& attr_key, AType& val) {
					H5::Attribute attr = _hdf5_io_mgr_sing.at(_path).openDataSet(ds_key).openAttribute(attr_key);
					H5::StrType dtype(PredType::C_S1, attr.getStorageSize());
					H5std_string buf; 
					attr.read(dtype, buf);
					val = buf;
				}

			template <typename AType, typename ... Types>
				void _read_attr(const key_type& ds_key, const key_type& attr_key, AType& val, Types& ... others) {
					_read_attr(ds_key, attr_key, val);
					_read_attr(ds_key, others ...);
				}

			template <typename ... Types>
				void read_attr(const key_type& ds_key, Types& ... x) {
					_read_attr(ds_key, x ...);
				}

		private:
			template <typename DType>
				const H5::DataType& getH5PredType(const DType& val) {
					return h5_dtype_dict.at(typeid(val));
				}

		private:
			hdf5_io_mgr& _hdf5_io_mgr_sing = hdf5_io_mgr::getInstance();
			path_type _path;
			ios::openmode _mode;


	}; // class h5file_t
};

#endif // ifdef IOER_WITH_HDF5
#endif
