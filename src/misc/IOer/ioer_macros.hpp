#ifndef _IOER_MACROS_HPP
#define _IOER_MACROS_HPP


#if __cplusplus < 201100L
	#error "C++11 surpport is required."
#endif 

namespace ioer 
{
	enum class file_specification {
		UNKNOWN,
		STREAM,
		HDF5
	};
};


#endif
