#ifndef _IOER_OUTPUT_T_AUX_HPP
#define _IOER_OUTPUT_T_AUX_HPP
/* aux module for output utilities in ioer
 *
 * Gaohan
 */
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include "../type_traiter.hpp"
#include "ioer_macros.hpp"
#include "stream_io_mgr.hpp"
#include "ioer_exceptions.hpp"

namespace ioer 
{
    using std::string;
    using std::size_t;
    using std::iostream;
    using type_traiter::is_direct_outputable;
    using type_traiter::is_sequence_container;

    class _pair_output_functor_t 
    {
		private:
			stream_io_mgr& stream_io_mgr_sing = stream_io_mgr::getInstance();

        public:
            _pair_output_functor_t(const string&, const string&, size_t, size_t, bool, bool, bool) noexcept;

        public:
            void config(const string&, const string&, size_t, size_t, bool, bool, bool) noexcept;

            template<typename KeyType, typename ValType>
                typename enable_if< is_direct_outputable< ValType >::value, 
                         const _pair_output_functor_t&>::type
                operator()(const KeyType& key, const ValType& val) const;

            template<typename KeyType, typename ValType>
                typename enable_if< is_sequence_container<ValType>::value,
                         const _pair_output_functor_t& >::type
                operator()(const KeyType& key, const ValType& val) const;

        protected:
            string _path; 
            string _dlm;
            size_t _width;
            size_t _precision;
            bool _keyfirst;
            bool _nonewline;
            bool _flush;

    }; // class _pair_output_functor_t

    /* 
     * Implimentation
     */

    /*****************************************************************************/

    _pair_output_functor_t::_pair_output_functor_t(const string& path, 
            const string& dlm, size_t width, size_t precision, 
            bool keyfirst, bool nonewline, 
            bool flush) noexcept : 
        _path(path), _dlm(dlm), 
        _width(width), _precision(precision), 
        _keyfirst(keyfirst), 
        _nonewline(nonewline), _flush(flush) 
        {}

    /*****************************************************************************/

    void _pair_output_functor_t::config(const string& path, 
            const string& dlm, size_t width, size_t precision,
            bool keyfirst, bool nonewline,
            bool flush) noexcept
    {
        _path = path;
        _dlm = dlm;
        _width = width;
        _precision = precision;
        _keyfirst = keyfirst;
        _nonewline = nonewline;
        _flush = flush;
    }

    /*****************************************************************************/

    template<typename KeyType, typename ValType>
        typename enable_if< is_direct_outputable< ValType >::value, 
                 const _pair_output_functor_t&>::type
        _pair_output_functor_t::operator()(const KeyType& key, const ValType& val) const
        {
            iostream& dest = stream_io_mgr_sing.at(_path);
            if (_keyfirst) 
                dest << setw(_width) << setprecision(_precision) << key << _dlm 
					<< setw(_width) << setprecision(_precision) << val;
            else 
                dest << setw(_width) << setprecision(_precision) << val << _dlm 
					<< setw(_width) << setprecision(_precision) << key;
            if (!_nonewline) dest <<  "\n";
            if (_flush) dest << flush;
            return *this;
        }

    /*****************************************************************************/

    template<typename KeyType, typename ValType>
        typename enable_if< is_sequence_container<ValType>::value,
                 const _pair_output_functor_t& >::type
        _pair_output_functor_t::operator()(const KeyType& key, const ValType& val) const
			{
				iostream& dest = stream_io_mgr_sing.at(_path);
				if(_keyfirst) {
					dest <<  setw(_width) << setprecision(_precision) << key << _dlm;
					for(auto& it : val) dest <<  setw(_width) << setprecision(_precision) << it;
				}
				else{
					for(auto& it : val) dest <<  setw(_width) << setprecision(_precision) << it;
					dest <<  _dlm << setw(_width) << setprecision(_precision) << key ;
				}
				if(!_nonewline) dest <<  "\n";
				if (_flush) dest << flush;
				return *this;
			}

    /*****************************************************************************/

};

#endif // _IOER_OUTPUT_T_HPP
