#ifndef _VECTOR_HPP
#define _VECTOR_HPP
/* module for std::vector operations
 * require C++11
 * written based on <vector>
 *
 * Gaohan
 */
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <stdexcept>
#include "types.hpp"
#include "type_traiter.hpp"
#include "crasher.hpp"

// helpers

template <typename T1, typename T2>
inline void assertsize(const std::vector<T1>& v1, const std::vector<T2>& v2)
{
	misc::confirm<std::runtime_error>(v1.size() == v2.size(), 
			"assertsize: v1 & v2 have different sizes!");
}

template <typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, T >::type
conj(const T& a) {
	return a;
}

// element-wise logical operations

template <typename T>
std::vector<BOOL_T> logical_not(const std::vector<T>& v) 
{
	const SIZE_T N(v.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = static_cast<BOOL_T>(!v[j]);
	return rst;
}

template <typename T1, typename T2>
std::vector<BOOL_T> logical_and(const std::vector<T1>& v1, const std::vector<T2>& v2) 
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = static_cast<BOOL_T>(v1[j] && v2[j]);
	return rst;
}

template <typename T1, typename T2>
std::vector<BOOL_T> logical_and(const std::vector<T1>& v1, const T2& a) 
{
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = static_cast<BOOL_T>(v1[j] && a);
	return rst;
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> logical_and(const T2& a, const std::vector<T1>& v1) 
{
	return logical_and(v1, a);
}

template <typename T1, typename T2>
std::vector<BOOL_T> logical_or(const std::vector<T1>& v1, const std::vector<T2>& v2) 
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = static_cast<BOOL_T>(v1[j] || v2[j]);
	return rst;
}

template <typename T1, typename T2>
std::vector<BOOL_T> logical_or(const std::vector<T1>& v1, const T2& a) 
{
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = static_cast<BOOL_T>(v1[j] || a);
	return rst;
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> logical_or(const T2& a, const std::vector<T1>& v1) 
{
	return logical_or(v1, a);
}

// element-wise relational operations

template <typename T1, typename T2>
std::vector<BOOL_T> equal(const std::vector<T1>& v1, const std::vector<T2>& v2)
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = (v1[j] == v2[j]);
	return rst;
}

template <typename T1, typename T2>
std::vector<BOOL_T> equal(const std::vector<T1>& v1, const T2& a)
{
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = (v1[j] == a);
	return rst;
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> equal(const T2& a, const std::vector<T1>& v1)
{
	return equal(v1, a);
}

template <typename T1, typename T2>
std::vector<BOOL_T> not_equal(const std::vector<T1>& v1, const std::vector<T2>& v2)
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = (v1[j] != v2[j]);
	return rst;
}

template <typename T1, typename T2>
std::vector<BOOL_T> not_equal(const std::vector<T1>& v1, const T2& a)
{
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = (v1[j] != a);
	return rst;
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> not_equal(const T2& a, const std::vector<T1>& v1)
{
	return not_equal(v1, a);
}

template <typename T1, typename T2>
std::vector<BOOL_T> greater_equal(const std::vector<T1>& v1, const std::vector<T2>& v2)
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = (v1[j] >= v2[j]);
	return rst;
}

template <typename T1, typename T2>
std::vector<BOOL_T> greater_equal(const std::vector<T1>& v1, const T2& a)
{
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = (v1[j] >= a);
	return rst;
}

template <typename T1, typename T2>
std::vector<BOOL_T> less_equal(const std::vector<T1>& v1, const std::vector<T2>& v2)
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = (v1[j] <= v2[j]);
	return rst;
}

template <typename T1, typename T2>
std::vector<BOOL_T> less_equal(const std::vector<T1>& v1, const T2& a)
{
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = (v1[j] <= a);
	return rst;
}

template <typename T1, typename T2>
std::vector<BOOL_T> greater(const std::vector<T1>& v1, const std::vector<T2>& v2)
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = (v1[j] > v2[j]);
	return rst;
}

template <typename T1, typename T2>
std::vector<BOOL_T> greater(const std::vector<T1>& v1, const T2& a)
{
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = (v1[j] > a);
	return rst;
}

template <typename T1, typename T2>
std::vector<BOOL_T> less(const std::vector<T1>& v1, const std::vector<T2>& v2)
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = (v1[j] < v2[j]);
	return rst;
}

template <typename T1, typename T2>
std::vector<BOOL_T> less(const std::vector<T1>& v1, const T2& a)
{
	const SIZE_T N(v1.size());
	std::vector<BOOL_T> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = (v1[j] < a);
	return rst;
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> greater_equal(const T2& a, const std::vector<T1>& v1)
{
	return less_equal(v1, a);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> less_equal(const T2& a, const std::vector<T1>& v1)
{
	return greater_equal(v1, a);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> greater(const T2& a, const std::vector<T1>& v1)
{
	return less(v1, a);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> less(const T2& a, const std::vector<T1>& v1)
{
	return greater(v1, a);
}

// element-wise math operations

template <typename T>
typename std::enable_if< std::is_arithmetic<T>::value, std::vector<T> >::type
real(const std::vector< std::complex<T> >& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::real(v[j]);
	return rst;
}

template <typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, std::vector<T> >::type
real(const std::vector<T>& v)
{
	return v;
}

template <typename T>
typename std::enable_if< std::is_arithmetic<T>::value, std::vector<T> >::type
imag(const std::vector< std::complex<T> >& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::imag(v[j]);
	return rst;
}

template <typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, std::vector<T> >::type
imag(const std::vector<T>& v)
{
	return std::vector<T>(v.size(), static_cast<T>(0.0));
}

template <typename T>
typename std::enable_if< std::is_arithmetic<T>::value, std::vector< std::complex<T> > >::type
conj(const std::vector< std::complex<T> >& v)
{
	const SIZE_T N(v.size());
	std::vector< std::complex<T> > rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::conj(v[j]);
	return rst;
}

template <typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, std::vector<T> >::type
conj(const std::vector<T>& v)
{
	return v;
}

template <typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, std::vector<T> >::type
abs(const std::vector< std::complex<T> >& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::abs(v[j]);
	return rst;
}

template <typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, std::vector<T> >::type
abs(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::abs(v[j]);
	return rst;
}

template <typename T>
std::vector<T> sqrt(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::sqrt(v[j]);
	return rst;
}

template <typename T1, typename T2>
auto pow(const std::vector<T1>& v1, const std::vector<T2>& v2) -> std::vector<decltype(std::pow(v1.at(0), v2.at(0)))>
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<decltype(std::pow(v1.at(0), v2.at(0)))> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = std::pow(v1[j], v2[j]);
	return rst;
}

template <typename T1, typename T2>
auto pow(const std::vector<T1>& v1, const T2& a) -> std::vector<decltype(std::pow(v1.at(0), a))>
{
	const SIZE_T N(v1.size());
	std::vector<decltype(std::pow(v1.at(0), a))> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = std::pow(v1[j], a);
	return rst;
}

template <typename T1, typename T2>
auto pow(const T2& a, const std::vector<T1>& v1) -> std::vector<decltype(std::pow(a, v1.at(0)))>
{
	const SIZE_T N(v1.size());
	std::vector<decltype(a, std::pow(v1.at(0)))> rst(N);
	for(size_t j(0); j < N; ++j)
		rst[j] = std::pow(a, v1[j]);
	return rst;
}

template <typename T>
std::vector<T> exp(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::exp(v[j]);
	return rst;
}

template <typename T>
std::vector<T> log(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::log(v[j]);
	return rst;
}

template <typename T>
std::vector<T> log2(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::log2(v[j]);
	return rst;
}

template <typename T>
std::vector<T> log10(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::log10(v[j]);
	return rst;
}

template <typename T>
std::vector<T> sin(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::sin(v[j]);
	return rst;
}

template <typename T>
std::vector<T> cos(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::cos(v[j]);
	return rst;
}

template <typename T>
std::vector<T> tan(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::tan(v[j]);
	return rst;
}

template <typename T>
std::vector<T> asin(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::asin(v[j]);
	return rst;
}

template <typename T>
std::vector<T> acos(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::acos(v[j]);
	return rst;
}

template <typename T>
std::vector<T> atan(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::atan(v[j]);
	return rst;
}

template <typename T>
std::vector<T> sinh(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::sinh(v[j]);
	return rst;
}

template <typename T>
std::vector<T> cosh(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::cosh(v[j]);
	return rst;
}

template <typename T>
std::vector<T> tanh(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j) 
		rst[j] = std::tanh(v[j]);
	return rst;
}

// collective operations

template <typename T>
inline BOOL_T all(const std::vector<T>& v)
{
	// see numpy.all
	return std::all_of(v.begin(), v.end(), [](const T& vi) { return (static_cast<BOOL_T>(vi)); });
}

template <typename T>
inline BOOL_T any(const std::vector<T>& v)
{
	// see numpy.any
	return std::any_of(v.begin(), v.end(), [](const T& vi) { return (static_cast<BOOL_T>(vi)); });
}

template <typename T>
T sum(const std::vector<T>& v)
{
	// return v[0] + v[1] + ... + v[N-1]
	T rst = static_cast<T>(0);
	const SIZE_T N(v.size());
	for(SIZE_T j(0); j < N; ++j) 
		rst += v.at(j);
	return rst;
}

template <typename T>
inline T mean(const std::vector<T>& v)
{
	return sum(v) / v.size();
}

template <typename T>
T prod(const std::vector<T>& v)
{
	T rst = static_cast<T>(1);
	const SIZE_T N(v.size());
	for(SIZE_T j(0); j < N; ++j) 
		rst *= v[j];
	return rst;
}

template <typename T>
typename std::enable_if< std::is_integral<T>::value, DOUBLE_T>::type
norm2(const std::vector<T>& v)
{
	DOUBLE_T rst(0.0);
	const SIZE_T N(v.size());
	for(SIZE_T j(0); j < N; ++j)
		rst += v[j] * v[j];
	return rst;
}

template <typename T>
inline typename std::enable_if< std::is_integral<T>::value, DOUBLE_T>::type
norm(const std::vector<T>& v) 
{
	return std::sqrt(norm2(v));
}

template <typename T>
typename std::enable_if< std::is_floating_point<T>::value, T>::type
norm2(const std::vector<T>& v)
{
	T rst(0.0);
	const SIZE_T N(v.size());
	for(SIZE_T j(0); j < N; ++j)
		rst += v[j] * v[j];
	return rst;
}

template <typename T>
inline typename std::enable_if< std::is_floating_point<T>::value, T>::type
norm(const std::vector<T>& v) 
{
	return std::sqrt(norm2(v));
}

template <typename T>
typename std::enable_if< std::is_integral<T>::value, DOUBLE_T>::type
norm2(const std::vector< std::complex<T> >& v)
{
	DOUBLE_T rst(0.0);
	const SIZE_T N(v.size());
	for(SIZE_T j(0); j < N; ++j)
		rst += std::pow(v[j].real(), 2) + pow(v[j].imag(), 2);
	return rst;
}

template <typename T>
inline typename std::enable_if< std::is_integral<T>::value, DOUBLE_T>::type
norm(const std::vector< std::complex<T> >& v)
{
	return std::sqrt(norm2(v));
}

template <typename T>
typename std::enable_if< std::is_floating_point<T>::value, T>::type
norm2(const std::vector< std::complex<T> >& v)
{
	T rst(0.0);
	const size_t n(v.size());
	for(size_t j(0); j < n; ++j)
		rst += std::pow(v[j].real(), 2) + pow(v[j].imag(), 2);
	return rst;
}

template <typename T>
inline typename std::enable_if< std::is_floating_point<T>::value, T>::type
norm(const std::vector< std::complex<T> >& v) 
{
	return std::sqrt(norm2(v));
}

template <typename T>
T min(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	T rst = v[0];
	for(SIZE_T j(1); j < N; ++j)
		rst = (v[j] < rst)?v[j]:rst;
	return rst;
}

template <typename T>
T max(const std::vector<T>& v)
{
	const SIZE_T N(v.size());
	T rst = v[0];
	for(SIZE_T j(1); j < N; ++j)
		rst = (v[j] > rst)?v[j]:rst;
	return rst;
}

// arithmetic operators 

template <typename T1, typename T2>
auto operator+(const std::vector<T1>& v1, const std::vector<T2>& v2) -> std::vector<decltype(v1.at(0) + v2.at(0))>
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<decltype(v1.at(0) + v2.at(0))> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = v1[j] + v2[j];
	return rst;
}

template <typename T1, typename T2>
auto operator+(const std::vector<T1>& v1, const T2& a) -> std::vector<decltype(v1.at(0) + a)>
{
	const SIZE_T N(v1.size());
	std::vector<decltype(v1.at(0) + a)> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = v1[j] + a;
	return rst;
}

template <typename T1, typename T2>
auto operator+(const T2& a, const std::vector<T1>& v1) -> std::vector<decltype(a + v1.at(0))>
{
	return v1 + a;
}

template <typename T1, typename T2>
std::vector<T1>& operator+=(std::vector<T1>& v1, const T2& a)
{
	const SIZE_T N(v1.size());
	for(SIZE_T j(0); j < N; ++j)
		v1[j] += a;
	return v1;
}

template <typename T1, typename T2>
std::vector<T1>& operator+=(std::vector<T1>& v1, const std::vector<T2>& v2)
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	for(SIZE_T j(0); j < N; ++j)
		v1[j] += v2[j];
	return v1;
}

template <typename T1>
std::vector<T1>& operator++(std::vector<T1>& v1)
{
	const SIZE_T N(v1.size());
	for(SIZE_T j(0); j < N; ++j)
		v1[j] += static_cast<T1>(1.0);
	return v1;
}

template <typename T1>
inline const std::vector<T1> operator++(std::vector<T1>& v1, int)
{
	std::vector<T1> rst(v1);
	++v1;
	return rst;
}

template <typename T1, typename T2>
auto operator-(const std::vector<T1>& v1, const std::vector<T2>& v2) -> std::vector<decltype(v1.at(0) - v2.at(0))>
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<decltype(v1.at(0) - v2.at(0))> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = v1[j] - v2[j];
	return rst;
}

template <typename T1, typename T2>
auto operator-(const std::vector<T1>& v1, const T2& a) -> std::vector<decltype(v1.at(0) - a)>
{
	const SIZE_T N(v1.size());
	std::vector<decltype(v1.at(0) - a)> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = v1[j] - a;
	return rst;
}

template <typename T1, typename T2>
auto operator-(const T2& a, const std::vector<T1>& v1) -> std::vector<decltype(a - v1.at(0))>
{
	const SIZE_T N(v1.size());
	std::vector<decltype(a - v1.at(0))> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = a - v1[j];
	return rst;
}

template <typename T1, typename T2>
std::vector<T1>& operator-=(std::vector<T1>& v1, const T2& a)
{
	const SIZE_T N(v1.size());
	for(SIZE_T j(0); j < N; ++j)
		v1[j] -= a;
	return v1;
}

template <typename T1, typename T2>
std::vector<T1>& operator-=(std::vector<T1>& v1, const std::vector<T2>& v2)
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	for(SIZE_T j(0); j < N; ++j)
		v1[j] -= v2[j];
	return v1;
}

template <typename T1>
std::vector<T1>& operator--(std::vector<T1>& v1)
{
	const SIZE_T N(v1.size());
	for(SIZE_T j(0); j < N; ++j)
		v1[j] -= static_cast<T1>(1.0);
	return v1;
}

template <typename T1>
inline const std::vector<T1> operator--(std::vector<T1>& v1, int)
{
	std::vector<T1> rst(v1);
	--v1;
	return rst;
}

template <typename T>
std::vector<T> operator-(const std::vector<T>& v1)
{
	const SIZE_T N(v1.size());
	std::vector<T> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = -v1[j];
	return rst;
}

template <typename T1, typename T2>
auto operator*(const std::vector<T1>& v1, const std::vector<T2>& v2) -> std::vector<decltype(v1.at(0) * v2.at(0))>
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<decltype(v1.at(0) * v2.at(0))> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = v1[j] * v2[j];
	return rst;
}

template <typename T1, typename T2>
auto operator*(const std::vector<T1>& v1, const T2& a) -> std::vector<decltype(v1.at(0) * a)>
{
	const SIZE_T N(v1.size());
	std::vector<decltype(v1.at(0) * a)> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = v1[j] * a;
	return rst;
}

template <typename T1, typename T2>
auto operator*(const T2& a, const std::vector<T1>& v1) -> std::vector<decltype(a * v1.at(0))>
{
	return v1 * a;
}

template <typename T1, typename T2>
std::vector<T1>& operator*=(std::vector<T1>& v1, const T2& a)
{
	const SIZE_T N(v1.size());
	for(SIZE_T j(0); j < N; ++j)
		v1[j] *= a;
	return v1;
}

template <typename T1, typename T2>
std::vector<T1>& operator*=(std::vector<T1>& v1, const std::vector<T2>& v2)
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	for(SIZE_T j(0); j < N; ++j)
		v1[j] *= v2[j];
	return v1;
}

template <typename T1, typename T2>
auto operator/(const std::vector<T1>& v1, const std::vector<T2>& v2) -> std::vector<decltype(v1.at(0) / v2.at(0))>
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<decltype(v1.at(0) / v2.at(0))> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = v1[j] / v2[j];
	return rst;
}

template <typename T1, typename T2>
auto operator/(const std::vector<T1>& v1, const T2& a) -> std::vector<decltype(v1.at(0) / a)>
{
	return v1 * (static_cast<decltype(v1.at(0) / a)>(1.0) / a);
}

template <typename T1, typename T2>
auto operator/(const T2& a, const std::vector<T1>& v1) -> std::vector<decltype(a / v1.at(0))>
{
	const SIZE_T N(v1.size());
	std::vector<decltype(a / v1.at(0))> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = a / static_cast<decltype(rst.at(0))>(v1[j]);
	return rst;
}

template <typename T1, typename T2>
std::vector<T1>& operator/=(std::vector<T1>& v1, const T2& a)
{
	const auto a_inv(static_cast<decltype(v1.at(0) / a)>(1.0) / a);
	const SIZE_T N(v1.size());
	for(SIZE_T j(0); j < N; ++j)
		v1[j] *= a_inv;
	return v1;
}

template <typename T1, typename T2>
std::vector<T1>& operator/=(std::vector<T1>& v1, const std::vector<T2>& v2)
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	for(SIZE_T j(0); j < N; ++j)
		v1[j] /= v2[j];
	return v1;
}

template <typename T1, typename T2>
auto operator%(const std::vector<T1>& v1, const std::vector<T2>& v2) -> std::vector<decltype(v1.at(0) % v2.at(0))>
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	std::vector<decltype(v1.at(0) % v2.at(0))> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = v1[j] % v2[j];
	return rst;
}

template <typename T1, typename T2>
inline auto operator%(const std::vector<T1>& v1, const T2& a) -> std::vector<decltype(v1.at(0) % a)>
{
	const SIZE_T N(v1.size());
	std::vector<decltype(v1.at(0) % a)> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = v1[j] % a;
	return rst;
}

template <typename T1, typename T2>
auto operator%(const T2& a, const std::vector<T1>& v1) -> std::vector<decltype(a % v1.at(0))>
{
	const SIZE_T N(v1.size());
	std::vector<decltype(a % v1.at(0))> rst(N);
	for(SIZE_T j(0); j < N; ++j)
		rst[j] = a % v1[j];
	return rst;
}

template <typename T1, typename T2>
std::vector<T1>& operator%=(std::vector<T1>& v1, const T2& a)
{
	const SIZE_T N(v1.size());
	for(SIZE_T j(0); j < N; ++j)
		v1[j] %= a;
	return v1;
}

template <typename T1, typename T2>
std::vector<T1>& operator%=(std::vector<T1>& v1, const std::vector<T2>& v2)
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	for(SIZE_T j(0); j < N; ++j)
		v1[j] %= v2[j];
	return v1;
}

// relational operators

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator==(const std::vector<T1>& v1, const std::vector<T2>& v2)
{
	return equal(v1, v2);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator==(const std::vector<T1>& v1, const T2& a)
{
	return equal(a, v1);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator==(const T2& a, const std::vector<T1>& v1)
{
	return equal(v1, a);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator!=(const std::vector<T1>& v1, const std::vector<T2>& v2)
{
	return not_equal(v1, v2);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator!=(const std::vector<T1>& v1, const T2& a)
{
	return not_equal(v1, a);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator>=(const std::vector<T1>& v1, const std::vector<T2>& v2)
{
	return greater_equal(v1, v2);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator>=(const std::vector<T1>& v1, const T2& a)
{
	return greater_equal(v1, a);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator<=(const std::vector<T1>& v1, const std::vector<T2>& v2)
{
	return less_equal(v1, v2);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator<=(const std::vector<T1>& v1, const T2& a)
{
	return less_equal(v1, a);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator>(const std::vector<T1>& v1, const std::vector<T2>& v2)
{
	return greater(v1, v2);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator>(const std::vector<T1>& v1, const T2& a)
{
	return greater(v1, a);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator<(const std::vector<T1>& v1, const std::vector<T2>& v2)
{
	return less(v1, v2);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator<(const std::vector<T1>& v1, const T2& a)
{
	return less(v1, a);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator!=(const T2& a, const std::vector<T1>& v1)
{
	return not_equal(a, v1);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator>=(const T2& a, const std::vector<T1>& v1)
{
	return greater_equal(a, v1);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator<=(const T2& a, const std::vector<T1>& v1)
{
	return less_equal(a, v1);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator>(const T2& a, const std::vector<T1>& v1)
{
	return greater(a, v1);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator<(const T2& a, const std::vector<T1>& v1)
{
	return less(v1, a);
}

// logical operators

template <typename T>
inline std::vector<BOOL_T> operator!(const std::vector<T>& v) 
{
	return logical_not(v);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator&&(const std::vector<T1>& v1, const std::vector<T2>& v2) 
{
	// here overloading operator&& (as well as || below) shoule be fine since it does not return bool
	return logical_and(v1, v2);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator&&(const std::vector<T1>& v1, const T2& a) 
{
	return logical_and(v1, a);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator&&(const T2& a, const std::vector<T1>& v1) 
{
	return logical_and(a, v1);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator||(const std::vector<T1>& v1, const std::vector<T2>& v2) 
{
	return logical_or(v1, v2);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator||(const std::vector<T1>& v1, const T2& a) 
{
	return logical_or(v1, a);
}

template <typename T1, typename T2>
inline std::vector<BOOL_T> operator||(const T2& a, const std::vector<T1>& v1) 
{
	return logical_or(a, v1);
}

// other operators

template <typename T1, typename T2>
auto operator^(const std::vector<T1>& v1, const std::vector<T2>& v2) -> std::vector<decltype(std::pow(v1.at(0), v2.at(0)))>
{
	// pay attention to operator priority when using!
	return pow(v1, v2);
}

template <typename T1, typename T2>
auto operator^(const std::vector<T1>& v1, const T2& a) -> std::vector<decltype(std::pow(v1.at(0), a))>
{
	return pow(v1, a);
}

template <typename T1, typename T2>
auto operator^(const T2& a, const std::vector<T1>& v1) -> std::vector<decltype(std::pow(a, v1.at(0)))>
{
	return pow(a, v1);
}

template <typename T1, typename T2>
std::vector<T1>& operator^=(std::vector<T1>& v1, const std::vector<T2>& v2)
{
	assertsize(v1, v2);
	const SIZE_T N(v1.size());
	for(size_t j(0); j < N; ++j)
		v1[j] = std::pow(v1[j], v2[j]);
	return v1;
}

template <typename T1, typename T2>
std::vector<T1> operator^=(std::vector<T1>& v1, const T2& a)
{
	const SIZE_T N(v1.size());
	for(size_t j(0); j < N; ++j)
		v1[j] = std::pow(v1[j], a);
	return v1;
}

template <typename T1>
std::vector<T1> operator|(const std::vector<T1>& v1, const std::vector<T1>& v2)
{
	// connect vector elements, pay attention to operator priority when using!
	std::vector<T1> rst(v1);
	rst.insert(rst.end(), v2.begin(), v2.end());
	return rst;
}

template <typename T1, typename T2>
std::vector<T1> operator|(const std::vector<T1>& v1, const T2& a)
{
	std::vector<T1> rst;
	rst.reserve(v1.size() + 1);
	rst.insert(rst.end(), v1.begin(), v1.end());
	rst.push_back(static_cast<T1>(a));
	return rst;
}

template <typename T1, typename T2>
std::vector<T1> operator|(const T2& a, const std::vector<T1>& v1)
{
	std::vector<T1> rst(1, static_cast<T1>(a));
	rst.reserve(v1.size() + 1);
	rst.insert(rst.end(), v1.begin(), v1.end());
	return rst;
}

template <typename T1>
inline std::vector<T1>& operator|=(std::vector<T1>& v1, const std::vector<T1>& v2)
{
	v1.insert(v1.end(), v2.begin(), v2.end());
	return v1;
}

template <typename T1, typename T2>
inline std::vector<T1>& operator|=(std::vector<T1>& v1, const T2& a)
{
	v1.push_back(static_cast<T1>(a));
	return v1;
}

template <typename T1, typename T2>
typename std::enable_if< std::is_integral<T2>::value, std::vector<T1> >::type
operator<<(const std::vector<T1>& v1, const T2& a)
{
	// element rotation pay attention to operator priority when using!
	const SIZE_T N(v1.size());
	std::vector<T1> rst(v1);
	if (N >= 2) {
		if (a > 0) {
			std::rotate(rst.begin(), rst.begin() + a % N, rst.end());
		}
		else if (a < 0) {
			std::rotate(rst.begin(), rst.begin() + N - (-a) % N, rst.end());
		}
	}
	return rst;
}

template <typename T1, typename T2>
typename std::enable_if< std::is_integral<T2>::value, std::vector<T1>& >::type
operator<<=(std::vector<T1>& v1, const T2& a)
{
	const SIZE_T N(v1.size());
	if (N >= 2) {
		if (a > 0) {
			std::rotate(v1.begin(), v1.begin() + a % N, v1.end());
		}
		else if (a < 0) {
			std::rotate(v1.begin(), v1.begin() + N - (-a) % N, v1.end());
		}
	}
	return v1;
}

template <typename T1, typename T2>
inline typename std::enable_if< std::is_integral<T2>::value, std::vector<T1> >::type
operator>>(const std::vector<T1>& v1, const T2& a)
{
	return v1 << -a;
}

template <typename T1, typename T2>
inline typename std::enable_if< std::is_integral<T2>::value, std::vector<T1>& >::type
operator>>=(std::vector<T1>& v1, const T2& a)
{
	v1<<=-a;
	return v1;
}

template <typename T1>
auto operator~(const std::vector<T1>& v1) -> std::vector<decltype(v1.at(0) / norm(v1))>
{
	// extract unit vector along v1 direction
	// pay attention to operator priority when using!
	const auto v1nrm(norm(v1));
	const SIZE_T N(v1.size());
	std::vector<decltype(v1[0] / v1nrm)> rst(N);
	for (SIZE_T j(0); j < N; ++j) {
		rst[j] = v1[j] / v1nrm;
	}
	return rst;
}

template <typename T1, typename T2>
typename std::enable_if< std::is_floating_point<T1>::value and std::is_arithmetic<T2>::value, std::vector<T1> >::type
operator&(const std::vector<T1>& v, const std::vector<T2>& n)
{
	// vec&vec => get projection of v along the direction n
	// pay attention to operator priority when using!
	assertsize(v, n);

	const auto nnrm2(pow(norm(n), 2));

	misc::confirm( nnrm2 > 0.0, "operator&: n is zero vector!" );

	const SIZE_T N(v.size());
	const T1 vdotn_nnrm2(inner_product(v.begin(), v.end(), n.begin(), 0.0) / nnrm2);
	std::vector<T1> rst(N, 0.0);

	if (std::abs(vdotn_nnrm2) > 0.0) {
		for (SIZE_T j(0); j < N; ++j) {
			rst[j] = vdotn_nnrm2 * static_cast<T1>(n[j]); 
		}
	}
	return rst;
}

template <typename T1, typename T2>
typename std::enable_if< std::is_floating_point<T1>::value and std::is_arithmetic<T2>::value, std::vector<T1>& >::type
operator&=(std::vector<T1>& v, const std::vector<T2>& n)
{
	assertsize(v, n);

	const auto nnrm2(pow(norm(n), 2));

	misc::confirm( nnrm2 > 0.0, "operator&: n is zero vector!" );

	const SIZE_T N(v.size());
	const T1 vdotn_nnrm2(inner_product(v.begin(), v.end(), n.begin(), 0.0) / nnrm2);

	if (std::abs(vdotn_nnrm2) > 0.0) {
		for (SIZE_T j(0); j < N; ++j) {
			v[j] = vdotn_nnrm2 * static_cast<T1>(n[j]); 
		}
	}
	return v;
}

// other utilities

template <typename T>
inline std::vector<T> subvec(const std::vector<T>& v, SIZE_T start, SIZE_T N)
{
	// get sub-vector (v[start], v[start+1], ..., v[start + N - 1])
	return std::vector<T>(v.begin() + start, v.begin() + start + N);
}

template <typename T>
std::vector<T> subvec(const std::vector<T>& v, const std::vector<BOOL_T>& flags)
{
	// get sub-vector (v[k]) for all flags[k] == true
	assertsize(v, flags);
	const SIZE_T N(v.size());
	std::vector<T> rst;
	for (SIZE_T j(0); j < N; ++j) {
		if (flags[j]) {
			rst.push_back(v[j]);
		}
	}
	return rst;
}

template <typename T, typename T1>
typename std::enable_if< std::is_integral<T1>::value, std::vector<T> >::type
subvec(const std::vector<T>& v, const std::vector<T1>& idx)
{
	// get subvector (v[k]) for all k in std::vector idx
	const SIZE_T N(idx.size());
	std::vector<T> rst(N);
	for (T1 j(0); j < N; ++j) {
		rst[j] = v[idx[j]];
	}
	return rst;
}


template <typename T> 
std::vector<T> linspace(const T& xmin, const T& xmax, SIZE_T Nx, BOOL_T endpoint = true) 
{
	// see numpy linspace
	std::vector<T> x(Nx);
	const T dx(endpoint ? (xmax - xmin) / (Nx - 1) : (xmax - xmin) / Nx);
	for (SIZE_T j(0); j < Nx; ++j) {
		x[j] = xmin + j * dx;
	}
	return x;
}

template <typename T> 
std::vector<T> linspace(const T& xmin, const T& xmax, SIZE_T Nx, T& dx, BOOL_T endpoint = true) 
{
	// see numpy linspace, dx requested 
	std::vector<T> x(Nx);
	dx = (endpoint ? (xmax - xmin) / (Nx - 1) : (xmax - xmin) / Nx);
	for (SIZE_T j(0); j < Nx; ++j) {
		x[j] = xmin + j * dx;
	}
	return x;
}


template <typename T>
std::vector<T> arange(T start, T end, T step = 1) 
{
	// see numpy.arange
	std::vector<T> rst;
	rst.reserve(static_cast<SIZE_T>((end - start) / step + 1));

	T x(start);
	while (x < end) {
		rst.push_back(x);
		x += step;
	}
	return rst;
}

template <typename T>
inline std::vector<T> arange(T end) 
{
	// see numpy.arange
	return arange(static_cast<T>(0), end);
}


template <typename T>
void meshgrid(const std::vector<T>& X, const std::vector<T>& Y, std::vector<T>& meshX, std::vector<T>& meshY) 
{
	// see numpy.meshgrid
	const SIZE_T Nx(X.size());
	const SIZE_T Ny(Y.size());

	meshX.resize(Nx * Ny);
	meshY.resize(Nx * Ny);

	for (SIZE_T j(0); j < Nx; ++j) {
		for (SIZE_T k(0); k < Ny; ++k) {
			meshX[k + j * Ny] = X[j];
			meshY[k + j * Ny] = Y[k];
		}
	}
}


template <typename T1, typename T2>
typename std::enable_if< std::is_floating_point<T1>::value and std::is_arithmetic<T2>::value, std::vector<T1> >::type
inline component(const std::vector<T1>& v, const std::vector<T2>& n)
{
	// extract projection of v on n direction
	return v&n;
}

template <typename T>
inline auto unit(const std::vector<T>& v) -> std::vector<decltype(v.at(0) / norm(v))>
{
	// return unit vector on v direction
	return ~v; 
}

template <typename T>
std::vector<T> outer_product(const std::vector<T>& v1, const std::vector<T>& v2)
{
	// get outer product of v1 and v2
	const SIZE_T Nrow(v1.size());
	const SIZE_T Ncol(v2.size());
	std::vector<T> rst(Nrow * Ncol);
	for (SIZE_T r(0); r < Nrow; ++r) {
        for (SIZE_T c(0); c < Ncol; ++c) {
            rst[r+c*Nrow] = v1[r] * v2[c];
        }
	}
	return rst;
}

#endif
