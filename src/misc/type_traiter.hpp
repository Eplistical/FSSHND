#ifndef _TYPE_TRAITER_HPP
#define _TYPE_TRAITER_HPP
/* module for TypeTraiter
 *
 */
#include <complex>
#include <string>
#include <array>
#include <vector>
#include <deque>
#include <list>
#include <forward_list>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <typeindex>
#include <type_traits>


namespace type_traiter
{
    using namespace std;

    // -- is_bool -- //
    template<typename T>
        struct is_bool 
        {
            static constexpr bool value 
                = is_same<typename remove_cv<T>::type, bool>::value; 
        };

    // -- is_string -- //
    template<typename ... Types>
        struct is_string 
        {
            static constexpr bool value = false;
        };

    template<typename ... Types>
        struct is_string<basic_string<Types...> > 
        {
            static constexpr bool value = true;
        };

    // -- is_c_string -- //
    template<typename ... Types>
        struct is_c_string 
        {
            static constexpr bool value = false;
            static constexpr size_t size = 0;
        };

    template<typename T, size_t N>
        struct is_c_string<T[N]> 
        {
            // char[N]
            static constexpr bool value = is_same<char, typename remove_cv<T>::type>::value;
            static constexpr size_t size = N;
        };

    template<typename T>
        struct is_c_string<T> 
        {
            // char*
            static constexpr bool value = is_same<char*, typename remove_cv<T>::type>::value ||
                is_same<const char*, typename remove_cv<T>::type>::value;
            static constexpr size_t size = 0;
        };

    // -- is_complex -- //
    template<typename T, typename ... Types>
        struct is_complex 
        {
            static constexpr bool value = false;
        };

    template<typename ... Types>
        struct is_complex<complex<Types...> > 
        {
            static constexpr bool value = true;
        };

    // -- is_array_container -- //
    template<typename T, typename ... Types>
        struct is_array_container 
        {
            static constexpr bool value = false;
        };

    template<typename ... Types>
        struct is_array_container<std::array<Types...> > 
        {
            static constexpr bool value = true;
        };

    // -- is_direct_outputable -- //
    template <typename T>
        struct is_direct_outputable 
        {
            static constexpr bool value
                =  is_fundamental<T>::value
                || is_string<T>::value
                || is_c_string<T>::value
                || is_complex<T>::value
                || is_enum<T>::value
                || is_same<char*, typename remove_cv<T>::type>::value;
        };

    // -- is_deque -- //
    template<typename T, typename ... Types>
        struct is_deque 
        {
            static constexpr bool value = false;
        };

    template<typename ... Types>
        struct is_deque<deque<Types...> > 
        {
            static constexpr bool value = true;
        };


    // -- is_forward_list -- // 
    template<typename T, typename ... Types>
        struct is_forward_list 
        {
            static constexpr bool value = false;
        };

    template<typename ... Types>
        struct is_forward_list<forward_list<Types...> > 
        {
            static constexpr bool value = true;
        };


    // -- is_list -- // 
    template<typename T, typename ... Types>
        struct is_list 
        {
            static constexpr bool value = false;
        };

    template<typename ... Types>
        struct is_list<list<Types...> > 
        {
            static constexpr bool value = true;
        };


    // -- is_vector -- // 
    template<typename T, typename ... Types>
        struct is_vector 
        {
            static constexpr bool value = false;
        };

    template<typename ... Types>
        struct is_vector<vector<Types...> > 
        {
            static constexpr bool value = true;
        };


    // -- is_map -- // 
    template<typename T, typename ... Types>
        struct is_map 
        {
            static constexpr bool value = false;
        };

    template<typename ... Types>
        struct is_map<map<Types...> > 
        {
            static constexpr bool value = true;
        };


    // -- is_set -- // 
    template<typename T, typename ... Types>
        struct is_set {
            static constexpr bool value = false;
        };

    template<typename ... Types>
        struct is_set<set<Types...> > 
        {
            static constexpr bool value = true;
        };


    // -- is_unordered_map -- // 
    template<typename T, typename ... Types>
        struct is_unordered_map 
        {
            static constexpr bool value = false;
        };

    template<typename ... Types>
        struct is_unordered_map<unordered_map<Types...> > 
        {
            static constexpr bool value = true;
        };

    // -- is_unordered_set -- // 
    template<typename T, typename ... Types>
        struct is_unordered_set 
        {
            static constexpr bool value = false;
        };

    template<typename ... Types>
        struct is_unordered_set<unordered_set<Types...> > 
        {
            static constexpr bool value = true;
        };


    // -- is_sequence_container -- // 
    template <typename T>
        struct is_sequence_container 
        {
            static constexpr bool value
                =  is_deque<T>::value
                || is_forward_list<T>::value
                || is_list<T>::value
                || is_vector<T>::value
                || is_array_container<T>::value;
        };


    // -- is_associative_container -- // 
    template <typename T>
        struct is_associative_container 
        {
            static constexpr bool value
                =  is_map<T>::value
                || is_set<T>::value;
        };


    // -- is_unordered_asscociative_container -- // 
    template <typename T>
        struct is_unordered_associative_container 
        {
            static constexpr bool value
                =  is_unordered_map<T>::value
                || is_unordered_set<T>::value;
        };


    // -- is_container -- // 
    template <typename T>
        struct is_container 
        {
            static constexpr bool value
                =  is_sequence_container<T>::value
                || is_associative_container<T>::value
                || is_unordered_associative_container<T>::value;
        };
};

#endif
