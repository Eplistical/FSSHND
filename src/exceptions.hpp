#ifndef _EXCEPTIONS_HPP
#define _EXCEPTIONS_HPP
// module for exception handling

#include <string>
#include <stdexcept>

namespace {


    // --- Exceptions Definition --- //


    class Error : public std::exception
    {
        protected:
            std::string m_title;
            std::string m_msg;
        public:
            explicit Error(const std::string& msg = "") {
                m_title = "Error"; 
                m_msg = (msg == "") ? m_title : (m_title + ":" + msg);
            }
            // default destructor
            virtual ~Error() = default;
            virtual const char* what() const noexcept{ return m_msg.c_str(); }
    };

    class IndexError final : public Error
    {
        public:
            explicit IndexError(const std::string& msg = "")
            {
                m_title = "IndexError"; 
                m_msg = (msg == "") ? m_title : (m_title + ":" + msg);
            } 
    };


    // --- Handeler --- //


    template <typename ExceptionType = Error>
        void confirm(bool condition, const std::string& msg) {
            if (not condition) {
                throw(ExceptionType(msg));
            }
        }


    template <typename ExceptionType = Error>
        void crash(bool condition, const std::string& msg) {
            confirm(not condition, msg);
        }

};
#endif
