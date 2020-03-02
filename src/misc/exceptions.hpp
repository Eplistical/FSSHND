#ifndef _EXCEPTIONS_HPP
#define _EXCEPTIONS_HPP
// module for exception handling

#include <string>
#include <stdexcept>
#include "types.hpp"

namespace misc {
    

        // --- Exceptions Definitions --- //


        class Error : public std::exception {
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

        class IndexError final : public Error {
            public:
                explicit IndexError(const std::string& msg = "")
                {
                    m_title = "IndexError";
                    m_msg = (msg == "") ? m_title : (m_title + ":" + msg);
                }
        };

        class ValueError final : public Error {
            public:
                explicit ValueError(const std::string& msg = "")
                {
                    m_title = "ValueError";
                    m_msg = (msg == "") ? m_title : (m_title + ":" + msg);
                }
        };

        class StateError final : public Error {
            public:
                explicit StateError(const std::string& msg = "")
                {
                    m_title = "StateError";
                    m_msg = (msg == "") ? m_title : (m_title + ":" + msg);
                }
        };
        
        class RuntimeError final : public Error {
            public:
                explicit RuntimeError(const std::string& msg = "")
                {
                    m_title = "RuntimeError";
                    m_msg = (msg == "") ? m_title : (m_title + ":" + msg);
                }
        };

};
#endif
