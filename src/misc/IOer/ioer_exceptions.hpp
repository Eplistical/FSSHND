#ifndef _IOER_EXCEPTIONS_HPP
#define _IOER_EXCEPTIONS_HPP
// module for ioer exceptions

#include <string>
#include <stdexcept>

namespace ioer
{
    // Exceptions
    class IOError : public std::exception
    {
        protected:
            std::string _title;
            std::string _msg;
        public:
            explicit IOError(const std::string& msg = "") :
                _title("IOError"), 
                _msg((msg == "")?_title:(_title + ":" + msg)) {}
            // default destructor
            virtual ~IOError() = default;
            virtual const char* what() const noexcept{ return _msg.c_str(); }
    };

    class FileOpenedError : public IOError
    {
        public:
            explicit FileOpenedError(const std::string& msg = "")
            {
                _title = "FileOpenedError"; 
                _msg = (msg == "")?_title:(_title + ":" + msg);
            } 
            ~FileOpenedError() = default;
    };

    class FileNotOpenedError : public IOError
    {
        public:
            explicit FileNotOpenedError(const std::string& msg = "")
            {
                _title = "FileNotOpenedError"; 
                _msg = (msg == "")?_title:(_title + ":" + msg);
            } 
            ~FileNotOpenedError() = default;
    };

    class FileNotExistError : public IOError
    {
        public:
            explicit FileNotExistError(const std::string& msg = "")
            {
                _title = "FileNotExistError"; 
                _msg = (msg == "")?_title:(_title + ":" + msg);
            } 
            ~FileNotExistError() = default;
    };

    class FileModeIncompatibleError : public IOError
    {
        public:
            explicit FileModeIncompatibleError(const std::string& msg = "")
            {
                _title = "FileModeIncompatibleError"; 
                _msg = (msg == "")?_title:(_title + ":" + msg);
            } 
            ~FileModeIncompatibleError() = default;
    };

};
#endif
