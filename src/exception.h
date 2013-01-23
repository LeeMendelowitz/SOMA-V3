#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <exception>
#include <string>

// Generic Exception class
class Exception : public std::exception
{ 
    public:
    Exception(const std::string& msg) throw(): msg_(msg) {};
    ~Exception() throw() {};

    virtual const char * what() const throw()
    {
        return msg_.c_str();
    }

    private:
        std::string msg_;
};

#endif
