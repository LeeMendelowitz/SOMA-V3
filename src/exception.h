#ifndef EXCEPTION_H
#define EXCEPTION_H
#include <exception>

using namespace std;

// Generic Exception class
class Exception : public exception
{ 
    public:
    Exception(const string& msg) throw(): msg_(msg) {};
    ~Exception() throw() {};

    virtual const char * what() const throw()
    {
        return msg_.c_str();
    }

    private:
        string msg_;
};

#endif
