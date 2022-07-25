/***********************************************/
/**
* @file exception.h
*
* @brief Exception handling.
*
* @author Torsten Mayer-Guerr
* @date 2001-06-16
*
*/
/***********************************************/

#ifndef __GROOPS_EXCEPTION__
#define __GROOPS_EXCEPTION__

/***** CLASS ***********************************/

/**
* @brief Exception handling.
*/
class Exception : public std::exception
{
  std::string message;

public:

/**
* Exception handling.
* @param msg error message
*/
Exception(const std::string &msg) throw() : message(msg) {}
virtual ~Exception() throw()  {}

/**
* Returns a C-style character string
* describing the general cause of the current error.
*/
virtual const char* what() const throw() {return message.c_str();}
};

/***********************************************/

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)
#define ERRORLINE (std::string("in ")+__FILE__+":"+QUOTEME(__LINE__)+" ("+__func__+")")
#define rethrow(e) throw(Exception(ERRORLINE+"\n"+e.what()))

/***********************************************/

#endif /* __GROOPS_Exception__ */
