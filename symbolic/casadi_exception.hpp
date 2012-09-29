/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef CASADI_EXCEPTION_HPP
#define CASADI_EXCEPTION_HPP

#include <exception>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>

namespace CasADi{

/** \brief  Casadi exception class
	\author Joel Andersson 
	\date 2010
	Example for simple exception throwing:
	\code
		throw CasadiException("This is a nasty error");
	\endcode
	Example for exception chaining:
	\code
		try {
			throw CasadiException("This is a nasty error");
		catch (CasadiException &e) {
			throw CasadiException("Serious error.") << e;
		}
	\endcode
*/
class CasadiException : public std::exception{
  public:
  //! \brief Default constructor
  CasadiException(){
  }
    
  //! \brief Form message string    
  explicit CasadiException(const std::string& msg) : msg_(msg){}

  //! \brief Destructor
  ~CasadiException() throw(){}
    
  //! \brief Display error
  virtual const char* what() const throw(){
    return msg_.c_str();
  }
  
  //! \brief Append a message
  CasadiException& operator<<(const std::string& msg){
    msg_ += msg;
    return *this;
  }

  //! \brief Append an exception
  CasadiException& operator<<(const std::exception& ex){
    msg_ += " => ";
    msg_ += ex.what();
    return *this;
  }

  protected:
  std::string msg_;
};


// Assertion similar to the standard C assert statement, with the difference that it throws an exception with the same information
#ifdef CASADI_NDEBUG
// Release mode
#define casadi_assert(x)
#define casadi_assert_message(x,msg)
#define casadi_assert_warning(x,msg)
#define casadi_warning(msg)
#define casadi_error(msg)
 
#else // CASADI_NDEBUG
// Debug mode
// Convert to string
#define CASADI_ASSERT_STR1(x) #x
#define CASADI_ASSERT_STR(x) CASADI_ASSERT_STR1(x)

// String denoting where the assertation is situated
#define CASADI_ASSERT_WHERE " on line " CASADI_ASSERT_STR(__LINE__) " of file " CASADI_ASSERT_STR(__FILE__)

#define casadi_log(msg) \
  if(verbose()){ \
    std::cout << "CasADi log message: " << msg << std::endl; \
  }
  
#define casadi_error(msg) \
 {\
  std::stringstream ss_internal_; \
  ss_internal_ << CASADI_ASSERT_WHERE << std::endl << msg  <<  std::endl; \
  throw CasADi::CasadiException(ss_internal_.str()); \
 }

// This assertion checks for illegal user inputs. It will not be checked if CASADI_NDEBUG is defined
#define casadi_assert_message(x,msg) \
{ \
  bool is_ok; \
  try{ \
    is_ok = x; \
  } catch(std::exception& ex){ \
      throw CasADi::CasadiException(std::string("When trying to check the assertion \"" CASADI_ASSERT_STR(x) "\"" CASADI_ASSERT_WHERE ", caught: \n")+ex.what());\
  } \
 if(!is_ok) { \
  std::stringstream ss_internal_; \
  ss_internal_ << "The assertion \"" CASADI_ASSERT_STR(x) "\"" CASADI_ASSERT_WHERE " failed. " << std::endl << msg  <<  std::endl; \
  throw CasADi::CasadiException(ss_internal_.str()); \
 }\
} \

// This assersion if for errors caused by bugs in CasADi, use it instead of C:s assert(), but never in destructors
#define casadi_assert(x) casadi_assert_message(x,"(Hint for developers: CasadiOptions.setCatchErrorsPython(False) to obtain gdb stacktrace in python.)" << std::endl << "Please notify the CasADi developers.")

// This is for warnings to be issued when casadi is not in release mode and an assertion fails
#define casadi_assert_warning(x,msg) \
if((x)==false){ \
  std::cerr << "CasADi warning: \"" << msg << "\" (assertion \"" CASADI_ASSERT_STR(x) "\"" CASADI_ASSERT_WHERE " failed.)" << std::endl;\
}

// This is for warnings to be issued when casadi is not in release mode
#define casadi_warning(msg) \
std::cerr << "CasADi warning: \"" << msg << "\" issued " CASADI_ASSERT_WHERE ". " << std::endl;

#endif // CASADI_NDEBUG
  
} // namespace CasADi

#endif // CASADI_EXCEPTION_HPP
