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
#ifdef NDEBUG
// Release mode
#define casadi_assert(x)
#define casadi_assert_message(x,msg)
#else // NDEBUG
// Debug mode
// Convert to string
#define CASADI_ASSERT_STR1(x) #x
#define CASADI_ASSERT_STR(x) CASADI_ASSERT_STR1(x)
// This assersion if for errors caused by bugs in CasADi
#define casadi_assert(x) \
if(!(x)) throw CasadiException("The assertion " CASADI_ASSERT_STR(x) " on line " CASADI_ASSERT_STR(__LINE__) " of file " CASADI_ASSERT_STR(__FILE__) " failed. Please notify the CasADi developers.")

// This assersion if for illigal user inputs that should not be checked in the release version for effiency resonds, for example out of bounds
#define casadi_assert_message(x,msg) \
if(!(x)) throw CasadiException("The assertion " CASADI_ASSERT_STR(x) " on line " CASADI_ASSERT_STR(__LINE__) " of file " CASADI_ASSERT_STR(__FILE__) " failed. " msg)

#endif // NDEBUG
  
} // namespace CasADi

#endif // CASADI_EXCEPTION_HPP
