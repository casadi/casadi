/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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
    
  //! \brief Copy constructor
  CasadiException(const CasadiException& ex){
    errbuf << ex.errbuf.str();
  }
  
  //! \brief Form message string    
  explicit CasadiException(const std::string& msg){
    errbuf << msg;
  }

  //! \brief Destructor
  ~CasadiException() throw(){}
    
  //! \brief Display error
  virtual const char* what() const throw(){
     return errbuf.str().c_str();
  }
  
  //! \brief Append a message
  CasadiException& operator<<(const std::string& msg){
    errbuf << msg;
    return *this;
  }

  //! \brief Append an exception
  CasadiException& operator<<(const std::exception& ex){
    errbuf << " => " << ex.what();
    return *this;
  }

  protected:
  std::stringstream errbuf;
};



} // namespace CasADi

#endif // CASADI_EXCEPTION_HPP
