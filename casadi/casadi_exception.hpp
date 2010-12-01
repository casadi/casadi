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
