#ifndef PRINTABLE_OBJECT_HPP
#define PRINTABLE_OBJECT_HPP

#include <iostream>
#include <string>

namespace CasADi{

/** \brief Base class for objects that have a natural string representation
  \author Joel Andersson 
  \date 2010
*/
class PrintableObject{
  public:
    /// Print a destription of the object
    virtual void print(std::ostream &stream=std::cout) const;

    /// Print a representation of the object
    virtual void repr(std::ostream &stream) const;

    /// Print a representation of the object to a stream
    friend std::ostream& operator<<(std::ostream &stream, const PrintableObject& obj);
    
    /// Return a string with a representation (for SWIG)
    std::string getRepresentation() const;

    /// Return a string with a destription (for SWIG)
    std::string getDescription() const;
};

} // namespace CasADi


#endif // PRINTABLE_OBJECT_HPP

