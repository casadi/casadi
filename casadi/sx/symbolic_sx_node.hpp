#ifndef SYMBOLIC_SCALAR_HPP
#define SYMBOLIC_SCALAR_HPP

#include "sx_node.hpp"

namespace CasADi{

/** \brief Represents a scalar symbolic expression
  \author Joel Andersson 
  \date 2010
  A regular user is not supposed to work with this Node class.
  This user can call SX(name) instead.
*/
class SymbolicSXNode : public SXNode{
public:

explicit SymbolicSXNode(const std::string &name) : name(name){}

virtual bool isSymbolic() const{ 
  return true; 
}

virtual const std::string& getName() const{ 
  return name; 
}

/** \brief  Name */
std::string name;

protected:

/** \brief  print */
virtual void print(std::ostream &stream) const{
  stream << name;
}

};

} // namespace CasADi

#endif // SYMBOLIC_SCALAR_HPP
