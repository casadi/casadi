#ifndef SYMBOLIC_MATRIX_HPP
#define SYMBOLIC_MATRIX_HPP

#include "mx_node.hpp"
#include "../sx/sx_matrix.hpp"

namespace CasADi{
/** \brief Represents a symbolic MX
  \author Joel Andersson 
  \date 2010
  A regular user is not supposed to work with this Node class.
  This user can call MX(name,n,m) directly.
*/
class SymbolicMatrix : public MXNode{
public:

/** \brief  Constructors */
explicit SymbolicMatrix(const std::string& name, int n=1, int m=1);

/** \brief  Clone function */
virtual SymbolicMatrix* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
//  virtual void evaluateAdj();

/** \brief  Is symbolic */
  virtual bool isSymbolic() const;
  
/** \brief  Get the name */
  virtual const std::string& getName() const;
  
protected:
  std::string name;
};

} // namespace CasADi


#endif // SYMBOLIC_MATRIX_HPP
