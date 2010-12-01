#ifndef MATRIX_SCALAR_OP_HPP
#define MATRIX_SCALAR_OP_HPP

#include "mx_node.hpp"

namespace CasADi{
/** \brief Represents a general matrix scalar opertion on MXes
  \author Joel Andersson 
  \date 2010
*/	
class MatrixScalarOp : public MXNode{
public:

/** \brief  Constructor */
MatrixScalarOp (OPERATION op, const MX& x, const MX& y);

/** \brief  Clone function */
virtual MatrixScalarOp * clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
//  virtual void evaluateAdj();
  
protected:

  OPERATION op;
};

} // namespace CasADi


#endif // MATRIX_SCALAR_OP_HPP
