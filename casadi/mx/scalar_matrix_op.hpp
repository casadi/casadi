#ifndef SCALAR_MATRIX_OP_HPP
#define SCALAR_MATRIX_OP_HPP

#include "mx_node.hpp"

namespace CasADi{
/** Represents any binary operation that involves a matrix and a scalar
  \author Joel Andersson 
  \date 2010
*/
class ScalarMatrixOp : public MXNode{
public:

/** \brief  Constructor */
ScalarMatrixOp (OPERATION op, const MX& x, const MX& y);

/** \brief  Clone function */
virtual ScalarMatrixOp * clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
/*  virtual void evaluateAdj();*/
  

  //! \brief Operation
  OPERATION op;
};

} // namespace CasADi


#endif // SCALAR_MATRIX_OP_HPP
