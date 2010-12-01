#ifndef MATRIX_MATRIX_OP_HPP
#define MATRIX_MATRIX_OP_HPP

#include "mx_node.hpp"

namespace CasADi{
/** \brief Represents any binary operation that involves two matrices 
  \author Joel Andersson 
  \date 2010	
*/
class MatrixMatrixOp : public MXNode{
public:

/** \brief  Constructor */
MatrixMatrixOp (OPERATION op, const MX& x, const MX& y);

/** \brief  Clone function */
virtual MatrixMatrixOp * clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
//  virtual void evaluateAdj();
  
protected:
  //! \brief Operation
  OPERATION op;
};

} // namespace CasADi


#endif // MATRIX_MATRIX_OP_HPP
