#ifndef UNARY_OP_HPP
#define UNARY_OP_HPP

#include "mx_node.hpp"

namespace CasADi{
/** \brief Represents a general unary operation on an MX
  \author Joel Andersson 
  \date 2010
*/	
class UnaryOp : public MXNode{
public:

/** \brief  Constructor */
UnaryOp(OPERATION op, const MX& x);

/** \brief  Clone function */
virtual UnaryOp * clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
//  virtual void evaluateAdj();
  
protected:
  //! \brief operation
  OPERATION op;
};

} // namespace CasADi


#endif // UNARY_OP_HPP
