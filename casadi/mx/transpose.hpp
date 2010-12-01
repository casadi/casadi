#ifndef TRANSPOSE_HPP
#define TRANSPOSE_HPP

#include "mx_node.hpp"

namespace CasADi{

/** Represents a transposition of an MX
  \author Joel Andersson 
  \date 2010
*/
class Transpose : public MXNode{
friend class MX;

public:

/** \brief  Constructor */
Transpose(const MX& x);

/** \brief  Clone function */
virtual Transpose* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
//  virtual void evaluateAdj();

};

} // namespace CasADi

#endif // TRANSPOSE_HPP
