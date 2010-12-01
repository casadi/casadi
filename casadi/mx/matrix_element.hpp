#ifndef MATRIX_ELEMENT_HPP
#define MATRIX_ELEMENT_HPP

#include "mx_node.hpp"

namespace CasADi{
/** \brief  Element of the matrix which is allowed to change the object (matrix style index)
  \author Joel Andersson 
  \date 2010
  \see Element	
*/
class MatrixElement : public MXNode{
public:

/** \brief  Constructor */
MatrixElement(const MX& x, int i, int j);

/** \brief  Clone function */
virtual MatrixElement* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
//  virtual void evaluateAdj();

protected:
int i, j;

};

} // namespace CasADi


#endif // MATRIX_ELEMENT_HPP
