#ifndef IF_ELSE_NODE_HPP
#define IF_ELSE_NODE_HPP

#include "mx_node.hpp"

namespace CasADi{

/** \brief Represents a branch in an MX tree
  \author Joel Andersson 
  \date 2010
*/
class IfElseNode : public MXNode{
public:

/** \brief  Constructors */
IfElseNode(const MX& cond, const MX& if_true, const MX& if_false);
//IfElseNode(const MX& cond, const FX& fcn_true, const FX& fcn_false, const std::vector<MX>& arg);

/** \brief  Clone function */
virtual IfElseNode* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
//  virtual void evaluateAdj();
  
protected:

  // Functions to evaluate if true/false
//  FX fcn_true, fcn_false;
  
  //! \brief tolerance
  static double tol;
};

} // namespace CasADi


#endif // IF_ELSE_NODE_HPP
