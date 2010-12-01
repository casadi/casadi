#ifndef BINARY_SCALAR_HPP
#define BINARY_SCALAR_HPP

#include "sx_node.hpp"
#include "binary_functions.hpp"

namespace CasADi{

/** \brief Represents a basic binary operation on two SX nodes
  \author Joel Andersson 
  \date 2010
*/
class BinarySXNode : public SXNode{
public:

/** \brief  Constructors */
BinarySXNode(OPERATION op_, const SX& child1_){
 op = op_;
 child[0] = child1_;
 child[1] = 0;
}

BinarySXNode(OPERATION op_, const SX& child1_, const SX& child2_){
 op = op_;
 child[0] = child1_;
 child[1] = child2_;
}

virtual bool isSmooth() const;

virtual bool isBinary() const{ return true; }

/** \brief  get the reference of a child */
virtual const SX& dependent(int i) const;

/** \brief  Get the operation */
virtual int getOp() const{ return op;}

/** \brief  Data members */
OPERATION  op;
SX      child[2];

virtual void print(std::ostream &stream) const;


};

} // namespace CasADi


#endif // BINARY_SCALAR_HPP
