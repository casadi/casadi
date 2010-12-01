#ifndef SX_NODE_HPP
#define SX_NODE_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

/** \brief  Scalar expression (which also works as a smart pointer class to this class) */
#include "sx.hpp"

namespace CasADi{

/** \brief  Interal node class for SX
  \author Joel Andersson 
  \date 2010
*/
class SXNode{
friend class SX;

public:

/** \brief  constructor */
SXNode();

/** \brief  destructor  */
~SXNode(); // non-virtual to avoid stack overflow upon destruction!
//@{
/** \brief  check properties of a node */
virtual bool isConstant() const; // check if constant
virtual bool isInteger() const; // check if integer
virtual bool isSymbolic() const; // check if symbolic
virtual bool isBinary() const; // check if binary
virtual bool isZero() const; // check if zero
virtual bool isOne() const; // check if one
virtual bool isMinusOne() const; // check if minus one
virtual bool isNan() const; // check if not a number
virtual bool isInf() const; // check if infinity
virtual bool isMinusInf() const; // check if minus infinity
//@}

//@{
/** \brief  Get value of a constant node */
virtual double getValue() const;  // only works for constant nodes
virtual int getIntValue() const;  // only works for integer nodes
//@}

virtual const std::string& getName() const; // get the name
/** \brief get the operation 
only for binary nodes
*/
virtual int getOp() const; // get the operation (only for binary nodes)
/// comparison
bool isEqual(const SXNode& node) const; // comparison
/// comparison
bool isEqual(const SX& scalar) const; // comparison

/** \brief  get the reference of a child */
virtual const SX& dependent(int i) const;

/** \brief  Check if smooth */
virtual bool isSmooth() const;

/** \brief  print */
virtual void print(std::ostream &stream) const = 0;

/** Temporary variables to be used in user algorithms like sorting, 
 the user is resposible of making sure that use is thread-safe
 The variable is initialized to zero
*/
int temp;
// int temp2;

protected:
// Reference counter -- counts the number of parents of the node
unsigned int count;

};

} // namespace CasADi


/** \brief  Derived classes */
#include "constant_sx_node.hpp"
#include "symbolic_sx_node.hpp"
#include "binary_sx_node.hpp"

#endif // SX_NODE_HPP
