#include "sx_node.hpp"
#include <limits>
#include <typeinfo>
#include <cassert>

using namespace std;
namespace CasADi{

SXNode::SXNode(){
  count = 0;
  temp = 0;
/*  temp2 = 0;*/
}

SXNode::~SXNode(){
  // Make sure that this is there are no scalar expressions pointing to it when it is destroyed
  assert(count==0);
}

double SXNode::getValue() const{
  return numeric_limits<double>::quiet_NaN();
/*  std::cerr << "getValue() not defined for class " << typeid(*this).name() << std::endl;
  throw "SXNode::getValue()";*/
}

int SXNode::getOp() const{
  throw CasadiException(string("getOp() not defined for class ") + typeid(*this).name());
}

int SXNode::getIntValue() const{
  throw CasadiException(string("getIntValue() not defined for class ") + typeid(*this).name());
}

bool SXNode::isZero() const{
  return false;
}

bool SXNode::isOne() const{
  return false;
}

bool SXNode::isMinusOne() const{
  return false;
}

bool SXNode::isNan() const{
  return false;
}

bool SXNode::isInf() const{
  return false;
}

bool SXNode::isMinusInf() const{
  return false;
}

bool SXNode::isConstant() const{
  return false;
}

bool SXNode::isInteger() const{
  return false;
}

bool SXNode::isSymbolic() const{
  return false;
}

bool SXNode::isBinary() const{
  return false;
}

bool SXNode::isEqual(const SXNode& node) const{
  // check if the objects are the same or if both are integers with the same values
  return (this == &node) || (isConstant() && node.isConstant() && getValue() == node.getValue());
}

bool SXNode::isEqual(const SX& scalar) const{
  return isEqual(*scalar.get());
}

const std::string& SXNode::getName() const{
  throw CasadiException("SXNode::getName failed, the node must be symbolic");
}

const SX& SXNode::dependent(int i) const{
  stringstream ss;
  ss << "child() not defined for class " << typeid(*this).name() << std::endl;
  throw CasadiException(ss.str());
}


bool SXNode::isSmooth() const{
  return true; // nodes are smooth by default
}

} // namespace CasADi
