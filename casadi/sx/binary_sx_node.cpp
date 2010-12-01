#include "binary_sx_node.hpp"
#include <cassert>

using namespace std;
namespace CasADi{

void BinarySXNode::print(ostream &stream) const{
  stringstream s0,s1;
  s0 << child[0];
  s1 << child[1];
  print_c[op](stream,s0.str(),s1.str());
}

bool BinarySXNode::isSmooth() const{
  if(op == STEP_NODE || op == FLOOR_NODE)
    return false;
  else
    return true;
}

const SX& BinarySXNode::dependent(int i) const{
  assert(i==0 || i==1);
  return child[i];
}

} // namespace CasADi
