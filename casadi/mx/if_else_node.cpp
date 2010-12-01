#include "if_else_node.hpp"
#include "../stl_vector_tools.hpp"
#include <cassert>
#include <vector>
#include <sstream>

using namespace std;

namespace CasADi{

double IfElseNode::tol = 1e-6;
  
IfElseNode::IfElseNode(const MX& cond, const MX& if_true, const MX& if_false) : MXNode(cond,if_true,if_false){
  assert(cond.numel()==1);
  assert(if_true.size1()==if_false.size1() || if_true.size2()==if_false.size2());  
  sz = if_true.size();
}

IfElseNode* IfElseNode::clone() const{
  return new IfElseNode(*this);
}

void IfElseNode::print(std::ostream &stream) const{
  stream << "(" << dep(0) << "?" <<  dep(1) << ":" << dep(2) << ")";
}

void IfElseNode::evaluate(int fsens_order, int asens_order){
  assert(fsens_order==0 || asens_order==0);
  
if(fsens_order==0){
  
  const vector<double>& cond = dep(0)->val(0);  // condition
  const vector<double>& if_true = dep(1)->val(0);  // if condition true
  const vector<double>& if_false = dep(2)->val(0);  // if condition false
  vector<double>& res = val(0);

  res = fabs(cond[0])>tol ? if_true : if_false;
} else {
  const vector<double>& cond = dep(0)->val(0);  // condition
  const vector<double>& if_true_der = dep(1)->val(1);  // if condition true derivative
  const vector<double>& if_false_der = dep(2)->val(1);  // if condition false derivative
  vector<double>& res_der = val(1);

  res_der = fabs(cond[0])>tol ? if_true_der : if_false_der;
}

if(asens_order>0){
  const vector<double>& cond = dep(0)->val(0);  // condition
  vector<double>& if_true_der = dep(1)->val(1);  // if condition true derivative
  vector<double>& if_false_der = dep(2)->val(1);  // if condition false derivative
  const vector<double>& res_der = val(1);

  if(fabs(cond[0])>tol){ // if true
    for(int i=0; i<res_der.size(); ++i)
      if_true_der[i] += res_der[i];
  } else { // if false
    for(int i=0; i<res_der.size(); ++i)
      if_false_der[i] += res_der[i];
  }
}
}

} // namespace CasADi

