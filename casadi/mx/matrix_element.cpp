#include "matrix_element.hpp"
#include <cassert>
using namespace std;

namespace CasADi{

MatrixElement::MatrixElement(const MX& x, int i_, int j_) : MXNode(x), i(i_), j(j_){
  assert(i>=0 && i<x.size1());
  assert(j>=0 && j<x.size2());
  sz.nrow = 1;
  sz.ncol = 1;
}

MatrixElement* MatrixElement::clone() const{
  return new MatrixElement(*this);
}

void MatrixElement::print(std::ostream &stream) const{
  if(dep(0).size2()==1)
    stream << dep(0) << "[" << i << "]";
  else
    stream << dep(0) << "(" << i << "," << j << ")";
}

void MatrixElement::evaluate(int fsens_order, int asens_order){
  assert(fsens_order==0 || asens_order==0);
  
  if(fsens_order==0){
  // Get references to the terms
  const vector<double>& arg = dep(0)->val(0); // first term
  vector<double>& res = val(0);
  
  // carry out the assignment
  res[0] = arg[j+i*sz.ncol];
  } else {
    // Get references to the terms
    const vector<double>& arg = dep(0)->val(1); // first term
    vector<double>& res = val(1);
  
    // carry out the assignment
    res[0] = arg[j+i*sz.ncol];
  }
  
  if(asens_order>0){
    // Get references to the terms
    vector<double>& arg = dep(0)->val(1); // first term
    const vector<double>& res = val(1);
  
    // carry out the addition
    arg[j+i*sz.ncol] = res[0];
    
  }
}

} // namespace CasADi
