#include "transpose.hpp"
#include <cassert>

using namespace std;

namespace CasADi{

Transpose::Transpose(const MX& x) : MXNode(x){
  sz.nrow = x.size2();
  sz.ncol = x.size1();
}

Transpose* Transpose::clone() const{
  return new Transpose(*this);
}

void Transpose::print(std::ostream &stream) const{
  stream << "trans(" << dep(0) << ")";
}

void Transpose::evaluate(int fsens_order, int asens_order){
  if(fsens_order==0 || asens_order==0);
  
  if(fsens_order==0){
  // Get references to the terms
  const vector<double>& arg = dep(0)->val(0);
  vector<double>& res = val(0);
  
  // carry out the transpose
  for(int i=0; i<sz.nrow; ++i)
    for(int j=0; j<sz.ncol; ++j)
      res[j + i*sz.ncol] = arg[i + j*sz.nrow];
  } else {

    // Get references to the terms
    const vector<double>& arg = dep(0)->val(1);
    vector<double>& res = val(1);
  
    // carry out the transpose
    for(int i=0; i<sz.nrow; ++i)
      for(int j=0; j<sz.ncol; ++j)
        res[j + i*sz.ncol] = arg[i + j*sz.nrow];
  }
  
  if(asens_order>0){
    // Get references to the terms
    vector<double>& arg = dep(0)->val(1);
    const vector<double>& res = val(1);
  
    // carry out the transpose
    for(int i=0; i<sz.nrow; ++i)
      for(int j=0; j<sz.ncol; ++j)
        arg[i + j*sz.nrow] += res[j + i*sz.ncol];
    
  }
}

} // namespace CasADi
