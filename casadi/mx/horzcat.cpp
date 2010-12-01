#include "horzcat.hpp"
#include "../stl_vector_tools.hpp"
#include <cassert>
#include <iterator>

using namespace std;

namespace CasADi{

// Constructor
Horzcat::Horzcat(const vector<MX>& dep__) : MXNode(dep__){
  assert(!dep_.empty());
  int sz1=dep_[0].size1();
  int sz2=0;
  for(vector<MX>::const_iterator it=dep_.begin(); it!=dep_.end(); ++it){
    assert(sz1==it->size1());
    sz2 += it->size2();
  }  
  sz.nrow = sz1;
  sz.ncol = sz2;
}

Horzcat* Horzcat::clone() const{
  return new Horzcat(dep_);
}

void Horzcat::print(ostream &stream) const{
  stream << "[";
  copy(dep_.begin(), dep_.end(), ostream_iterator<MX>(stream, ","));
  stream << "]";
}

void Horzcat::evaluate(int fsens_order, int asens_order){
  assert(0);
}

} // namespace CasADi
