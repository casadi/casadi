#include "vertcat.hpp"
#include "../stl_vector_tools.hpp"
#include <cassert>
#include <iterator>

using namespace std;

namespace CasADi{

// Constructor
Vertcat::Vertcat(const vector<MX>& dep__) : MXNode(dep__){
  assert(!dep_.empty());
  int sz1=0;
  int sz2=dep(0).size2();
  for(vector<MX>::const_iterator it=dep_.begin(); it!=dep_.end(); ++it){
    sz1 += it->size1();
    assert(sz2==it->size2());
  }  
  sz.nrow = sz1;
  sz.ncol = sz2;
}

Vertcat* Vertcat::clone() const{
  return new Vertcat(dep_);
}

void Vertcat::print(ostream &stream) const{
  stream << "[";
  copy(dep_.begin(), dep_.end(), ostream_iterator<MX>(stream, ";"));
  stream << "]";
}

void Vertcat::evaluate(int fsens_order, int asens_order){
  assert(fsens_order==0 || asens_order==0);
  
  if(fsens_order==0){
    int i = 0;
    for(vector<MX>::const_iterator it=dep_.begin(); it!=dep_.end(); ++it){
      copy((*it)->val(0).begin(),(*it)->val(0).end(),&val(0)[i]);
      i += it->numel();
    }
  } else {
    int i = 0;
    for(vector<MX>::const_iterator it=dep_.begin(); it!=dep_.end(); ++it){
      copy((*it)->val(1).begin(),(*it)->val(1).end(),&val(1)[i]);
      i += it->numel();
    }
  }
  
  if(asens_order>0){
    int i = 0;
    for(vector<MX>::iterator it=dep_.begin(); it!=dep_.end(); ++it){
      copy(&val(1)[i],&val(1)[i] + it->numel(), (*it)->val(1).begin());
      i += it->numel();
    }
  }
}

// void Vertcat::setOutput(const vector<double>& x, int ord){
//   int i = 0;
//   for(vector<MX>::iterator it=dep_.begin(); it!=dep_.end(); ++it){
//     (*it)->setOutput(&x[i],ord);
//     i += it->size();
//   }
// }
// 
// void Vertcat::getOutput(vector<double>& x, int ord) const{
//   int i = 0;
//   for(vector<MX>::const_iterator it=dep_.begin(); it!=dep_.end(); ++it){
//     (*it)->getOutput(&x[i],ord);
//     i += it->size();
//   }
// }



} // namespace CasADi
