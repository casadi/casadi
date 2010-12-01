#include "evaluation.hpp"
#include "../stl_vector_tools.hpp"
#include <cassert>

using namespace std;

namespace CasADi{

// Constructor
Evaluation::Evaluation(const FX& fcn, const vector<MX>& dep, int oind_) : MXNode(dep), fcn_(fcn), oind(oind_) {
  sz.nrow = fcn_->output_[oind].size1();
  sz.ncol = fcn_->output_[oind].size2();
}

void Evaluation::print(ostream &stream) const{
  stream << fcn_ << "[" << dep() << "]";
}

void Evaluation::evaluate(int fsens_order, int asens_order){
  assert(fsens_order==0 || asens_order==0);
  
  // Pass the input to the function
  for(int i=0; i<ndep(); ++i)
    fcn_.input(i).set(dep(i)->val(0));

  // Give the forward seed to the function
  if(fsens_order>0) 
    for(int i=0; i<ndep(); ++i)
      fcn_.input(i).setF(dep(i)->val(1));
  
  // Evaluate
  fcn_.evaluate(fsens_order, asens_order);
  
  // Get the results
  fcn_.output(oind).get(&val(0)[0]);
  if(fsens_order>0)
    fcn_.output(oind).getF(val(1));

  // Adjoints
  if(asens_order>0){
    // Pass the adjoint seed to the function
    fcn_.output(oind).setA(val(1));
    
    // Evaluate
    fcn_->evaluate(0,1);

    // Get the results
    for(int i=0; i<ndep(); ++i)
      fcn_.input(i).getA(dep(i)->val(1));
  }
}

} // namespace CasADi
