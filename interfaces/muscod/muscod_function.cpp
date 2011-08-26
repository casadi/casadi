#include "muscod_function.hpp"
#include <casadi/stl_vector_tools.hpp>
#include <casadi/expression_tools.hpp>
#include <cassert>
#include <limits>

using namespace std;

namespace CasADi{

// Get the previous pointers recursively
template<int num_ptrs>
void getFunctionPtrs(muscodFunctionPtr* v){
  getFunctionPtrs<num_ptrs-1>(v);
  
  // Add the last pointer
  v[num_ptrs-1] = fcn_template<num_ptrs-1>;
}

// Stop the recursion when num_ptrs reaches zero
template<>
void getFunctionPtrs<0>(muscodFunctionPtr* v){
}


vector<muscodFunctionPtr> MuscodFunction::generate_functions(){
  // Generate list of pointers
  vector<muscodFunctionPtr> v(MAX_NUM_INSTANCES);
  getFunctionPtrs<MAX_NUM_INSTANCES>(vecptr(v));
  return v;
}

vector<muscodFunctionPtr> MuscodFunction::functions_ = generate_functions();
vector<MuscodFunction*> MuscodFunction::instances_;
stack<int> MuscodFunction::free_instances_;
      
void MuscodFunction::fcn(double  *t, double *xd, double *xa, double *u, double *p, double *rhs,double *rwh, long *iwh,  long *info){
  // Pass inputs
  f_.setInput(t,MUSCOD_FCN_T);
  f_.setInput(xd,MUSCOD_FCN_XD);
  f_.setInput(xa,MUSCOD_FCN_XA);
  f_.setInput(u,MUSCOD_FCN_U);
  f_.setInput(p,MUSCOD_FCN_P);

  // Evaluate
  f_.evaluate();
  
  // Get output
  f_.getOutput(rhs,MUSCOD_FCN_RHS);
//  f_.getOutput(rwh,MUSCOD_FCN_RES);
}

MuscodFunction::MuscodFunction(const FX& f) : f_(f){
  // Check if there are indices that can be reused
  if(free_instances_.empty()){
    instance_no_ = instances_.size();
    if(instance_no_>=MAX_NUM_INSTANCES)
      throw CasadiException("MuscodFunction::MuscodFunction: too many instances");
  
    // Add a pointer to the object
    instances_.push_back(this);
  } else {
    // Use a free index
    instance_no_ = free_instances_.top();
    free_instances_.pop();
    instances_[instance_no_] = this;
  }
}

MuscodFunction::MuscodFunction(const MuscodFunction& fcn) : f_(fcn.f_){
  // Check if there are indices that can be reused
  if(free_instances_.empty()){
    instance_no_ = instances_.size();
    if(instance_no_>=MAX_NUM_INSTANCES)
      throw CasadiException("MuscodFunction::MuscodFunction: too many instances");
  
    // Add a pointer to the object
    instances_.push_back(this);
  } else {
    // Use a free index
    instance_no_ = free_instances_.top();
    free_instances_.pop();
    instances_[instance_no_] = this;
  }
}

MuscodFunction& MuscodFunction::operator=(const MuscodFunction& fcn){
  f_ = fcn.f_;
}


MuscodFunction::~MuscodFunction(){
  free_instances_.push(instance_no_);
  instances_[instance_no_] = 0;
}

muscodFunctionPtr MuscodFunction::getPtr(){
  if(f_.isNull())
    return 0;
  else
    return functions_[instance_no_];
}


} // namespace CasADi

