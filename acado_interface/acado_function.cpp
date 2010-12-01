#include "acado_function.hpp"
#include <casadi/stl_vector_tools.hpp>
#include <casadi/expression_tools.hpp>
#include <cassert>
#include <limits>

using namespace std;
namespace CasADi{

AcadoFunction::AcadoFunction(const FX& f) : f_(f){
  // Set pointer to null
  fcn_ = 0;
}
void AcadoFunction::init(){
  if(f_.isNull())
    throw CasadiException("AcadoFunction::init: function is null");

  // Initialize the casadi function
  f_.init();

  // Get dimensions
  dim_.resize(f_.getNumInputs());
  for(int i=0; i<dim_.size(); ++i)
    dim_[i] = f_.input(i).numel();

  // Get number of equations
  neq_ = f_.output().numel();
  
  // Create an acado function
  fcn_ = new ACADO::CFunction( neq_, fcn_wrapper, fcn_fwd_wrapper, fcn_adj_wrapper);
  fcn_->setUserData(this);
}
     
AcadoFunction::~AcadoFunction(){
  if(fcn_) delete fcn_;
}

void AcadoFunction::fcn_wrapper( double *x, double *res, void *user_data ){
  AcadoFunction *this_ = (AcadoFunction*)user_data;
  this_->fcn(x,res);
}

void AcadoFunction::fcn_fwd_wrapper(int number, double *x, double *seed, double *f, double *df, void *user_data){
  AcadoFunction *this_ = (AcadoFunction*)user_data;
  this_->fcn_fwd(number,x,seed,f,df);
}

void AcadoFunction::fcn_adj_wrapper(int number, double *x, double *seed, double *f, double *df, void *user_data){
  AcadoFunction *this_ = (AcadoFunction*)user_data;
  this_->fcn_adj(number,x,seed,f,df);
}

void AcadoFunction::fcn( double *x, double *res){
  int ind = 0;
  for(int i=0; i<dim_.size(); ++i)
    if(dim_[i]>0){
      f_.setInput(&x[ind],i); 
      ind += dim_[i];
    }
  
  // Evaluate
  f_.evaluate();
  
  // Get result
  f_.getOutput(res);
}
 
void AcadoFunction::fcn_fwd(int number, double *x, double *seed, double *f, double *df){
    // Pass inputs and forward seeds
  int ind = 0;
  for(int i=0; i<dim_.size(); ++i)
    if(dim_[i]>0){
      f_.setInput(&x[ind],i); 
      f_.setFwdSeed(&seed[ind],i);
      ind += dim_[i];
    }
  
  // Evaluate with forward sensitivities
  f_.evaluate(1,0);
  
  // Get result
  f_.getOutput(f);

  // Get sensitivities
  f_.getFwdSens(df);
}

void AcadoFunction::fcn_adj(int number, double *x, double *seed, double *f, double *df){
  // Pass inputs
  int ind = 0;
  for(int i=0; i<dim_.size(); ++i)
    if(dim_[i]>0){
      f_.setInput(&x[ind],i); 
      ind += dim_[i];
    }

  // Pass adjoint seed
  f_.setAdjSeed(seed);

  // Evaluate with adjoint sensitivities
  f_.evaluate(0,1);
  
  // Get result
  f_.getOutput(f);

  // Get adjoint sensitivities
  ind = 0;
  for(int i=0; i<dim_.size(); ++i)
    if(dim_[i]>0){
      f_.getAdjSens(&df[ind],i); 
      ind += dim_[i];
    }
}




} // namespace CasADi

