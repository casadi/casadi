/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "control_simulator_internal.hpp"
#include "integrator_internal.hpp"
#include "../stl_vector_tools.hpp"
#include "sx_function.hpp"
#include "../sx/sx_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "fx_tools.hpp"
#include <utility>
#include <string>

INPUTSCHEME(ControlSimulatorInput)

using namespace std;
namespace CasADi{

  
ControlSimulatorInternal::ControlSimulatorInternal(const FX& control_dae, const FX& output_fcn, const vector<double>& gridc) : control_dae_(control_dae), orig_output_fcn_(output_fcn), gridc_(gridc){
  setOption("name","unnamed controlsimulator");
  addOption("nf",OT_INTEGER,1,"Number of minor grained integration steps per major interval. nf>0 must hold. This option is not used when 'minor_grid' is provided.");
  addOption("minor_grid",OT_INTEGERVECTOR,GenericType(),"The local grid used on each major interval, with time normalized to 1. By default, option 'nf' is used to construct a linearly spaced grid.");
  addOption("integrator",               OT_INTEGRATOR, GenericType(), "An integrator creator function");
  addOption("integrator_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the integrator");
  addOption("simulator_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the simulator");
  addOption("control_interpolation",   OT_STRING,     "none", "none|nearest|linear");
  addOption("control_endpoint",        OT_BOOLEAN,       false, "Include a control value at the end of the simulation domain. Used for interpolation.");
  
  inputScheme_ = SCHEME_ControlSimulatorInput;
}
  
ControlSimulatorInternal::~ControlSimulatorInternal(){
}


void ControlSimulatorInternal::init(){

  bool control_endpoint = getOption("control_endpoint");
  
  if (!control_dae_.isInit()) control_dae_.init();
  
  casadi_assert_message(!gridc_.empty(),"The supplied time grid must not be empty.");
    
  casadi_assert_message(isIncreasing(gridc_),"The supplied time grid must be strictly increasing. Notably, you cannot have a time instance repeating."); 

  if (control_dae_.getNumInputs()==DAE_NUM_IN) {
    vector<MX> control_dae_in_(CONTROL_DAE_NUM_IN);
    vector<MX> dae_in_ = control_dae_.symbolicInput();
    control_dae_in_[CONTROL_DAE_T]    = dae_in_[DAE_T];
    control_dae_in_[CONTROL_DAE_X]    = dae_in_[DAE_X];
    control_dae_in_[CONTROL_DAE_Z]    = dae_in_[DAE_Z];
    control_dae_in_[CONTROL_DAE_P]    = dae_in_[DAE_P];
    control_dae_ = MXFunction(control_dae_in_,control_dae_.call(dae_in_));
    control_dae_.init();
  }
  casadi_assert_message(control_dae_.getNumInputs()==CONTROL_DAE_NUM_IN,"ControlSimulatorInternal::init: supplied control_dae does not conform to the CONTROL_DAE or DAE input scheme.");
  
  // Cast control_dae in a form that integrator can manage
  vector<MX> dae_in_(DAE_NUM_IN);
  
  dae_in_[DAE_T]    = MX("tau",control_dae_.input(CONTROL_DAE_T).sparsity());
  dae_in_[DAE_X]    = MX("x",control_dae_.input(CONTROL_DAE_X).sparsity());
  dae_in_[DAE_Z]    = MX("z",control_dae_.input(CONTROL_DAE_Z).sparsity());
  
  np_ = control_dae_.input(CONTROL_DAE_P).size();
  
  CRSSparsity u_sparsity = control_dae_.input(CONTROL_DAE_U).sparsity();
  if (control_dae_.input(CONTROL_DAE_U_INTERP).size()!=0) {
    u_sparsity = control_dae_.input(CONTROL_DAE_U_INTERP).sparsity();
  }
  
  nu_ = control_dae_.input(CONTROL_DAE_U).size();
  
  if (!control_dae_.input(CONTROL_DAE_U_INTERP).empty() && nu_!=0) {
    casadi_assert_message(control_dae_.input(CONTROL_DAE_U).sparsity()==control_dae_.input(CONTROL_DAE_U_INTERP).sparsity(),"You specfified both U and U_INTERP, but the sparsities do not match: " << control_dae_.input(CONTROL_DAE_U).dimString() << "  <-> " << control_dae_.input(CONTROL_DAE_U_INTERP).sparsity());
  }
  
  int nu_end   = control_dae_.input(CONTROL_DAE_U_INTERP).size();
  nu_ = std::max(nu_end,nu_);

  ny_ = control_dae_.input(CONTROL_DAE_X).size();
  
  // Structure of DAE_P : T0 TF P Ustart Uend Y_MAJOR 
  
  dae_in_[DAE_P]    = MX("P",2+np_+nu_+nu_end+ny_,1);

  int iT0 = 0;
  int iTF = 1;
  IMatrix iP  = 2+IMatrix(control_dae_.input(CONTROL_DAE_P).sparsity(),range(np_));
  IMatrix iUstart  = 2+np_ + IMatrix(u_sparsity,range(nu_));
  IMatrix iUend    = 2+np_ + nu_ + IMatrix(control_dae_.input(CONTROL_DAE_U_INTERP).sparsity(),range(nu_end));
  IMatrix iYM;
  if (!control_dae_.input(CONTROL_DAE_X_MAJOR).empty()) {
    iYM = 2+np_ + nu_ + nu_end + IMatrix(control_dae_.input(CONTROL_DAE_X_MAJOR).sparsity(),range(ny_));
  }
  
  vector<MX> control_dae_in_(CONTROL_DAE_NUM_IN);
  
  if (!control_dae_.input(CONTROL_DAE_T).empty())
    control_dae_in_[CONTROL_DAE_T]        = dae_in_[DAE_P](iT0) + (dae_in_[DAE_P](iTF)-dae_in_[DAE_P](iT0))*dae_in_[DAE_T];
  if (!control_dae_.input(CONTROL_DAE_T0).empty())
    control_dae_in_[CONTROL_DAE_T0]       = dae_in_[DAE_P](iT0);
  if (!control_dae_.input(CONTROL_DAE_TF).empty())
    control_dae_in_[CONTROL_DAE_TF]       = dae_in_[DAE_P](iTF);
  control_dae_in_[CONTROL_DAE_X]        = dae_in_[DAE_X];
  if (!control_dae_.input(CONTROL_DAE_Z).empty())
    control_dae_in_[CONTROL_DAE_Z]        = dae_in_[DAE_Z];
  control_dae_in_[CONTROL_DAE_P]        = dae_in_[DAE_P](iP);
  if (!control_dae_.input(CONTROL_DAE_U).empty())
    control_dae_in_[CONTROL_DAE_U]        = dae_in_[DAE_P](iUstart);
  if (!control_dae_.input(CONTROL_DAE_U_INTERP).empty()) {
    MX tau = (dae_in_[DAE_P](iTF)-dae_in_[DAE_T])/(dae_in_[DAE_P](iTF)-dae_in_[DAE_P](iT0));
    control_dae_in_[CONTROL_DAE_U_INTERP] = dae_in_[DAE_P](iUstart) * (1-tau) + tau* dae_in_[DAE_P](iUend);
  }
  if (!control_dae_.input(CONTROL_DAE_X_MAJOR).empty())
    control_dae_in_[CONTROL_DAE_X_MAJOR] = dae_in_[DAE_P](iYM);

  std::vector<MX> control_dae_call = control_dae_.call(control_dae_in_);
  
  
  std::vector<MX> dae_out(daeOut("ode",(dae_in_[DAE_P](iTF)-dae_in_[DAE_P](iT0))*control_dae_call[DAE_ODE]));

  int i=1;
  while( control_dae_call.size()>i && dae_out.size()>i) {dae_out[i] = control_dae_call[i];i++;}
   
  dae_ = MXFunction(dae_in_,dae_out);
  
  dae_.init();
 
  // Create an integrator instance
  integratorCreator integrator_creator = getOption("integrator");
  integrator_ = integrator_creator(dae_,FX());
  if(hasSetOption("integrator_options")){
    integrator_.setOption(getOption("integrator_options"));
  }
  
  // Size of the coarse grid
  ns_ = gridc_.size();
    
  // Number of fine-grained steps
  nf_ = getOption("nf");
  
  casadi_assert_message(nf_>0,"Option 'nf' must be greater than zero.");
  
  if (hasSetOption("minor_grid")) {
    gridlocal_ = getOption("minor_grid");
    nf_ = gridlocal_.size()-1;
    casadi_assert_message(gridlocal_.size()>1,"Option 'minor_grid' must have more then one element.");
  } else {
    gridlocal_.resize(nf_+1);
    linspace(gridlocal_,0,1);
  }

  
  // Populate the fine-grained grid_
  if (nf_==1) { 
           // The default case: don't change the grid 
           grid_ = gridc_; 
         } else { 
           // Interpolate the grid. 
           grid_.resize((gridc_.size()-1)*nf_+1);                       
           for (int k=0;k<gridc_.size()-1;++k) { 
      for (int i=0;i<gridlocal_.size()-1;++i) {
        grid_[k*nf_+i] = gridc_[k] + gridlocal_[i]*(gridc_[k+1]-gridc_[k]);
      }
           } 
         
           grid_[grid_.size()-1] = gridc_[gridc_.size()-1]; 
         } 
  
  // Let the integration time start from the np_first point of the time grid.
  if (!gridc_.empty()) integrator_.setOption("t0",gridc_[0]);

  // Initialize the integrator
  integrator_.init();
  
  // Generate an output function if there is none (returns the whole state)
  if(orig_output_fcn_.isNull()){
    
    SXMatrix t        = ssym("t",control_dae_.input(CONTROL_DAE_T).sparsity());
    SXMatrix t0       = ssym("t0",control_dae_.input(CONTROL_DAE_T0).sparsity());
    SXMatrix tf       = ssym("tf",control_dae_.input(CONTROL_DAE_TF).sparsity());
    SXMatrix x        = ssym("x",control_dae_.input(CONTROL_DAE_X).sparsity());
    SXMatrix z        = ssym("z",control_dae_.input(CONTROL_DAE_Z).sparsity());
    SXMatrix p        = ssym("p",control_dae_.input(CONTROL_DAE_P).sparsity());
    SXMatrix u        = ssym("u",control_dae_.input(CONTROL_DAE_U).sparsity());
    SXMatrix u_interp = ssym("u_interp",control_dae_.input(CONTROL_DAE_U_INTERP).sparsity());
    SXMatrix x0       = ssym("x0",control_dae_.input(CONTROL_DAE_X_MAJOR).sparsity());
  
    vector<SXMatrix> arg(CONTROL_DAE_NUM_IN);
    arg[CONTROL_DAE_T]  = t;
    arg[CONTROL_DAE_X]  = x;
    arg[CONTROL_DAE_P]  = p;
    arg[CONTROL_DAE_Z] = z;
    arg[CONTROL_DAE_X_MAJOR] = x0;
    arg[CONTROL_DAE_U] = u;
    arg[CONTROL_DAE_U_INTERP] = u_interp;
    arg[CONTROL_DAE_T0] = t0;
    arg[CONTROL_DAE_TF] = tf;
    
    vector<SXMatrix> out(INTEGRATOR_NUM_OUT);
    out[INTEGRATOR_XF] = x;

    // Create the output function
    output_fcn_ = SXFunction(arg,out);
    output_fcn_.setOption("name","output");
    outputScheme_ = SCHEME_IntegratorOutput;
  } else {
    output_fcn_ = orig_output_fcn_;
  }

 
 
  // Initialize the output function
  output_fcn_.init();
  
  casadi_assert_message(output_fcn_.getNumInputs()==CONTROL_DAE_NUM_IN,"Output function, if supplied, must adhere to ControlledDAEInput scheme.");
    
  // Extend the output function two extra outputs at the start: DAE_X and DAE_XDOT
  vector<MX> output_fcn_in_ = output_fcn_.symbolicInput();
  
  vector<MX> output_fcn_out_(2 + output_fcn_.getNumOutputs());
  output_fcn_out_[0] = output_fcn_in_[CONTROL_DAE_X];
  
  vector<MX> output_fcn_call_ = output_fcn_.call(output_fcn_in_);
    
  copy(output_fcn_call_.begin(),output_fcn_call_.end(),output_fcn_out_.begin()+2);
  
  output_fcn_ = MXFunction(output_fcn_in_,output_fcn_out_);

  // Initialize the output function again
  output_fcn_.init();
  
  if (!output_fcn_.input(CONTROL_DAE_U).empty() && control_dae_.input(CONTROL_DAE_U).empty()) {
    casadi_error("ControlSimulatorInternal::init: output function has CONTROL_DAE_U input. The supplied DAE should have it as well.");
  }
  if (!output_fcn_.input(CONTROL_DAE_U_INTERP).empty() && control_dae_.input(CONTROL_DAE_U_INTERP).empty()) {
    casadi_error("ControlSimulatorInternal::init: output function has CONTROL_DAE_U_INTERP input. The supplied DAE should have it as well.");
  }
  if (!output_fcn_.input(CONTROL_DAE_X_MAJOR).empty() && control_dae_.input(CONTROL_DAE_X_MAJOR).empty()) {
    casadi_error("ControlSimulatorInternal::init: output function has CONTROL_DAE_X_MAJOR input. The supplied DAE should have it as well.");
  }
   
  output_fcn_in_ = vector<MX>(CONTROL_DAE_NUM_IN);

  dae_in_[DAE_T] = msym("tau");
  if (!output_fcn_.input(CONTROL_DAE_T).empty()) {
    output_fcn_in_[CONTROL_DAE_T]        = dae_in_[DAE_P](iT0) + (dae_in_[DAE_P](iTF)-dae_in_[DAE_P](iT0))*dae_in_[DAE_T];
  }
  if (!output_fcn_.input(CONTROL_DAE_T0).empty())
    output_fcn_in_[CONTROL_DAE_T0]       = dae_in_[DAE_P](iT0);
  if (!output_fcn_.input(CONTROL_DAE_TF).empty())
    output_fcn_in_[CONTROL_DAE_TF]       = dae_in_[DAE_P](iTF);
  output_fcn_in_[CONTROL_DAE_X]        = dae_in_[DAE_X];
  output_fcn_in_[CONTROL_DAE_P]        = dae_in_[DAE_P](iP);
  if (!output_fcn_.input(CONTROL_DAE_U).empty())
    output_fcn_in_[CONTROL_DAE_U]        = dae_in_[DAE_P](iUstart);
  if (!output_fcn_.input(CONTROL_DAE_U_INTERP).empty()) {
    output_fcn_in_[CONTROL_DAE_U_INTERP] = dae_in_[DAE_P](iUstart) * (1-dae_in_[DAE_T]) + dae_in_[DAE_T]* dae_in_[DAE_P](iUend);
  }
  if (!output_fcn_.input(CONTROL_DAE_X_MAJOR).empty())
    output_fcn_in_[CONTROL_DAE_X_MAJOR] = dae_in_[DAE_P](iYM);

  // Transform the output_fcn_ with CONTROL_DAE input scheme to a DAE input scheme
  output_fcn_ = MXFunction(dae_in_,output_fcn_.call(output_fcn_in_));
  
  // Initialize the output function again
  output_fcn_.init();
  
  // Create the simulator
  simulator_ = Simulator(integrator_,output_fcn_,gridlocal_);
  if(hasSetOption("simulator_options")){
    simulator_.setOption(getOption("simulator_options"));
  }
  simulator_.init();
  
  // Allocate inputs
  input_.resize(CONTROLSIMULATOR_NUM_IN);
  input(CONTROLSIMULATOR_X0)  = DMatrix(dae_.input(DAE_X));
  input(CONTROLSIMULATOR_P)   = control_dae_.input(CONTROL_DAE_P);
  input(CONTROLSIMULATOR_U)   = DMatrix(ns_ - 1 + (control_endpoint ? 1 : 0) ,nu_,0);

  // Allocate outputs
  output_.resize(output_fcn_->output_.size()-2);
  for(int i=0; i<output_.size(); ++i)
    output(i) = Matrix<double>((ns_-1)*nf_+1,output_fcn_.output(i+2).numel(),0);

  // Call base class method
  FXInternal::init();
  
  // Variables on which the chain of simulator calls (all_output_) depend
  MX Xk("Xk", input(CONTROLSIMULATOR_X0).size());
  MX P("P",input(CONTROLSIMULATOR_P).size());
  MX U("U",input(CONTROLSIMULATOR_U).sparsity());
 
  // Group these variables as an input list for all_output_
  vector<MX> all_output_in(CONTROLSIMULATOR_NUM_IN);
  all_output_in[CONTROLSIMULATOR_X0] = Xk;
  all_output_in[CONTROLSIMULATOR_P] = P;
  all_output_in[CONTROLSIMULATOR_U] = U;
  
  // Placeholder with which simulator.input(INTEGRATOR_P) will be fed [t0 tf P Ustart Uend Y_MAJOR]
  vector<MX> P_eval(6);
  P_eval[2] = P; // We can already set the fixed part in advance.
  
  // Placeholder to collect the outputs of all simulators (but not those 2 extra we introduced)
  vector< vector<MX> > simulator_outputs(output_fcn_.getNumOutputs()-2);

  // Input arguments to simulator.call
  vector<MX> simulator_in(INTEGRATOR_NUM_IN);
  // Output of simulator.call
  vector<MX> simulator_out;
  
  for(int k=0; k<ns_-1; ++k){
    // Set the appropriate inputs for the k'th simulator call
    simulator_in[INTEGRATOR_X0] = Xk;
    P_eval[0] = MX(gridc_[k]);
    P_eval[1] = MX(gridc_[k+1]);
    if (nu_>0) {
      P_eval[3] = trans(U(k,range(nu_)));
    }
    if (control_dae_.input(CONTROL_DAE_U_INTERP).size()>0) {
      if (k+1==U.size1()) {
        P_eval[4] = P_eval[3];
      } else {
        P_eval[4] = trans(U(k+1,range(nu_)));
      }
    }
    P_eval[5] = Xk;
    
    simulator_in[INTEGRATOR_P] = vertcat(P_eval);
 
    simulator_out = simulator_.call(simulator_in);
    
    // Remember the end state and dstate for next iteration in this loop
    Xk = trans(simulator_out[0](simulator_out[0].size1()-1,ALL));
    
    // Copy all the outputs (but not those 2 extra we introduced)
    for (int i=0;i<simulator_out.size()-2;++i) {
      if(simulator_out[i+2].isNull()) continue; // NOTE: Joel: quick-fix
      
      simulator_outputs[i].push_back(simulator_out[i+2](range(nf_),ALL));
      if (k+1==ns_-1) {  // Output of the last minor step of the last major step
        simulator_outputs[i].push_back(simulator_out[i+2](std::vector<int>(1,nf_),ALL));
      }
    }
    
  }
  

  // Concatenate the results of all simulator calls
  vector<MX> all_output_out(simulator_outputs.size());
  for (int k=0;k<all_output_out.size();++k) {
      all_output_out[k] = vertcat(simulator_outputs[k]);
  }
  
  // Finally, construct all_output_
  all_output_ = MXFunction(all_output_in,all_output_out);
  all_output_.init();
  
  ControlSimulatorInternal::updateNumSens(false);
  
}

void ControlSimulatorInternal::evaluate(int nfdir, int nadir){

  // Copy all inputs
  for (int i=0;i<input_.size();++i) {
    all_output_.input(i).set(input(i));
    // Pass forward seeds
    for(int dir=0; dir<nfdir; ++dir){
      all_output_.fwdSeed(i,dir).set(fwdSeed(i,dir));
    }
  }
  
  for (int i=0;i<output_.size();++i) {
    // Pass adjoint seeds
    for(int dir=0; dir<nadir; ++dir){
      all_output_.adjSeed(i,dir).set(adjSeed(i,dir));
    }
  }
  
  all_output_.evaluate(nfdir,nadir);
  
  // Copy all outputs
  for (int i=0;i<output_.size();++i) {
    output(i).set(all_output_.output(i));
    // Copy all forward sensitivities
    for(int dir=0; dir<nfdir; ++dir){
      fwdSens(i,dir).set(all_output_.fwdSens(i,dir));
    }
  }
  
  for (int i=0;i<input_.size();++i) {
    // Copy all adjoint sensitivities
    for(int dir=0; dir<nadir; ++dir){
      adjSens(i,dir).set(all_output_.adjSens(i,dir));
    }
  }
  
  
}

Matrix<double> ControlSimulatorInternal::getVFine() const {
           Matrix<double> ret(grid_.size()-1,nu_,0);
           for (int i=0;i<ns_-1;++i) {
             for (int k=0;k<nf_;++k) {
               copy(input(CONTROLSIMULATOR_U).data().begin()+i*nu_,input(CONTROLSIMULATOR_U).data().begin()+(i+1)*nu_,ret.begin()+i*nu_*nf_+k*nu_);
             }
           }
           return ret;
}

std::vector< int > ControlSimulatorInternal::getCoarseIndex() const {
   return range(0,grid_.size(),nf_);
}

void ControlSimulatorInternal::updateNumSens(bool recursive){
  if (recursive) {
    FXInternal::updateNumSens(recursive);
  }
  
  if (!output_fcn_.isNull()) {
    output_fcn_.setOption("number_of_fwd_dir",getOption("number_of_fwd_dir"));
    output_fcn_.updateNumSens();
  }
  
  if (!integrator_.isNull()) {
    integrator_.setOption("number_of_fwd_dir",getOption("number_of_fwd_dir"));
    integrator_.updateNumSens();
  }
  
  if (!simulator_.isNull()) {
    simulator_.setOption("number_of_fwd_dir",getOption("number_of_fwd_dir"));
    simulator_.updateNumSens();
  }
  
  if (!all_output_.isNull()) {
    all_output_.setOption("number_of_fwd_dir",getOption("number_of_fwd_dir"));
    all_output_.updateNumSens();
  }
  
}

         

} // namespace CasADi


