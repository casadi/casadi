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

INPUTSCHEME(ControlSimulatorInput)

using namespace std;
namespace CasADi{

  
ControlSimulatorInternal::ControlSimulatorInternal(const FX& dae, const FX& output_fcn, const vector<double>& gridc) : dae_(dae), output_fcn_(output_fcn), gridc_(gridc){
  setOption("name","unnamed piecewise simulator");
  addOption("np",OT_INTEGER,GenericType(),"The number of parameters. If this option is not set, all of input(INTEGRATOR_P) is considered static parameters. The remainder of nv = input(INTEGRATOR_P) is considered to be varying parameters.");
  addOption("nf",OT_INTEGER,1,"Number of fine grained integration steps.");
  addOption("integrator",               OT_INTEGRATOR, GenericType(), "An integrator creator function");
  addOption("integrator_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the integrator");
}
  
ControlSimulatorInternal::~ControlSimulatorInternal(){
}

void ControlSimulatorInternal::init(){

  // Initialize dae if necessary
  if (!dae_.isInit()) dae_.init();
 
  // Create an integrator instance
  integratorCreator integrator_creator = getOption("integrator");
  integrator_ = integrator_creator(parameterizeTime(dae_),FX());
  if(hasSetOption("integrator_options")){
    integrator_.setOption(getOption("integrator_options"));
  }
  
  
  // Size of the coarse grid
  ns_ = gridc_.size();
    
  // Number of fine-grained steps
  nf_ = getOption("nf");
  
  // Populate the fine-grained grid_
  if (nf_==1) { 
 	  // The default case: don't change the grid 
 	  grid_ = gridc_; 
 	} else { 
 	  // Interpolate the grid. 
 	  grid_.resize((gridc_.size()-1)*nf_+1); 	     
 	  std::vector< double > refined(nf_+1,0); 	 
 	  for (int k=0;k<gridc_.size()-1;++k) { 
      linspace(refined,gridc_[k],gridc_[k+1]); 
      std::copy(refined.begin(),refined.end()-1,grid_.begin()+k*nf_); 
 	  } 
         
 	  grid_[grid_.size()-1] = gridc_[gridc_.size()-1]; 
 	} 
  
  casadi_assert_message(nf_>=0,"Invalid parameter. nf must be at least 1.");
  
  // Retrieve np option, do some checking and assign np_, nv_
  if (hasSetOption("np")) {
    np_ = getOption("np");
    casadi_assert_message(np_<=dae_.input(DAE_P).size(),"Invalid parameter. np (" << np_ << ") cannot be greater that dae.input(DAE_P), which is of size " << dae_.input(DAE_P).size() << ".");
    casadi_assert_message(np_>=0,"Invalid parameter. np (" << np_ << ") must be greater than or equal to zero.");
    nv_ = dae_.input(DAE_P).size() - np_;
  } else {
    np_ = dae_.input(DAE_P).size();
    nv_ = 0;
  }
  
  // Let the integration time start from the np_first point of the time grid.
  if (!gridc_.empty()) integrator_.setOption("t0",gridc_[0]);

  // Initialize the integrator
  integrator_.init();
  
  // Generate an output function if there is none (returns the whole state)
  if(output_fcn_.isNull()){
    SXMatrix t = ssym("t");
    SXMatrix x = ssym("x",dae_.input(DAE_Y).sparsity());
    SXMatrix xdot = ssym("xp",dae_.input(DAE_YDOT).sparsity());
    SXMatrix p = ssym("p",dae_.input(DAE_P).sparsity());

    vector<SXMatrix> arg(DAE_NUM_IN);
    arg[DAE_T] = t;
    arg[DAE_Y] = x;
    arg[DAE_P] = p;
    arg[DAE_YDOT] = xdot;

    vector<SXMatrix> out(INTEGRATOR_NUM_OUT);
    out[INTEGRATOR_XF] = x;
    out[INTEGRATOR_XPF] = xdot;

    // Create the output function
    output_fcn_ = SXFunction(arg,out);
  }

  // Initialize the output function
  output_fcn_.init();
  
    
  // Extend the output function two extra outputs at the start: DAE_Y and DAE_YDOT
  vector<MX> output_fcn_in_ = output_fcn_.symbolicInput();
  
  vector<MX> output_fcn_out_(2 + output_fcn_.getNumOutputs());
  output_fcn_out_[0] = output_fcn_in_[DAE_Y];
  output_fcn_out_[1] = output_fcn_in_[DAE_YDOT];
  
  vector<MX> output_fcn_call_ = output_fcn_.call(output_fcn_in_);
  copy(output_fcn_call_.begin(),output_fcn_call_.end(),output_fcn_out_.begin()+2);
  
  output_fcn_ = MXFunction(output_fcn_in_,output_fcn_out_);

  // Initialize the output function again
  output_fcn_.init();
  
  // Make a local grid with non-dimensional time
  gridlocal_.resize(nf_+1);
  linspace(gridlocal_,0,1); 
  
  // Create the simulator
  simulator_ = Simulator(integrator_,parameterizeTimeOutput(output_fcn_),gridlocal_);
  simulator_.init();
  
  // Allocate inputs
  input_.resize(CONTROLSIMULATOR_NUM_IN);
  input(CONTROLSIMULATOR_X0)  = dae_.input(DAE_Y);
  input(CONTROLSIMULATOR_P)   = dae_.input(DAE_P)[range(np_)];
  input(CONTROLSIMULATOR_V)   = trans(repmat(integrator_.input(INTEGRATOR_P)[range(np_,np_+nv_)],1,ns_-1));
  input(CONTROLSIMULATOR_XP0) = dae_.input(DAE_YDOT);

  // Allocate outputs
  output_.resize(output_fcn_->output_.size()-2);
  for(int i=0; i<output_.size(); ++i)
    output(i) = Matrix<double>((ns_-1)*nf_+1,output_fcn_.output(i+2).numel(),0);

  // Call base class method
  FXInternal::init();
  
  
  // Variables on which the chain of simulator calls (all_output_) depend
  MX Xk("Xk", input(CONTROLSIMULATOR_X0).size());
  MX XPk("XPk", input(CONTROLSIMULATOR_XP0).size());
  MX P("P",input(CONTROLSIMULATOR_P).size());
  MX V("V",input(CONTROLSIMULATOR_V).sparsity());
 
  // Group these variables as an input list for all_output_
  vector<MX> all_output_in(CONTROLSIMULATOR_NUM_IN);
  all_output_in[CONTROLSIMULATOR_X0] = Xk;
  all_output_in[CONTROLSIMULATOR_XP0] = XPk;
  all_output_in[CONTROLSIMULATOR_P] = P;
  all_output_in[CONTROLSIMULATOR_V] = V;
  
  // Placeholder with which simulator.input(INTEGRATOR_P) will be fed [t0 tf pfixed pvariable]
  vector<MX> P_eval(4);
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
    simulator_in[INTEGRATOR_XP0] = XPk;
    P_eval[0] = MX(gridc_[k]);
    P_eval[1] = MX(gridc_[k+1]);
    P_eval[3] = trans(V(k,range(nv_)));
    simulator_in[INTEGRATOR_P] = vertcat(P_eval);
    
    simulator_out = simulator_.call(simulator_in);
    
    // Remember the end state and dstate for next iteration in this loop
    Xk = trans(simulator_out[0](simulator_out[0].size1()-1,range(simulator_out[0].size2())));
    XPk = trans(simulator_out[1](simulator_out[1].size1()-1,range(simulator_out[1].size2())));  
    
    // Copy all the outputs (but not those 2 extra we introduced)
    for (int i=0;i<simulator_out.size()-2;++i) {
      if (k==0) {  // We only need the initial state for the very first simulator
        simulator_outputs[i].push_back(simulator_out[i+2]);
      } else {  // For the rest, we can chop the initial state of
        simulator_outputs[i].push_back(simulator_out[i+2](range(1,nf_+1),range(simulator_out[i+2].size2())));
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
 	  Matrix<double> ret(grid_.size()-1,nv_,0);
 	  for (int i=0;i<ns_-1;++i) {
 	    for (int k=0;k<nf_;++k) {
 	      copy(input(CONTROLSIMULATOR_V).data().begin()+i*nv_,input(CONTROLSIMULATOR_V).data().begin()+(i+1)*nv_,ret.begin()+i*nv_*nf_+k);
 	    }
 	  }
 	  return ret;
}

std::vector< int > ControlSimulatorInternal::getCoarseIndex() const {
   return range(0,grid_.size(),nf_);
}
 	

} // namespace CasADi


