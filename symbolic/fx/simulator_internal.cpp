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

#include "simulator_internal.hpp"
#include "integrator_internal.hpp"
#include "../stl_vector_tools.hpp"
#include "sx_function.hpp"
#include "../sx/sx_tools.hpp"

INPUTSCHEME(IntegratorInput)

using namespace std;
namespace CasADi{

  
SimulatorInternal::SimulatorInternal(const Integrator& integrator, const FX& output_fcn, const vector<double>& grid) : integrator_(integrator), output_fcn_(output_fcn), grid_(grid){
  setOption("name","unnamed simulator");
  addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "", "initial|step", true);
  
  inputScheme_ = SCHEME_IntegratorInput;
}
  
SimulatorInternal::~SimulatorInternal(){
}

void SimulatorInternal::init(){
  // Let the integration time start from the first point of the time grid.
  if (!grid_.empty()) integrator_.setOption("t0",grid_[0]);
  // Let the integration time stop at the last point of the time grid.
  if (!grid_.empty()) integrator_.setOption("tf",grid_[grid_.size()-1]);
  
  casadi_assert_message(isNonDecreasing(grid_),"The supplied time grid must be non-decreasing."); 
  
  // Initialize the integrator
  integrator_.init();
  
  // Generate an output function if there is none (returns the whole state)
  if(output_fcn_.isNull()){
    SXMatrix t = ssym("t");
    SXMatrix x = ssym("x",integrator_.input(INTEGRATOR_X0).sparsity());
    SXMatrix p = ssym("p",integrator_.input(INTEGRATOR_P).sparsity());

    vector<SXMatrix> arg(DAE_NUM_IN);
    arg[DAE_T] = t;
    arg[DAE_X] = x;
    arg[DAE_P] = p;

    vector<SXMatrix> out(INTEGRATOR_NUM_OUT);
    out[INTEGRATOR_XF] = x;

    // Create the output function
    output_fcn_ = SXFunction(arg,out);
    
    outputScheme_ = SCHEME_IntegratorOutput;
  }

  // Initialize the output function
  output_fcn_.init();
  
  SimulatorInternal::updateNumSens(false);
  
  // Allocate inputs
  input_.resize(INTEGRATOR_NUM_IN);
  for(int i=0; i<INTEGRATOR_NUM_IN; ++i){
    input(i) = integrator_.input(i);
  }

  // Allocate outputs
  output_.resize(output_fcn_->output_.size());
  for(int i=0; i<output_.size(); ++i) {
    output(i) = Matrix<double>(grid_.size(),output_fcn_.output(i).numel(),0);
    if (!output_fcn_.output(i).empty()) {
      casadi_assert_message(output_fcn_.output(i).size2()==1,"SimulatorInternal::init: Output function output #" << i << " has shape " << output_fcn_.output(i).dimString() << ", while a column-matrix shape is expected.");
    }
  }
  
  casadi_assert_message( output_fcn_.input(DAE_T).numel() <=1, "SimulatorInternal::init: output_fcn DAE_T argument must be scalar or empty, but got " << output_fcn_.input(DAE_T).dimString());

  casadi_assert_message( output_fcn_.input(DAE_P).empty() || integrator_.input(INTEGRATOR_P).sparsity()== output_fcn_.input(DAE_P).sparsity(), "SimulatorInternal::init: output_fcn DAE_P argument must be empty or have dimension " << integrator_.input(INTEGRATOR_P).dimString() << ", but got " << output_fcn_.input(DAE_P).dimString());

  casadi_assert_message( output_fcn_.input(DAE_X).empty() || integrator_.input(INTEGRATOR_X0).sparsity()== output_fcn_.input(DAE_X).sparsity(), "SimulatorInternal::init: output_fcn DAE_X argument must be empty or have dimension " << integrator_.input(INTEGRATOR_X0).dimString() << ", but got " << output_fcn_.input(DAE_X).dimString());
  
  // Call base class method
  FXInternal::init();
  
  states_.resize(grid_.size());
  for (int k = 0; k < grid_.size(); ++k) {
    states_[k]=Matrix<double>::zeros(integrator_.input(INTEGRATOR_X0).size1());
  }
    
}

void SimulatorInternal::evaluate(int nfdir, int nadir){
  casadi_assert_message(nadir==0, "Not implemented");
  
  // Pass the parameters and initial state
  integrator_.setInput(input(INTEGRATOR_X0),INTEGRATOR_X0);
  integrator_.setInput(input(INTEGRATOR_P),INTEGRATOR_P);
  
  if (monitored("initial")) {
    std::cout << "SimulatorInternal::evaluate: initial condition:" << std::endl;
    std::cout << " y0     = "  << input(INTEGRATOR_X0) << std::endl;
    std::cout << " p      = "   << input(INTEGRATOR_P) << std::endl;
  }
  
  // Pass sensitivities if fsens 
  for(int dir=0; dir<nfdir; ++dir){ 
    integrator_.setFwdSeed(fwdSeed(INTEGRATOR_X0,dir),INTEGRATOR_X0,dir); 
    integrator_.setFwdSeed(fwdSeed(INTEGRATOR_P,dir),INTEGRATOR_P,dir); 
  }
    
  // Reset the integrator_
  integrator_.reset(nfdir);
  
  // Advance solution in time
  for(int k=0; k<grid_.size(); ++k){

    if (monitored("step")) {
      std::cout << "SimulatorInternal::evaluate: integrating up to: " <<  grid_[k] << std::endl;
      std::cout << " y0       = "  << integrator_.input(INTEGRATOR_X0) << std::endl;
      std::cout << " p        = "   << integrator_.input(INTEGRATOR_P) << std::endl;
    }
  
    // Integrate to the output time
    integrator_.integrate(grid_[k]);

    if (monitored("step")) {
      std::cout << " y_final  = "  << integrator_.output(INTEGRATOR_XF) << std::endl;
    }
    
    // Pass integrator output to the output function
    if(output_fcn_.input(DAE_T).size()!=0)
      output_fcn_.setInput(grid_[k],DAE_T);
    if(output_fcn_.input(DAE_X).size()!=0)
      output_fcn_.setInput(integrator_.output(INTEGRATOR_XF),DAE_X);
    if(output_fcn_.input(DAE_P).size()!=0)
      output_fcn_.setInput(input(INTEGRATOR_P),DAE_P);
      
    // Save the states for use in backwards sensitivities
    states_[k].set(integrator_.output(INTEGRATOR_XF));
    
    for(int dir=0; dir<nfdir; ++dir){ 
      // Pass the forward seed to the output function 
      output_fcn_.setFwdSeed(0.0,DAE_T,dir); 
      output_fcn_.setFwdSeed(integrator_.fwdSens(INTEGRATOR_XF,dir),DAE_X,dir); 
      output_fcn_.setFwdSeed(fwdSeed(INTEGRATOR_P,dir),DAE_P,dir); 
    }

    // Evaluate output function
    output_fcn_.evaluate(nfdir);

    // Save the output of the function
    for(int i=0; i<output_.size(); ++i){
      const Matrix<double> &res = output_fcn_.output(i);
      Matrix<double> &ores = output(i);
      for(int j=0; j<res.numel(); ++j){
        ores(k,j) = res(j); // NOTE: inefficient implementation
      }
    
      // Save the forward sensitivities
      for(int dir=0; dir<nfdir; ++dir){
        const Matrix<double> &fres = output_fcn_.fwdSens(i,dir); 
        Matrix<double> &ofres = fwdSens(i,dir); 
        for(int j=0; j<fres.numel(); ++j){ 
          ofres(k,j) = fres(j); // NOTE: inefficient implementation 
        }
      }     
    }
  }
}

void SimulatorInternal::updateNumSens(bool recursive){

  if (recursive) {
    FXInternal::updateNumSens(recursive);
  }
  
  if (!output_fcn_.isNull()) {
    output_fcn_.setOption("number_of_fwd_dir",getOption("number_of_fwd_dir"));
    output_fcn_.setOption("number_of_adj_dir",getOption("number_of_adj_dir"));
    output_fcn_.updateNumSens();
  }
    
  if (!integrator_.isNull()) {
    integrator_.setOption("number_of_fwd_dir",getOption("number_of_fwd_dir"));
    integrator_.setOption("number_of_adj_dir",getOption("number_of_adj_dir"));
    integrator_.updateNumSens();
  }
}

} // namespace CasADi


