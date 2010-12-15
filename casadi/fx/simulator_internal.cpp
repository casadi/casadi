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
#include <cassert>

using namespace std;
namespace CasADi{
  
SimulatorInternal::SimulatorInternal(const Integrator& integrator, const FX& output_fcn, const vector<double>& grid) :
  integrator_(integrator), output_fcn_(output_fcn), grid_(grid){

  setOption("name","unnamed simulator");
      
  // Allocate inputs
  input_.resize(SIMULATOR_NUM_IN);
  input_.at(SIMULATOR_X0).setSize(integrator_->nx_);
  input_.at(SIMULATOR_P).setSize(integrator_->np_);
  
  // Allocate outputs
  output_.resize(output_fcn_->output_.size());
  for(int i=0; i<output_.size(); ++i)
    output_.at(i).setSize(grid_.size(),output_fcn_->output_.at(i).numel());
}
  
SimulatorInternal::~SimulatorInternal(){ 
}

void SimulatorInternal::init(){
  // Call base class method
  FXNode::init();
  
  // Initialize the integrator_ and output functions
  integrator_.init();
  output_fcn_.init();
}

void SimulatorInternal::evaluate(int fsens_order, int asens_order){

  // Pass the parameters and initial state
  integrator_.input(INTEGRATOR_X0).set(input(SIMULATOR_X0).data());
  integrator_.input(INTEGRATOR_P).set(input(SIMULATOR_P).data());

  // Pass sensitivities if fsens
  if(fsens_order>0){
    integrator_.input(INTEGRATOR_X0).setF(input(SIMULATOR_X0).dataF());
    integrator_.input(INTEGRATOR_P).setF(input(SIMULATOR_P).dataF());
  }
  
  // Reset the integrator_
  integrator_.reset(fsens_order, asens_order);
  
  // Advance solution in time
  for(int k=0; k<grid_.size(); ++k){
    
    // Integrate to the output time
    integrator_.integrate(grid_[k]);
    
    // Pass integrator_ output to the output function
    output_fcn_.input(OUTPUT_T).set(grid_[k]);
    output_fcn_.input(OUTPUT_X).set(integrator_.output().data());
    output_fcn_.input(OUTPUT_P).set(input(SIMULATOR_P).data());

    // Evaluate output function
    output_fcn_.evaluate();
    
    // Save the output of the function
    for(int i=0; i<output_.size(); ++i){
      const vector<double> &res = output_fcn_.output(i).data();
      copy(res.begin(),res.end(),&output(i).data()[k*res.size()]);
    }
    
    if(fsens_order>0){
      
      // Pass the forward seed to the output function
      output_fcn_.input(OUTPUT_X).setF(integrator_.output().dataF());
      output_fcn_.input(OUTPUT_P).setF(input(SIMULATOR_P).dataF());
      
      // Evaluate output function
      output_fcn_.evaluate(1,0);

      // Save the output of the function
      for(int i=0; i<output_.size(); ++i){
        const vector<double> &res = output_fcn_.output(i).dataF();
        copy(res.begin(),res.end(),&output(i).dataF()[k*res.size()]);
      }
    }
  }
  
  // Adjoint sensitivities
  if(asens_order>0){

    #if 0
          // Clear the seeds (TODO: change this when XF is included as output of the simulator!)
    vector<double> &xfs = integrator_.output(INTEGRATOR_XF).aseed();
    fill(xfs.begin(),xfs.end(),0);
    vector<double> &x0s = integrator_.input(INTEGRATOR_X0,1);
    fill(x0s.begin(),x0s.end(),0);
    vector<double> &ps = integrator_.input(INTEGRATOR_P,1);
    fill(ps.begin(),ps.end(),0);
    vector<double> &ps_sim = input(SIMULATOR_P).data(1);
    fill(ps_sim.begin(),ps_sim.end(),0);
    
    // Reset the integrator for backward integration
    integrator_->resetAdj();

    // Seeds from the output function
    const vector<double> &xf_seed = output(0).data(1); // TODO: output is here assumed to be trivial, returning the state
    assert(xf_seed.size() == grid_.size()*xfs.size());

    // Integrate backwards
    for(int k=grid_.size()-1; k>=0; --k){
      // Integrate back to the previous grid point
      integrator_->integrateAdj(grid_[k]);
      
      // Pass adjoint seeds to integrator
      for(int i=0; i<xfs.size(); ++i)
        xfs.at(i) = x0s.at(i) + xf_seed.at(k*xfs.size() + i);
      
      // Add the contribution to the parameter sensitivity
      for(int i=0; i<ps.size(); ++i)
        ps_sim[i] += ps[i];

      // Reset the integrator to deal with the jump in seeds
      integrator_->resetAdj();
    }
    
    // Save
    vector<double> &x0_sim = input(SIMULATOR_X0).data(1);
    copy(x0s.begin(),x0s.end(),x0_sim.begin());

    #endif
    
  }

}

} // namespace CasADi


