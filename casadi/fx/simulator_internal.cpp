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

  
SimulatorInternal::SimulatorInternal(const Integrator& integrator, const FX& output_fcn, const vector<double>& grid) : integrator_(integrator), output_fcn_(output_fcn), gridr_(grid){
  
  setOption("name","unnamed simulator");
  addOption("np",OT_INTEGER,GenericType(),"The number of parameters. If this option is not set, all of input(INTEGRATOR_P) is considered static parameters. The remainder of nv = input(INTEGRATOR_P) is considered to be varying parameters.");
  addOption("nf",OT_INTEGER,1,"Number of fine grained integration steps.");
}
  
SimulatorInternal::~SimulatorInternal(){
}

void SimulatorInternal::init(){
  // Number of fine-grained steps
  nf_ = getOption("nf");
  
  casadi_assert_message(nf_>=0,"Invalid parameter. nf must be at least 1.");
  
  if (nf_==1) {
    // The default case: don't change the grid
    grid_ = gridr_;
  } else {
    // Interpolate the grid.
    grid_.resize((gridr_.size()-1)*nf_+1);
    
    std::vector< double > refined(nf_+1,0);

    for (int k=0;k<gridr_.size()-1;++k) {
      linspace(refined,gridr_[k],gridr_[k+1]);
      std::copy(refined.begin(),refined.end()-1,grid_.begin()+k*nf_);
    }
    
    grid_[grid_.size()-1] = gridr_[gridr_.size()-1];
  }
  
  // Let the integration time start from the first point of the time grid.
  if (!grid_.empty()) integrator_.setOption("t0",grid_[0]);

  // Initialize the integrator
  integrator_.init();
  
  if (hasSetOption("np")) {
    np_ = getOption("np");
    casadi_assert_message(np_<=integrator_.input(INTEGRATOR_P).size(),"Invalid parameter. np (" << np_ << ") cannot be greater that input(INTEGRATOR_P), which is of size " << integrator_.input(INTEGRATOR_P).size() << ".");
    casadi_assert_message(np_>0,"Invalid parameter. np (" << np_ << ") must be greater than zero.");
    nv_ = integrator_.input(INTEGRATOR_P).size() - np_;
  } else {
    np_ = integrator_.input(INTEGRATOR_P).size();
    nv_ = 0;
  }
  
  // Cache some ranges
  np_i = range(np_);
  nv_i = range(np_,np_+nv_);

  // Generate an output function if there is none (returns the whole state)
  if(output_fcn_.isNull()){
    SXMatrix t = ssym("t");
    SXMatrix x = ssym("x",integrator_.input(INTEGRATOR_X0).sparsity());
    SXMatrix xdot = ssym("xp",integrator_.input(INTEGRATOR_XP0).sparsity());
    SXMatrix p = ssym("p",integrator_.input(INTEGRATOR_P).sparsity());

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
  
  // Allocate inputs
  input_.resize(SIMULATOR_NUM_IN);
  input(SIMULATOR_X0)  = integrator_.input(INTEGRATOR_X0);
  input(SIMULATOR_P)   = integrator_.input(INTEGRATOR_P)[np_i];
  input(SIMULATOR_XP0) = integrator_.input(INTEGRATOR_XP0);
  input(SIMULATOR_V)   = repmat(integrator_.input(INTEGRATOR_P)[nv_i],gridr_.size()-1,1);


  // Allocate outputs
  output_.resize(output_fcn_->output_.size());
  for(int i=0; i<output_.size(); ++i)
    output(i) = Matrix<double>(grid_.size(),output_fcn_.output(i).numel(),0);

  // Call base class method
  FXInternal::init();
  
}

void SimulatorInternal::evaluate(int nfdir, int nadir){
  // Pass the parameters and initial state
  integrator_.setInput(input(SIMULATOR_X0),INTEGRATOR_X0);
  integrator_.setInput(input(SIMULATOR_XP0),INTEGRATOR_XP0);
  integrator_.input(INTEGRATOR_P)[np_i] = input(SIMULATOR_P);
  integrator_.input(INTEGRATOR_P)[nv_i] = input(SIMULATOR_V)(0,ALL);
  
  // Pass sensitivities if fsens
  for(int dir=0; dir<nfdir; ++dir){
    integrator_.setFwdSeed(fwdSeed(SIMULATOR_X0,dir),INTEGRATOR_X0,dir);
    integrator_.setFwdSeed(fwdSeed(SIMULATOR_XP0,dir),INTEGRATOR_XP0,dir);
    integrator_.fwdSeed(INTEGRATOR_P,dir)[np_i] = fwdSeed(SIMULATOR_P);
    integrator_.fwdSeed(INTEGRATOR_P,dir)[nv_i] = fwdSeed(SIMULATOR_V)(0,ALL);
  }
  
  // Reset the integrator_
  integrator_.reset(nfdir, nadir);
  
  // An index that only increments on coarse time grid steps
  int k_coarse = -1;
  
  // Advance solution in time
  for(int k=0; k<grid_.size(); ++k){
   
    if (k % nf_==0 && nv_ > 0) {
      k_coarse++;
      integrator_.input(INTEGRATOR_P)[nv_i] = input(SIMULATOR_V)(k_coarse,ALL);
      // @TODO: reset the integration somehow.
      // http://sundials.2283335.n4.nabble.com/ReInit-functions-in-sundialsTB-td3239946.html
    }

    // Integrate to the output time
    integrator_.integrate(grid_[k]);

    // Pass integrator output to the output function
    output_fcn_.setInput(grid_[k],DAE_T);
    output_fcn_.setInput(integrator_.output(INTEGRATOR_XF),DAE_Y);
    if(output_fcn_.input(DAE_YDOT).size()!=0)
      output_fcn_.setInput(integrator_.output(INTEGRATOR_XPF),DAE_YDOT);
    output_fcn_.input(DAE_P) = integrator_.input(INTEGRATOR_P);

    for(int dir=0; dir<nfdir; ++dir){
      // Pass the forward seed to the output function
      output_fcn_.setFwdSeed(0.0,DAE_T,dir);
      output_fcn_.setFwdSeed(integrator_.fwdSens(INTEGRATOR_XF,dir),DAE_Y,dir);
      if(output_fcn_.input(DAE_YDOT).size()!=0)
        output_fcn_.setFwdSeed(integrator_.fwdSens(INTEGRATOR_XPF,dir),DAE_YDOT,dir);
      output_fcn_.fwdSeed(DAE_P,dir)[np_i] = fwdSeed(SIMULATOR_P,dir);
      output_fcn_.fwdSeed(DAE_P,dir)[nv_i] = fwdSeed(SIMULATOR_V,dir)(k_coarse,ALL);
    }
      
    // Evaluate output function
    output_fcn_.evaluate(nfdir,0);

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
  
  // Adjoint sensitivities
  if(nadir>0){
    casadi_assert_message(0, "not implemented");

    // Clear the seeds
    for(int dir=0; dir<nadir; ++dir){
      integrator_.adjSeed(INTEGRATOR_XF,dir).setAll(0);
    }

    // Reset the integrator for backward integration
    integrator_.resetAdj();

    // Integrate backwards
    for(int k=grid_.size()-1; k>=0; --k){
      
      // Integrate back to the previous grid point
      integrator_.integrateAdj(grid_[k]);
      
      // HERE I NEED THE STATE VECTOR!
      
      // Pass the state to the 
      
      // Pass adjoint seeds to integrator
/*      for(int i=0; i<xfs.size(); ++i)
        xfs.at(i) = x0s.at(i) + xf_seed.at(k*xfs.size() + i);
      
      // Add the contribution to the parameter sensitivity
      for(int i=0; i<ps.size(); ++i)
        ps_sim[i] += ps[i];

      // Reset the integrator to deal with the jump in seeds
      integrator_->resetAdj();*/
    }
    
    // Save
/*    vector<double> &x0_sim = input(SIMULATOR_X0).data(1);
    copy(x0s.begin(),x0s.end(),x0_sim.begin());*/
    
  }

}

} // namespace CasADi


