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

#include "integrator_internal.hpp"
#include <cassert>
#include "../stl_vector_tools.hpp"
#include "jacobian.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "mx_function.hpp"
#include "sx_function.hpp"

INPUTSCHEME(IntegratorInput)
OUTPUTSCHEME(IntegratorOutput)

using namespace std;
namespace CasADi{

IntegratorInternal::IntegratorInternal(const FX& f, const FX& q) : f_(f), q_(q){
  // set default options
  setOption("name","unnamed_integrator"); // name of the function 
  
  // Additional options
  addOption("print_stats",                 OT_BOOLEAN,  false, "Print out statistics after integration");
  addOption("nrhs",                        OT_INTEGER, 1); // number of right hand sides
  addOption("t0",                          OT_REAL, 0.0); // start of the integration
  addOption("tf",                          OT_REAL, 1.0); // end of the integration
  
  nx_ = 0;
  np_ = 0;
  
}

IntegratorInternal::~IntegratorInternal(){ 
}

void IntegratorInternal::setDimensions(int nxd, int nxa, int nxq, int nyd, int nya, int nyq, int np){
  // Save dimensions
  nxd_ = nxd;
  nxa_ = nxa;
  nxq_ = nxq;
  nyd_ = nyd;
  nya_ = nya;
  nyq_ = nyq;
  np_  = np;
  
  // Allocate space for inputs
  input_.resize(NEW_INTEGRATOR_NUM_IN);
  input(NEW_INTEGRATOR_XD0) = DMatrix::zeros(nxd_,1);
  input(NEW_INTEGRATOR_XQ0) = DMatrix::zeros(nxq_,1);
  input(NEW_INTEGRATOR_XA0) = DMatrix::zeros(nxa_,1);
  input(NEW_INTEGRATOR_P) = DMatrix::zeros(np_,1);
  
  // Allocate space for outputs
  output_.resize(NEW_INTEGRATOR_NUM_OUT);
  output(NEW_INTEGRATOR_XDF) = input(NEW_INTEGRATOR_XD0);
  output(NEW_INTEGRATOR_XQF) = input(NEW_INTEGRATOR_XQ0);
  output(NEW_INTEGRATOR_XAF) = input(NEW_INTEGRATOR_XA0);
  output(NEW_INTEGRATOR_YD0) = DMatrix::zeros(nyd_,1);
  output(NEW_INTEGRATOR_YQ0) = DMatrix::zeros(nyq_,1);
  output(NEW_INTEGRATOR_YA0) = DMatrix::zeros(nya_,1);
}

void IntegratorInternal::setDimensions(int nx, int np){
  nx_ = nx;
  np_ = np;
  
  // Allocate space for inputs
  input_.resize(INTEGRATOR_NUM_IN);
  input(INTEGRATOR_X0)  = DMatrix(nx_,1,0); // initial state value
  input(INTEGRATOR_XP0) = DMatrix(nx_,1,0); // initial state derivative value
  input(INTEGRATOR_P)   = DMatrix(np_,1,0); // parameter
  
  // Allocate space for outputs
  output_.resize(INTEGRATOR_NUM_OUT);
  output(INTEGRATOR_XF) = DMatrix(nx_,1,0);
  output(INTEGRATOR_XPF)= DMatrix(nx_,1,0);
}

void IntegratorInternal::evaluate(int nfdir, int nadir){
  
  // Reset solver
  reset(nfdir, nadir);

  // Advance solution in time
  integrate(tf_);

  if(nadir>0){
    // Re-initialize backward problem
    resetAdj();

    // Integrate backwards to the beginning
    integrateAdj(t0_);
  }
  
  // Print statistics
  if(getOption("print_stats")) printStats(std::cout);
}

void IntegratorInternal::init(){
  // Call the base class method
  FXInternal::init();

  // read options
  nrhs_ = getOption("nrhs");
  
  // Give an intial value for the time horizon
  t0_ = getOption("t0");
  tf_ = getOption("tf");
}

void IntegratorInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
  FXInternal::deepCopyMembers(already_copied);
  f_ = deepcopy(f_,already_copied);
  q_ = deepcopy(q_,already_copied);
}

} // namespace CasADi


