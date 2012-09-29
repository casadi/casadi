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

#include "liftopt_internal.hpp"
using namespace std;

namespace CasADi{

LiftoptSolver::LiftoptSolver(){ 
}

LiftoptSolver::LiftoptSolver(const FX& F__, const FX& G__){
  // This only works if the arguments are MXFunctions
  MXFunction F = shared_cast<MXFunction>(F__);
  MXFunction G = shared_cast<MXFunction>(G__);
  casadi_assert_message(F.isNull()==F__.isNull(), "not a MX function");
  casadi_assert_message(G.isNull()==G__.isNull(), "not a MX function");
  
  // Get the variable
  MX U = F.inputMX();
  MX U1 = G.inputMX();
  casadi_assert(U.get()==U1.get());
  
  // Get objective and constraints
  MX f = F.outputMX();
  MX g = G.outputMX();
  
  // Lagrange multipliers
  MX lam("lambda",g.size());
  
  // Lagrangian
  MX lag = f - inner_prod(lam,g);

  // NLP cost and objective function argument
  vector<MX> fcn_in(LO_NUM_IN);
  fcn_in[LO_U] = U;
  fcn_in[LO_LAMBDA] = lam;
  
  // NLP cost and objective function
  vector<MX> fcn_out(LO_NUM_OUT);
  fcn_out[LO_OBJRES] = f;
  fcn_out[LO_EQ] = g;
  fcn_out[LO_INEQ] = MX(0,0); // FIXME: MX should be a 0-by-0 matrix by default, or?
  fcn_out[LO_OBJ] = f;
  fcn_out[LO_LAGFCN] = lag;
  MXFunction fcn(fcn_in,fcn_out);

  assignNode(new LiftoptInternal(fcn));
}

LiftoptSolver::LiftoptSolver(const MXFunction& fcn){
  assignNode(new LiftoptInternal(fcn));
}

LiftoptInternal* LiftoptSolver::operator->(){
  return (LiftoptInternal*)(FX::operator->());
}

const LiftoptInternal* LiftoptSolver::operator->() const{
  return (const LiftoptInternal*)(FX::operator->());
}

bool LiftoptSolver::checkNode() const{
  return dynamic_cast<const LiftoptInternal*>(get());
}

std::vector<double>& LiftoptSolver::nodeInit(){
  return (*this)->nodeInit;
}

const std::vector<double>& LiftoptSolver::nodeInit() const{
  return (*this)->nodeInit;
}

} // namespace CasADi

