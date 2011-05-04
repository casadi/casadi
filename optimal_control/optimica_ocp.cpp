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

#include "optimica_ocp.hpp"
#include <algorithm>
#include <set>

#include "../casadi/casadi_exception.hpp"
#include "../casadi/stl_vector_tools.hpp"
#include "variable_tools.hpp"
#include "../casadi/matrix/matrix_tools.hpp"
#include "../casadi/sx/sx_tools.hpp"
#include "../interfaces/csparse/csparse_tools.hpp"

using namespace std;
namespace CasADi{
  namespace OptimalControl{

OCP::OCP(){
  variables = Variable("variables",false);
  is_scaled_ = false;
}

void OCP::repr(ostream &stream) const{
  stream << "Optimal control problem (";
  stream << "#dae = " << dae.size() << ", ";
  stream << "#initeq = " << initeq.size() << ", ";
  stream << "#cfcn = " << cfcn.size() << ", ";
  stream << "#mterm = " << mterm.size() << ", ";
  stream << "#lterm = " << lterm.size() << ")";
}

void OCP::print(ostream &stream) const{
  // Variables in the class hierarchy
  stream << "Variables" << endl;
  stream << variables << endl;

  // Print the variables
/*  stream << "{" << endl;
  stream << "  t = " << t << endl;
  stream << "  x =  " << x << endl;
  stream << "  z =  " << z << endl;
  stream << "  u =  " << u << endl;
  stream << "  p =  " << p << endl;
  stream << "  c =  " << c << endl;
  stream << "  d =  " << d << endl;
  stream << "}" << endl;*/
  stream << "Dimensions: "; 
  stream << "#x = " << x_.size() << ", ";
  stream << "#z = " << z_.size() << ", ";
  stream << "#u = " << u_.size() << ", ";
  stream << "#p = " << p_.size() << ", ";
  stream << "#c = " << c_.size() << ", ";
  stream << "#d = " << d_.size() << ")";
  
  // Print the differential-algebraic equation
  stream << "Differential-Algebraic Equations" << endl;
  for(vector<SX>::const_iterator it=dae.begin(); it!=dae.end(); it++){
    stream << "0 == "<< *it << endl;
  }
  stream << endl;

  stream << "Initial equations" << endl;
  for(vector<SX>::const_iterator it=initeq.begin(); it!=initeq.end(); it++){
    stream << "0 == " << *it << endl;
  }
  stream << endl;

  // Print the explicit differential equations
  stream << "Differential equations (explicit)" << endl;
  for(vector<Variable>::const_iterator it=x_.begin(); it!=x_.end(); it++){
    SX de = it->rhs();
    if(!de->isNan())
      stream << "der(" << *it << ") == " << de << endl;
  }
  stream << endl;
  
  // Dependent equations
  stream << "Dependent equations" << endl;
  for(vector<Variable>::const_iterator it=d_.begin(); it!=d_.end(); it++)
    stream << *it << " == " << it->rhs() << endl;
  stream << endl;

  // Mayer terms
  stream << "Mayer objective terms" << endl;
  #ifdef WITH_TIMEDVARIABLE
  for(int i=0; i<mterm.size(); ++i)
    stream << mterm[i] << endl;
  #else // WITH_TIMEDVARIABLE
  for(int i=0; i<mterm.size(); ++i)
    stream << mterm[i] << " at time == " << mtp[i] << endl;
  #endif // WITH_TIMEDVARIABLE
  stream << endl;
  
  // Lagrange terms
  stream << "Lagrange objective terms" << endl;
  for(int i=0; i<lterm.size(); ++i)
    stream << lterm[i] << endl;
  stream << endl;
  
  // Constraint functions
  stream << "Constraint functions" << endl;
  for(int i=0; i<cfcn.size(); ++i)
    stream << cfcn_lb[i] << " <= " << cfcn[i] << " <= " << cfcn_ub[i] << endl;
  stream << endl;
  
  // Constraint functions
  stream << "Time horizon" << endl;
  stream << "t0 = " << t0 << endl;
  stream << "tf = " << tf << endl;
  
}

void OCP::sortType(){
  // Get all the variables
  vector<Variable> v = variables;
  
  // Clear variables
  t_ = Variable();
  x_.clear();
  z_.clear();
  u_.clear();
  p_.clear();
  c_.clear();
  d_.clear();
  
  // Loop over variables
  for(vector<Variable>::iterator it=v.begin(); it!=v.end(); ++it){
    // Make sure that the variable is initialized
    switch(it->getType()){
      case TYPE_INDEPENDENT:        casadi_assert(t_.isNull());     t_ = *it;  break;
      case TYPE_STATE:              x_.push_back(*it);  break;
      case TYPE_ALGEBRAIC:          z_.push_back(*it);  break;
      case TYPE_CONTROL:            u_.push_back(*it);  break;
      case TYPE_PARAMETER:          p_.push_back(*it);  break;
      case TYPE_CONSTANT:           c_.push_back(*it);  break;
      case TYPE_DEPENDENT:          d_.push_back(*it);  break;
      default: throw CasadiException("OCP::sortVariables: unknown type for " + it->getName());
    }
  }
}

void OCP::scale(){
  /// Make sure that the OCP has not already been scaled
  casadi_assert(!is_scaled_);
  
  // Sort the variables according to type
  sortType();
  
  // Variables
  Matrix<SX> t = var(t_);
  Matrix<SX> x = var(x_);
  Matrix<SX> xdot = der(x_);
  Matrix<SX> z = var(z_);
  Matrix<SX> p = var(p_);
  Matrix<SX> u = var(u_);
  
  // Get all the variables
  Matrix<SX> v;
  append(v,t);
  append(v,x);
  append(v,xdot);
  append(v,z);
  append(v,p);
  append(v,u);
  
  // Nominal values
  Matrix<SX> t_n = getNominal(t_);
  Matrix<SX> x_n = getNominal(x_);
  Matrix<SX> xdot_n = getNominal(x_);
  Matrix<SX> z_n = getNominal(z_);
  Matrix<SX> p_n = getNominal(p_);
  Matrix<SX> u_n = getNominal(u_);
  
  // Get all the old variables in expressed in the nominal ones
  Matrix<SX> v_old;
  append(v_old,t*t_n);
  append(v_old,x*x_n);
  append(v_old,xdot*xdot_n);
  append(v_old,z*z_n);
  append(v_old,p*p_n);
  append(v_old,u*u_n);
  
  // Temporary variable
  Matrix<SX> temp;

  // Substitute equations
  dae     = substitute(dae,v,v_old).data();
  initeq  = substitute(initeq,v,v_old).data();
  cfcn    = substitute(cfcn,v,v_old).data();
  cfcn_lb = substitute(cfcn_lb,v,v_old).data();
  cfcn_ub = substitute(cfcn_ub,v,v_old).data();
  mterm   = substitute(mterm,v,v_old).data();
  lterm   = substitute(lterm,v,v_old).data();
  
  is_scaled_ = true;
}

void OCP::sortBLT(){
  // State derivatives and algebraic variables
  vector<SX> xdot = der(x_);
  vector<SX> z = var(z_);

  // Create Jacobian in order to find the sparsity
  vector<SX> v;
  v.insert(v.end(),xdot.begin(),xdot.end());
  v.insert(v.end(),z.begin(),z.end());
  SXFunction fcn(v,dae);
  Matrix<SX> J = fcn.jac();
  
  // BLT transformation
  Interfaces::BLT blt(J.sparsity());

  // Permute equations
  vector<SX> dae_old = dae;
  for(int i=0; i<dae.size(); ++i){
    dae[i] = dae[blt.rowperm[i]];
  }
  
  // Permute variables
  for(int i=0; i<v.size(); ++i){
    int j = blt.colperm[i];
    Variable& vj = j<xdot.size() ? x_[j] : z_[j-xdot.size()];
    vj.setIndex(i);
    vj.setEquation(0, dae[i]);
  }
}

void OCP::makeExplicit(){
  // Sort the variables
  sortType();

  // Dynamic equation
  SXMatrix dae(this->dae);
  SXMatrix xdot(x_.size(),1,0);
  for(int i=0; i<x_.size(); ++i)
    xdot[i] = x_[i].der();
  
  // Take the Jacobian of the ode with respect to xdot
  SXMatrix J = jacobian(dae,xdot);
  
  // Make sure that J is invertable
  casadi_assert_message(!dependsOn(J,xdot),"OCP::makeExplicit:: Dynamic equation not affine in state derivative");
  //casadi_assert_message(!det(J).isZero(),"OCP::makeExplicit: Jacobian not invertable");

  // Write the differential equation in explicit form
  SXMatrix rhs = solve(J,prod(J,xdot)-dae);

  // Remove the xdots
  rhs = substitute(rhs,xdot,SXMatrix(xdot.sparsity(),0));
  
  // Save as explicit derivative
  for(int i=0; i<x_.size(); ++i)
    x_[i].setEquation(x_[i].der(),rhs(i,0));
}

void OCP::makeSemiExplicit(){
  throw CasadiException("OCP::makeSemiExplicit: Commented out");
#if 0  
  // Move the fully implicit dynamic equations to the list of algebraic equations
  algeq.insert(algeq.end(), dyneq.begin(), dyneq.end());
  dyneq.clear();
    
  // Introduce new explicit differential equations describing the relation between states and state derivatives
  xd.insert(xd.end(), x.begin(), x.end());
  diffeq.insert(diffeq.end(), xdot.begin(), xdot.end());
  
  // Put the state derivatives in the algebraic state category (inefficient!!!)
  xa.insert(xa.end(), xdot.begin(), xdot.end());

  // Remove from old location
  xdot.clear();
  x.clear();
#endif
}

  

  } // namespace OptimalControl
} // namespace CasADi

