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

using namespace std;
namespace CasADi{
  namespace OptimalControl{

OCP::OCP(){
  variables = Variable("variables");
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
  OCPVariables var(variables);
  stream << var << endl;
  
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
  for(vector<Variable>::const_iterator it=var.x.begin(); it!=var.x.end(); it++){
    SX de = it->getDifferentialEquation();
    if(!de->isNan())
      stream << "der(" << *it << ") == " << de << endl;
  }
  stream << endl;
  
  // Dependent equations
  stream << "Dependent equations" << endl;
  for(vector<Variable>::const_iterator it=var.d.begin(); it!=var.d.end(); it++)
    stream << *it << " == " << it->getBindingEquation() << endl;
  stream << endl;

  // Mayer terms
  stream << "Mayer objective terms" << endl;
  for(int i=0; i<mterm.size(); ++i)
    stream << mterm[i] << " at time == " << mtp[i] << endl;
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

OCP OCP::scale() const{
  // Return object
  OCP ret(*this);

  // Sort the variables according to type
  OCPVariables var(variables);
  
  // Variables
  Matrix<SX> t = sx(var.t);
  Matrix<SX> x = sx(var.x);
  Matrix<SX> xdot = der(var.x);
  Matrix<SX> z = sx(var.z);
  Matrix<SX> p = sx(var.p);
  Matrix<SX> u = sx(var.u);
  
  // Get all the variables
  Matrix<SX> v;
  append(v,t);
  append(v,x);
  append(v,xdot);
  append(v,z);
  append(v,p);
  append(v,u);
  
  // Nominal values
  Matrix<SX> t_n = getNominal(var.t);
  Matrix<SX> x_n = getNominal(var.x);
  Matrix<SX> xdot_n = getNominal(var.x);
  Matrix<SX> z_n = getNominal(var.z);
  Matrix<SX> p_n = getNominal(var.p);
  Matrix<SX> u_n = getNominal(var.u);
  
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
  ret.dae     = substitute(dae,v,v_old);
  ret.initeq  = substitute(initeq,v,v_old);
  ret.cfcn    = substitute(cfcn,v,v_old);
  ret.cfcn_lb = substitute(cfcn_lb,v,v_old);
  ret.cfcn_ub = substitute(cfcn_ub,v,v_old);
  ret.mterm   = substitute(mterm,v,v_old);
  ret.lterm   = substitute(lterm,v,v_old);
  return ret;
}

void OCP::makeExplicit(){
  // Sort the variables
  OCPVariables var(variables);

  // Dynamic equation
  SXMatrix dae(this->dae);
  SXMatrix xdot(var.x.size(),1,0);
  for(int i=0; i<var.x.size(); ++i)
    xdot[i] = var.x[i].getDerivative();
  
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
  for(int i=0; i<var.x.size(); ++i)
    var.x[i].setDifferentialEquation(rhs(i,0));
  
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

