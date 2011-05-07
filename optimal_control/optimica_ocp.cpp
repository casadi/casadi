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
  is_scaled_ = false;
  blt_sorted_ = false;
  t_ = SX("t");
}

void OCP::repr(ostream &stream) const{
  stream << "Optimal control problem (";
  stream << "#dae = " << dynamic_eq_.size() << ", ";
  stream << "#initial_eq_ = " << initial_eq_.size() << ", ";
  stream << "#cfcn = " << cfcn.size() << ", ";
  stream << "#mterm = " << mterm.size() << ", ";
  stream << "#lterm = " << lterm.size() << ")";
}

void OCP::print(ostream &stream) const{
  // Variables in the class hierarchy
  stream << "Variables" << endl;
  variables_.print(stream);

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
  stream << "#x = " << xd_.size() << ", ";
  stream << "#z = " << xa_.size() << ", ";
  stream << "#u = " << u_.size() << ", ";
  stream << "#p = " << p_.size() << ", ";
  stream << endl << endl;
  
  // Print the differential-algebraic equation
  stream << "Dynamic equations" << endl;
  for(vector<SX>::const_iterator it=dynamic_eq_.begin(); it!=dynamic_eq_.end(); it++){
    stream << "0 == "<< *it << endl;
  }
  stream << endl;

  stream << "Initial equations" << endl;
  for(vector<SX>::const_iterator it=initial_eq_.begin(); it!=initial_eq_.end(); it++){
    stream << "0 == " << *it << endl;
  }
  stream << endl;

  // Print the explicit differential equations
/*  stream << "Differential equations (explicit)" << endl;
  for(vector<Variable>::const_iterator it=xd_.begin(); it!=xd_.end(); it++){
    SX de = it->rhs();
    if(!de->isNan())
      stream << "der(" << *it << ") == " << de << endl;
  }
  stream << endl;*/
  
  // Dependent equations
  stream << "Dependent equations" << endl;
  for(int i=0; i<explicit_lhs_.size(); ++i)
    stream << explicit_lhs_[i] << " == " << explicit_rhs_[i] << endl;
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

void OCP::eliminateDependent(){
  Matrix<SX> v = explicit_lhs_;
  Matrix<SX> v_old = explicit_rhs_;
  
  dynamic_eq_= substitute(dynamic_eq_,v,v_old).data();
  initial_eq_= substitute(initial_eq_,v,v_old).data();
  cfcn    = substitute(cfcn,v,v_old).data();
  cfcn_lb = substitute(cfcn_lb,v,v_old).data();
  cfcn_ub = substitute(cfcn_ub,v,v_old).data();
  mterm   = substitute(mterm,v,v_old).data();
  lterm   = substitute(lterm,v,v_old).data();
}

void OCP::addExplicitEquation(const SX& var, const SX& bind_eq){
  // Eliminate previous binding equations from the expression
  SX bind_eq_eliminated = substitute(bind_eq, explicit_lhs_, explicit_rhs_).data().front();
  
  explicit_lhs_.push_back(var);
  explicit_rhs_.push_back(bind_eq_eliminated);
}

void OCP::sortType(){
  // Get all the variables
  vector<Variable> v;
  variables_.getAll(v,true);
  
  // Clear variables
  x_.clear();
  u_.clear();
  p_.clear();
    
  // Loop over variables
  for(vector<Variable>::iterator it=v.begin(); it!=v.end(); ++it){
    // If not dependent
    if(!it->getDependent()){
      // Try to determine the type
      if(it->getVariability() == PARAMETER){
        p_.push_back(*it);
      } else if(it->getVariability() == CONTINUOUS) {
        if(it->getCausality() == INTERNAL){
          x_.push_back(*it);
        } else if(it->getCausality() == INPUT){
          u_.push_back(*it);
        }
      } else if(it->getVariability() == CONSTANT){
        cout << *it << endl;
        casadi_assert(0);
/*        d_.push_back(*it);*/
      }
    }
  }
  
  sortState();
}

void OCP::sortState(){
  xd_.clear();
  xa_.clear();
  for(vector<Variable>::const_iterator it=x_.begin(); it!=x_.end(); ++it){
    if(it->isDifferential()){
      xd_.push_back(*it);
    } else {
      xa_.push_back(*it);
    }
  }
}

void OCP::scale(){
  /// Make sure that the OCP has not already been scaled
  casadi_assert(!is_scaled_);
  
  
  // Variables
  Matrix<SX> t = t_;
  Matrix<SX> x = var(xd_);
  Matrix<SX> xdot = der(xd_);
  Matrix<SX> z = var(xa_);
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
  Matrix<SX> t_n = 1.;
  Matrix<SX> x_n = getNominal(xd_);
  Matrix<SX> xdot_n = getNominal(xd_);
  Matrix<SX> z_n = getNominal(xa_);
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
  explicit_rhs_= substitute(explicit_rhs_,v,v_old).data();
  dynamic_eq_= substitute(dynamic_eq_,v,v_old).data();
  initial_eq_= substitute(initial_eq_,v,v_old).data();
  cfcn    = substitute(cfcn,v,v_old).data();
  cfcn_lb = substitute(cfcn_lb,v,v_old).data();
  cfcn_ub = substitute(cfcn_ub,v,v_old).data();
  mterm   = substitute(mterm,v,v_old).data();
  lterm   = substitute(lterm,v,v_old).data();
  
  is_scaled_ = true;
}

void OCP::sortBLT(){
  // Unknown variables in the dynamic equations
  vector<SX> v;
  v.reserve(x_.size());
  for(vector<Variable>::const_iterator it=x_.begin(); it!=x_.end(); ++it)
    v.push_back(it->highest());

  // Create Jacobian in order to find the sparsity
  SXFunction fcn(v,dynamic_eq_);
  Matrix<SX> J = fcn.jac();
  
  // BLT transformation
  Interfaces::BLT blt(J.sparsity());

  // Permute equations
  vector<SX> dynamic_eq_old = dynamic_eq_;
  for(int i=0; i<dynamic_eq_.size(); ++i){
    dynamic_eq_[i] = dynamic_eq_old[blt.rowperm[i]];
  }
  
  // Permute variables
  vector<Variable> x_old = x_;
  for(int i=0; i<x_.size(); ++i){
    x_[i] = x_old[blt.colperm[i]];
  }
  sortState();
  
  // Save blocks
  rowblock_ = blt.rowblock;
  colblock_ = blt.colblock;
  
  blt_sorted_ = true;

}

void OCP::makeExplicit(){
  casadi_assert_message(blt_sorted_,"OCP has not been BLT sorted, call sortBLT()");

  // Unknown variables in the dynamic equations
  vector<SX> v;
  v.reserve(x_.size());
  for(vector<Variable>::const_iterator it=x_.begin(); it!=x_.end(); ++it)
    v.push_back(it->highest());

  // Create Jacobian
  SXFunction fcn(v,dynamic_eq_);
  SXMatrix J = fcn.jac();
  // J.printDense();

  // Cumulative variables and definitions
  SXMatrix vb_cum;
  SXMatrix def_cum;

  // Block variables and equations
  vector<SX> vb, fb;
  
  // Loop over blocks
  int nb = rowblock_.size();
  for(int b=0; b<nb; ++b){
    // Block size
    int bs = rowblock_[b+1] - rowblock_[b];
    
    // Get local variables
    vb.clear();
    for(int i=colblock_[b]; i<colblock_[b+1]; ++i)
      vb.push_back(v[i]);

    // Get local equations
    fb.clear();
    for(int i=rowblock_[b]; i<rowblock_[b+1]; ++i)
      fb.push_back(dynamic_eq_[i]);

    // Get local Jacobian
    SXMatrix Jb = J[range(rowblock_[b],rowblock_[b+1]),range(colblock_[b],colblock_[b+1])];
   
    if(dependsOn(Jb,vb)){
      // Cannot solve for vb, add to list of implicit equations
      casadi_assert_message(0,"Not implemented");
/*      append(vb_cum,SXMatrix(vb));
      append(def_cum,SXMatrix(vb));*/
      
    } else {
      // Divide fb into a part which depends on vb and a part which doesn't according to "fb == prod(Jb,vb) + fb_res"
      SXMatrix fb_res = substitute(fb,vb,SXMatrix(vb.size(),1,0));
      SXMatrix fb_exp;
      
      // Solve for vb
      if (bs <= 3){
        // Calculate inverse and multiply for very small matrices
        fb_exp = prod(inv(Jb),-fb_res);
      } else {
        // QR factorization
        fb_exp = solve(Jb,-fb_res);
      }
        
      // Substitute variables that have already been defined
      if(!def_cum.empty()){
        fb_exp = substitute(fb_exp,vb_cum,def_cum);
      }
      
      casadi_assert(0);
      append(vb_cum,SXMatrix(vb));
      append(def_cum,SXMatrix(fb_exp));
    }
  }
  
  /*  // Sort the variables
  sortType();

  // Dynamic equation
  SXMatrix dae(this->dae);
  SXMatrix xdot(xd_.size(),1,0);
  for(int i=0; i<xd_.size(); ++i)
    xdot[i] = xd_[i].der();
  
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
  for(int i=0; i<xd_.size(); ++i)
    xd_[i].setEquation(xd_[i].der(),rhs(i,0));*/
}

void OCP::makeSemiExplicit(){
  
  
//     fcn = SXFunction([v_new],[dae_new])
//   J = fcn.jac()
//   #J.printDense()
// 
//   # Cumulative variables and definitions
//   vb_cum = SXMatrix()
//   def_cum = SXMatrix()
// 
//   for b in range(blt.nb):
//     Jb = J[blt.rowblock[b]:blt.rowblock[b+1],blt.colblock[b]:blt.colblock[b+1]]
//     vb = v_new[blt.colblock[b]:blt.colblock[b+1]]
//     fb = dae_new[blt.rowblock[b]:blt.rowblock[b+1]]
//     
//     # Block size
//     bs = blt.rowblock[b+1] - blt.rowblock[b]
// 
//     #print "block ", b,
// 
//     if dependsOn(Jb,vb):
//       # Cannot solve for vb, add to list of implicit equations
//       raise Exception("Not implemented")
//       vb_cum.append(vb)
//       def_cum.append(vb)
//       
//     else:
//       # Divide fb into a part which depends on vb and a part which doesn't according to "fb == prod(Jb,vb) + fb_res"
//       fb_res = substitute(fb,vb,SXMatrix(len(vb),1,SX(0)))
//       
//       # Solve for vb
//       if bs <= 3:
//         # Calculate inverse and multiply for very small matrices
//         fb_exp = dot(inv(Jb),-fb_res)
//       else:
//         # QR factorization
//         fb_exp = solve(Jb,-fb_res)
//         
//       # Substitute variables that have already been defined
//       if not def_cum.empty():
//         fb_exp = substitute(fb_exp,vb_cum,def_cum)
//       
//       append(vb_cum,vb)
//       append(def_cum,fb_exp)

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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


VariableTree& VariableTree::subByName(const string& name, bool allocate){
  // try to locate the variable
  map<string, int>::iterator it = name_part_.find(name);

  // check if the variable exists
  if(it==name_part_.end()){
    // Allocate variable
    if(!allocate){
      stringstream ss;
      ss << "No child \"" << name << "\"";
      throw CasadiException(ss.str());
    }
    children_.push_back(VariableTree());
    name_part_[name] = children_.size()-1;
    return children_.back();
  } else {
    // Variable exists
    return children_.at(it->second);  
  }
}

VariableTree& VariableTree::subByIndex(int ind, bool allocate){
  // Modelica is 1-based
  const int base = 1; 
  
  casadi_assert(ind-base>=0);
  if(ind-base<children_.size()){
    // VariableTree exists
    return children_[ind-base];
  } else {
    // Create VariableTree
    if(!allocate){
      stringstream ss;
      ss << "Index [" << ind << "] out of bounds";
      throw CasadiException(ss.str());
    }
    children_.resize(ind-base+1);
    return children_.back();
  }
}

void VariableTree::getAll(std::vector<Variable>& v, bool skip_dependent) const{
  /// Add variables
  for(vector<VariableTree>::const_iterator it=children_.begin(); it!=children_.end(); ++it){
    it->getAll(v);
  }
  
  /// Add variable, if any
  if(!var_.isNull() && !(skip_dependent && var_.getDependent()))
    v.push_back(var_);
}

void VariableTree::print(ostream &stream, int indent) const{
  // Print variable
  if(!var_.isNull()){
    for(int i=0; i<indent; ++i) stream << " ";
    stream << var_ << endl;
  }
 
  for(std::map<std::string,int>::const_iterator it=name_part_.begin(); it!=name_part_.end(); ++it){
    cout << it->first << ", " << it->second << endl;
  }
 
  
           ;

  
  for(vector<VariableTree>::const_iterator it = children_.begin(); it!=children_.end(); ++it){
    it->print(stream,indent+2);
  }
}

  } // namespace OptimalControl
} // namespace CasADi

