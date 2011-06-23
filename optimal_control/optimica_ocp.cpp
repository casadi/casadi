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
#include <ctime>

#include "../casadi/casadi_exception.hpp"
#include "../casadi/pre_c99_support.hpp"
#include "../casadi/stl_vector_tools.hpp"
#include "variable_tools.hpp"
#include "../casadi/matrix/matrix_tools.hpp"
#include "../casadi/sx/sx_tools.hpp"
#include "../interfaces/csparse/csparse_tools.hpp"
#include "../casadi/fx/integrator.hpp"

using namespace std;
namespace CasADi{
  namespace OptimalControl{

OCP::OCP(){
  scaled_variables_ = false;
  scaled_equations_ = false;
  eliminated_dependents_ = false;
  blt_sorted_ = false;
  t_ = SX("t");
  t0 = numeric_limits<double>::quiet_NaN();
  tf = numeric_limits<double>::quiet_NaN();
}

void OCP::repr(ostream &stream) const{
  stream << "Optimal control problem (";
  stream << "#dae = " << implicit_fcn_.size() << ", ";
  stream << "#initial_eq_ = " << initial_eq_.size() << ", ";
  stream << "#path_fcn_ = " << path_fcn_.size() << ", ";
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
  stream << "#x = " << x_.size() << ", ";
  stream << "#z = " << z_.size() << ", ";
  stream << "#u = " << u_.size() << ", ";
  stream << "#p = " << p_.size() << ", ";
  stream << endl << endl;
  
  // Print the differential-algebraic equation
  stream << "Dynamic equations" << endl;
  for(vector<SX>::const_iterator it=implicit_fcn_.begin(); it!=implicit_fcn_.end(); it++){
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
  for(vector<Variable>::const_iterator it=x_.begin(); it!=x_.end(); it++){
    SX de = it->rhs();
    if(!de->isNan())
      stream << "der(" << *it << ") == " << de << endl;
  }
  stream << endl;*/
  
  // Dependent equations
  stream << "Dependent equations" << endl;
  for(int i=0; i<explicit_var_.size(); ++i)
    stream << explicit_var_[i] << " == " << explicit_fcn_[i] << endl;
  stream << endl;

  // Mayer terms
  stream << "Mayer objective terms" << endl;
  for(int i=0; i<mterm.size(); ++i)
    stream << mterm[i] << endl;
  stream << endl;
  
  // Lagrange terms
  stream << "Lagrange objective terms" << endl;
  for(int i=0; i<lterm.size(); ++i)
    stream << lterm[i] << endl;
  stream << endl;
  
  // Constraint functions
  stream << "Constraint functions" << endl;
  for(int i=0; i<path_fcn_.size(); ++i)
    stream << path_min_[i] << " <= " << path_fcn_[i] << " <= " << path_max_[i] << endl;
  stream << endl;
  
  // Constraint functions
  stream << "Time horizon" << endl;
  stream << "t0 = " << t0 << endl;
  stream << "tf = " << tf << endl;
  
}

void OCP::eliminateDependent(){
  cout << "eliminateDependent ..." << endl;
  double time1 = clock();

  Matrix<SX> v = explicit_var_;
  Matrix<SX> v_old = explicit_fcn_;
  
  implicit_fcn_= substitute(implicit_fcn_,v,v_old).data();
  initial_eq_= substitute(initial_eq_,v,v_old).data();
  path_fcn_    = substitute(path_fcn_,v,v_old).data();
  mterm   = substitute(mterm,v,v_old).data();
  lterm   = substitute(lterm,v,v_old).data();
  eliminated_dependents_ = true;

  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... eliminateDependent complete after " << dt << " seconds." << endl;
}

void OCP::addExplicitEquation(const Matrix<SX>& var, const Matrix<SX>& bind_eq, bool to_front){
  if(to_front){
    // Eliminate expression from the current binding equations
    Matrix<SX> explicit_fcn_eliminated = substitute(explicit_fcn_,var, bind_eq);

    // Save in the old location
    explicit_fcn_.swap(explicit_fcn_eliminated.data());
    
    // Add to list of equations
    explicit_var_.insert(explicit_var_.end(),var.begin(),var.end());
    explicit_fcn_.insert(explicit_fcn_.end(),bind_eq.begin(),bind_eq.end());
    
  } else {
    // Eliminate previous binding equations from the expression
    Matrix<SX> bind_eq_eliminated = substitute(bind_eq, explicit_var_, explicit_fcn_);
    
    // Add to list of equations
    explicit_var_.insert(explicit_var_.end(),var.begin(),var.data().end());
    explicit_fcn_.insert(explicit_fcn_.end(),bind_eq_eliminated.begin(),bind_eq_eliminated.data().end());
  }
}

void OCP::sortType(){
  // Get all the variables
  vector<Variable> v;
  variables_.getAll(v);
  
  // Clear variables
  x_.clear();
  z_.clear();
  u_.clear();
  p_.clear();
  d_.clear();
  
  // Mark all dependent variables
  for(vector<SX>::iterator it=explicit_var_.begin(); it!=explicit_var_.end(); ++it){
    it->setTemp(1);
  }
  
  // Get implicit variables
  bool find_implicit = implicit_var_.size() != implicit_fcn_.size();
  if(find_implicit){
    implicit_var_.clear();
    implicit_var_.reserve(implicit_fcn_.size());
  }
  
  // Loop over variables
  for(vector<Variable>::iterator it=v.begin(); it!=v.end(); ++it){
    // If not dependent
    if(it->highest().getTemp()!=1){
      // Try to determine the type
      if(it->getVariability() == PARAMETER){
        if(it->getFree()){
          p_.push_back(*it); 
        } else {
          casadi_assert(0);
        }
      } else if(it->getVariability() == CONTINUOUS) {
        if(it->getCausality() == INTERNAL){
          if(it->isDifferential()){
            x_.push_back(*it);
          } else {
            z_.push_back(*it);
          }
          
          // Add to list of implicit variables
          if(find_implicit) implicit_var_.push_back(it->highest());

        } else if(it->getCausality() == INPUT){
          u_.push_back(*it);
        }
      } else if(it->getVariability() == CONSTANT){
        casadi_assert(0);
      }
    } else {
      d_.push_back(*it);
    }
  }

  // Unmark all dependent variables
  for(vector<SX>::iterator it=explicit_var_.begin(); it!=explicit_var_.end(); ++it){
    it->setTemp(0);
  }
  
  // Assert consistent number of equations and variables
/*  cout << implicit_var_ << endl;
  cout << implicit_fcn_ << endl;*/
  casadi_assert(implicit_var_.size() == implicit_fcn_.size());
}


void OCP::scaleVariables(){
  cout << "Scaling variables ..." << endl;
  double time1 = clock();
  
  //Make sure that the variables has not already been scaled
  casadi_assert(!scaled_variables_);

  // Variables
  Matrix<SX> t = t_;
  Matrix<SX> x = var(x_);
  Matrix<SX> xdot = der(x_);
  Matrix<SX> z = var(z_);
  Matrix<SX> p = var(p_);
  Matrix<SX> u = var(u_);
  
  // Collect all the variables
  Matrix<SX> v;
  append(v,t);
  append(v,x);
  append(v,xdot);
  append(v,z);
  append(v,p);
  append(v,u);
  
  // Nominal values
  Matrix<SX> t_n = 1.;
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
  explicit_fcn_= substitute(explicit_fcn_,v,v_old).data();
  implicit_fcn_= substitute(implicit_fcn_,v,v_old).data();
  initial_eq_= substitute(initial_eq_,v,v_old).data();
  path_fcn_    = substitute(path_fcn_,v,v_old).data();
  mterm   = substitute(mterm,v,v_old).data();
  lterm   = substitute(lterm,v,v_old).data();
  
  scaled_variables_ = true;
  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... variable scaling complete after " << dt << " seconds." << endl;
}
    
void OCP::scaleEquations(){
  
  // Make sure that the equations has not already been scaled
  casadi_assert(!scaled_equations_);
  
  // Make sure that the dependents have been eliminated
  casadi_assert(eliminated_dependents_);
  
  // Make sure that the variables have been scaled
  casadi_assert(scaled_variables_);
  
  // Quick return if no implicit equations
  if(implicit_fcn_.empty())
    return;

  cout << "Scaling equations ..." << endl;
  double time1 = clock();

  // Variables
  enum Variables{T,X,XDOT,Z,P,U,NUM_VAR};
  vector<Matrix<SX> > v(NUM_VAR); // all variables
  v[T] = t_;
  v[X] = var(x_);
  v[XDOT] = der(x_);
  v[Z] = var(z_);
  v[P] = var(p_);
  v[U] = var(u_);

  // Create the jacobian of the implicit equations with respect to [x,z,p,u] 
  Matrix<SX> xz;
  append(xz,v[X]);
  append(xz,v[Z]);
  append(xz,v[P]);
  append(xz,v[U]);
  SXFunction fcn = SXFunction(xz,implicit_fcn_);
  SXFunction J(v,fcn.jac());

  // Evaluate the Jacobian in the starting point
  J.init();
  J.setInput(0.0,T);
  J.setInput(getStart(x_,true),X);
  J.input(XDOT).setAll(0.0);
  J.setInput(getStart(z_,true),Z);
  J.setInput(getStart(p_,true),P);
  J.setInput(getStart(u_,true),U);
  J.evaluate();
  
  // Get the maximum of every row
  Matrix<double> &J0 = J.output();
  vector<double> scale(J0.size1(),0.0); // scaling factors
  for(int i=0; i<J0.size1(); ++i){
    // Loop over non-zero entries of the row
    for(int el=J0.rowind(i); el<J0.rowind(i+1); ++el){
      // Column
      //int j=J0.col(el);
      
      // The scaling factor is the maximum norm, ignoring not-a-number entries
      if(!isnan(J0.at(el))){
        scale[i] = max(scale[i],fabs(J0.at(el)));
      }
    }
    
    // Make sure 
    if(scale[i]==0){
      cout << "Warning: Could not generate a scaling factor for equation " << i << "(0 == " << implicit_fcn_[i] << "), selecting 1." << endl;
      scale[i]=1.;
    }
  }
  
  // Scale the equations
  for(int i=0; i<implicit_fcn_.size(); ++i){
    implicit_fcn_[i] /= scale[i];
  }
  
  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... equation scaling complete after " << dt << " seconds." << endl;
  scaled_equations_ = true;
}

void OCP::sortBLT(bool with_x){
  cout << "BLT sorting ..." << endl;
  double time1 = clock();
  
  // Sparsity pattern
  CRSSparsity sp;
  
  if(with_x){
    // inverse time constant
    SX invtau("invtau");

    // Replace x with invtau*xdot in order to get a Jacobian which also includes x
    SXMatrix implicit_fcn_with_x = substitute(implicit_var_,var(x_),invtau*SXMatrix(var(x_)));
    
    // Create Jacobian in order to find the sparsity
    SXFunction fcn(implicit_var_,implicit_fcn_with_x);
    Matrix<SX> J = fcn.jac();
    sp = J.sparsity();
  } else {
    // Create Jacobian in order to find the sparsity
    SXFunction fcn(implicit_var_,implicit_fcn_);
    Matrix<SX> J = fcn.jac();
    sp = J.sparsity();
  }
  
  // BLT transformation
  Interfaces::BLT blt(sp);

  // Permute equations
  vector<SX> implicit_fcn_new(implicit_fcn_.size());
  for(int i=0; i<implicit_fcn_.size(); ++i){
    implicit_fcn_new[i] = implicit_fcn_[blt.rowperm[i]];
  }
  implicit_fcn_new.swap(implicit_fcn_);
  
  // Permute variables
  vector<SX> implicit_var_new(implicit_var_.size());
  for(int i=0; i<implicit_var_.size(); ++i){
    implicit_var_new[i]= implicit_var_[blt.colperm[i]];
  }
  implicit_var_new.swap(implicit_var_);
  
  // Save blocks
  rowblock_ = blt.rowblock;
  colblock_ = blt.colblock;
  nb_ = blt.nb;
  
  blt_sorted_ = true;

  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... BLT sorting complete after " << dt << " seconds." << endl;
}

void OCP::makeExplicit(){
  casadi_assert_message(blt_sorted_,"OCP has not been BLT sorted, call sortBLT()");

  cout << "Making explicit..." << endl;
  double time1 = clock();

  
  // Create Jacobian
  SXFunction fcn(implicit_var_,implicit_fcn_);
  SXMatrix J = fcn.jac();

  // Mark the algebraic variables
  for(int i=0; i<z_.size(); ++i){
    z_[i].var().setTemp(i+1);
  }
  
  // Get initial values for all implicit variables
  vector<double> implicit_var_guess(implicit_var_.size(),0);
  for(int i=0; i<implicit_var_.size(); ++i){
    int ind = implicit_var_[i].getTemp()-1;
    if(ind>=0){
      implicit_var_guess[i] = z_[ind].getStart()/z_[ind].getNominal();
    }
  }
  
  // Unmark the algebraic variables
  for(int i=0; i<z_.size(); ++i){
    z_[i].var().setTemp(0);
  }
  
  // Block variables and equations
  vector<SX> vb, fb;

  // Save old number of explicit variables
  int old_nexp = explicit_var_.size();

  // New implicit equation and variables
  vector<SX> implicit_var_new, implicit_fcn_new;
    
  // Loop over blocks
  for(int b=0; b<nb_; ++b){
    
    // Block size
    int bs = rowblock_[b+1] - rowblock_[b];
    
    // Get local variables
    vb.clear();
    for(int i=colblock_[b]; i<colblock_[b+1]; ++i)
      vb.push_back(implicit_var_[i]);

    // Get local equations
    fb.clear();
    for(int i=rowblock_[b]; i<rowblock_[b+1]; ++i)
      fb.push_back(implicit_fcn_[i]);

    // Get local Jacobian
    SXMatrix Jb = J(range(rowblock_[b],rowblock_[b+1]),range(colblock_[b],colblock_[b+1]));
    if(dependsOn(Jb,vb)){
      
      // Guess for the solution
      SXMatrix x_k(vb.size(),1,0);
      int offset = rowblock_[b];
      for(int i=0; i<x_k.size(); ++i){
        x_k.at(i) = implicit_var_guess[offset+i];
      }
      
      // Make Newton iterations
      bool exact_newton = true;
      SXFunction newtonIter; // Newton iteration function
      
      if(exact_newton){
        // Use exact Newton (i.e. update Jacobian in every iteration)
        newtonIter = SXFunction(vb,vb-solve(Jb,SXMatrix(fb)));
      } else {
        // Evaluate the Jacobian
        SXMatrix Jb0 = substitute(Jb,vb,x_k);
        
        // Use quasi-Newton
        newtonIter = SXFunction(vb,vb-solve(Jb0,SXMatrix(fb)));
      }
      newtonIter.init();
      
      // Make Newton iterations
      const int n_newton = 3; // NOTE: the number of newton iterations is fixed!
      for(int k=0; k<n_newton; ++k){
        x_k = newtonIter.eval(x_k);
      }
      
      cout << "Using " << (exact_newton ? "an exact " : "a quasi-") << "Newton iteration with "<< n_newton << " iterations to solve block " << b << " for the " << vb.size() << " variables " << vb << endl;
      cout << "the implicit equation has " << countNodes(fb) << " nodes" << endl;
      cout << "the Newton algorithm has " << newtonIter.algorithm().size() << " nodes" << endl;
      cout << "the explicit expression has " << countNodes(x_k) << " nodes" << endl;

      // Add binding equations
      addExplicitEquation(vb,x_k);
      
      // TODO: implement tearing
      if(false){ // not ready (but make sure it compiles)
      
        // Find out which dependencies are nonlinear and which are linear
        Matrix<int> Jb_lin(Jb.sparsity());
        for(int i=0; i<Jb.size(); ++i){
          Jb_lin.at(i) = dependsOn(Jb.at(i),vb) ? 2 : 1;
        }
        
        // Loop over rows (equations)
        for(int i=0; i<Jb_lin.size1(); ++i){
          // Number of variables appearing linearly
          int n_lin = 0;
          int j_lin = -1; // index of a variable appearing linearily
          
          // Loop over non-zero elements
          for(int el=Jb_lin.rowind(i); el<Jb_lin.rowind(i+1); ++el){
            // Column (variable)
            int j = Jb_lin.col(el);
            
            // Count number of variables appearing linearly
            if(Jb_lin.at(el)==1){
              j_lin = j;
              n_lin ++;
            }
          }
          
          // Make causal if only one variable appears linearily
          if(n_lin==1){
            // 
            
          }
        }
      }
                  
      // Cannot solve for vb, add to list of implicit equations
//      implicit_var_new.insert(implicit_var_new.end(),vb.begin(),vb.end());
//      implicit_fcn_new.insert(implicit_fcn_new.end(),fb.begin(),fb.end());
      
//      cout << "added " << fb << " and " << vb << " to list of implicit equations (" << vb.size() << " equations)" << endl;
      
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

      // Add explicit equation
      addExplicitEquation(vb,fb_exp);
    }
  }

  // Update implicit equations
  implicit_var_new.swap(implicit_var_);
  implicit_fcn_new.swap(implicit_fcn_);

  // Mark the variables made explicit
  for(vector<SX>::iterator it=explicit_var_.begin()+old_nexp; it!=explicit_var_.end(); ++it){
    it->setTemp(1);
  }
  
  // New algebraic variables
  vector<Variable> z_new;

  // Loop over algebraic variables
  for(vector<Variable>::iterator it=z_.begin(); it!=z_.end(); ++it){
    // Check if marked
    if(it->var().getTemp()){
      // Make dependent
      d_.push_back(*it);
      
      // If upper or lower bounds are finite, add path constraint
      if(!isinf(it->getMin()) || !isinf(it->getMax())){
        path_fcn_.push_back(it->var());
        path_min_.push_back(it->getMin()/it->getNominal());
        path_max_.push_back(it->getMax()/it->getNominal());
      }
    } else {
      z_new.push_back(*it);
    }
  }
  
  // Update new z_
  z_.swap(z_new);

  // Unark the variables made explicit
  for(vector<SX>::iterator it=explicit_var_.begin()+old_nexp; it!=explicit_var_.end(); ++it){
    it->setTemp(0);
  }
  
  // Eliminate the dependents
  eliminateDependent();

  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... makeExplicit complete after " << dt << " seconds." << endl;
}

void OCP::createFunctions(bool create_dae, bool create_ode, bool create_quad){
  cout << "createFunctions ..." << endl;
  double time1 = clock();

  // no quad if no quadrature states
  if(lterm.empty()) create_quad=false;

  if(create_ode){
    // Create an explicit ODE
    sortBLT(false);

    // Make the OCP explicit by eliminated all the algebraic states
    makeExplicit();
    
    // The equation could be made explicit symbolically
    if(implicit_fcn_.empty()){
      casadi_assert(z_.empty());
      
      // Mark the explicit variables
      for(int i=0; i<explicit_var_.size(); ++i){
        explicit_var_[i].setTemp(i+1);
      }

      // Find the ODE rhs
      vector<SX> ode(x_.size());
      for(int i=0; i<x_.size(); ++i){
        int ind = x_[i].der().getTemp()-1;
        casadi_assert(ind>=0);
        ode[i] = explicit_fcn_[ind];
      }
  
      // Find the dependency equations
      vector<SX> dep(d_.size());
      for(int i=0; i<d_.size(); ++i){
        int ind = d_[i].var().getTemp()-1;
        casadi_assert(ind>=0);
        dep[i] = explicit_fcn_[ind];
      }
      
  
      // Unmark the explicit variables
      for(int i=0; i<explicit_var_.size(); ++i){
        explicit_var_.at(i).setTemp(0);
      }

      // Evaluate constant expressions
      Matrix<SX> ode_elim = evaluateConstants(ode);
      
      // ODE right hand side function
      vector<SXMatrix> ode_in(DAE_NUM_IN);
      ode_in[DAE_T] = t_;
      ode_in[DAE_Y] = var(x_);
      ode_in[DAE_P] = var(u_);
      oderhs_ = SXFunction(ode_in,ode_elim);

      // Dependency function
      output_fcn_ = SXFunction(ode_in,dep);
      
      // ODE quadrature function
      if(create_quad){
        quadrhs_ = SXFunction(ode_in,lterm[0]);
        
        // ODE right hand side with lterm
        append(ode_in[DAE_Y],symbolic("lterm"));
        append(ode_elim,SXMatrix(lterm[0]));
        oderhs_with_lterm_ = SXFunction(ode_in,ode_elim);
      }
      
    } else {
      cout << implicit_fcn_.size() << " implicit equations" << endl;
      
      
    // A Newton algorithm is necessary
    casadi_assert(0);
    
//     # Create an implicit function residual
//     impres_in = (DAE_NUM_IN+1) * [[]]
//     impres_in[0] = xdot+z
//     impres_in[1+DAE_T] = t
//     impres_in[1+DAE_Y] = x
//     impres_in[1+DAE_P] = u
//     impres = SXFunction(impres_in,[ocp.implicit_fcn_])
//     impres.setOption("number_of_fwd_dir",len(x)+1)
// 
//     impres.init()
//     impres.setInput(xdot0+z0, 0)
//     impres.setInput(0., 1+DAE_T)
//     impres.setInput(x0, 1+DAE_Y)
//     impres.setInput(u0, 1+DAE_P)
//     impres.evaluate()
// 
//     # Create an implicit function (KINSOL)
//     impsolver = KinsolSolver(impres)
//     linsol = CSparse(CRSSparsity())
//     impsolver.setLinearSolver(linsol)
// 
//     impsolver.setOption("linear_solver","user_defined")
//     impsolver.setOption("max_krylov",100)
//     impsolver.setOption("number_of_fwd_dir",len(x)+1)
//     impsolver.setOption("number_of_adj_dir",0)
//     impsolver.init()
// 
//     impsolver.setInput(0., DAE_T)
//     impsolver.setInput(x0, DAE_Y)
//     impsolver.setInput(u0, DAE_P)
//     impsolver.output().set(xdot0+z0)
//     impsolver.evaluate(1,0)
//     #print "implicit function residual", impsolver.output()
// 
//     # Create the ODE residual
//     ode_in_mx = DAE_NUM_IN * [[]]
//     ode_in_mx[DAE_T] = MX("T")
//     ode_in_mx[DAE_Y] = MX("X",len(x))
//     ode_in_mx[DAE_P] = MX("U",len(u))
//     [xdot_z] = impsolver.call(ode_in_mx)
//     ode = MXFunction(ode_in_mx,[xdot_z[0:len(x)]])
//     ode.setOption("number_of_fwd_dir",len(x)+1)
// 
//     # ODE quadrature function
//     ode_in = DAE_NUM_IN * [[]]
//     ode_in[DAE_T] = t
//     ode_in[DAE_Y] = x
//     ode_in[DAE_P] = u
//     ode_lterm = SXFunction(ode_in,[[ocp.lterm[0]/1e3]])
//     ode_lterm.setOption("number_of_fwd_dir",len(x)+1)
      
    }
  } else if(create_dae) {
    // BLT sorting
    sortBLT(true);

    // Scale the equations
    scaleEquations();
    
    // DAE residual arguments
    vector<SXMatrix> dae_in(DAE_NUM_IN);
    dae_in[DAE_T] = t_;
    dae_in[DAE_Y] = var(x_);
    dae_in[DAE_YDOT] = der(x_);
    casadi_assert(0); /*    dae_in[DAE_Z] = var(z_);*/
    dae_in[DAE_P] = var(u_);
    
    // DAE residual
    daeres_ = SXFunction(dae_in,implicit_fcn_);

    // DAE quadrature function
    if(create_quad){
      quadrhs_ = SXFunction(dae_in,lterm[0]/1e3);
    }
  }
  
  // Mayer objective function
  SXMatrix xf = symbolic("xf",x_.size()+lterm.size(),1);
  costfcn_ = SXFunction(xf, xf.data().back());

  // Path constraint function
  vector<SXMatrix> cfcn_in(DAE_NUM_IN);
  cfcn_in[DAE_T] = t_;
  cfcn_in[DAE_Y] = var(x_);
  append(cfcn_in[DAE_Y],symbolic("lterm")); // FIXME
  cfcn_in[DAE_P] = var(u_);
  pathfcn_ = SXFunction(cfcn_in,path_fcn_);
  
  
  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... createFunctions complete after " << dt << " seconds." << endl;
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

void OCP::makeAlgebraic(const Variable& v){
  // Find variable among the variables
  for(vector<Variable>::iterator it=x_.begin(); it!=x_.end(); ++it){
    if(it->get()==v.get()){
      // Add dependent equation, xdot==0
      addExplicitEquation(v.der(),0.0,true);
            
      // Add to list of algebraic variables
      z_.push_back(v);
      
      // Remove from list of differential variables
      x_.erase(it);

      // Replace xdot with x (z) in the list of implicitly defined variables
      for(vector<SX>::iterator jt=implicit_var_.begin(); jt!=implicit_var_.end(); ++jt){
        if(jt->get()==v.der().get()){
          *jt = v.var();
          return;
        }
      }

      // Error if this point reached
      throw CasadiException("variable could not be found in list of implicitly defined variables");
    }
  }
  
  // Error if this point reached
  throw CasadiException("v not in list of differential states");
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

void VariableTree::getAll(std::vector<Variable>& v) const{
  /// Add variables
  for(vector<VariableTree>::const_iterator it=children_.begin(); it!=children_.end(); ++it){
    it->getAll(v);
  }
  
  /// Add variable, if any
  if(!var_.isNull())
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

std::vector<std::string> VariableTree::getNames() const{
  std::vector<std::string> ret(children_.size());
  for(map<string,int>::const_iterator it=name_part_.begin(); it!=name_part_.end(); ++it){
    ret[it->second] = it->first;
  }
  return ret;
}

void OCP::findConsistentIC(){
  // Evaluate the ODE functions
  oderhs_.init();
  oderhs_.setInput(0.0,DAE_T);
  oderhs_.setInput(getStart(x_,true),DAE_Y);
  oderhs_.setInput(getStart(u_,true),DAE_P);
  oderhs_.evaluate();
  
  // Save to the variables
  for(int i=0; i<x_.size(); ++i){
    double xdot0 = oderhs_.output().at(i) * x_[i].getNominal();
    x_[i].setDerivativeStart(xdot0);
  }
  
  // Evaluate the output functions
  output_fcn_.init();
  output_fcn_.setInput(0.0,DAE_T);
  output_fcn_.setInput(getStart(x_,true),DAE_Y);
  output_fcn_.setInput(getStart(u_,true),DAE_P);
  output_fcn_.evaluate();
  
  // Save to the variables
  for(int i=0; i<d_.size(); ++i){
    double z0 = output_fcn_.output().at(i) * d_[i].getNominal();
    d_[i].setStart(z0);
  }
}

SX OCP::getExplicit(const SX& v) const{
  SXMatrix x = v;
  x = substitute(x,explicit_var_,explicit_fcn_);
  return x.toScalar();
}



  } // namespace OptimalControl
} // namespace CasADi

