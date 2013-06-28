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

#include "direct_collocation_internal.hpp"
#include "../symbolic/matrix/matrix_tools.hpp"
#include "../symbolic/sx/sx_tools.hpp"
#include "../symbolic/mx/mx_tools.hpp"
#include "../symbolic/fx/fx_tools.hpp"
#include "../symbolic/stl_vector_tools.hpp"
#include "../symbolic/fx/integrator.hpp"

using namespace std;
namespace CasADi{
    
DirectCollocationInternal::DirectCollocationInternal(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn) : OCPSolverInternal(ffcn, mfcn, cfcn, rfcn){
  addOption("nlp_solver",                       OT_NLPSOLVER,  GenericType(), "An NLPSolver creator function");
  addOption("nlp_solver_options",               OT_DICTIONARY, GenericType(), "Options to be passed to the NLP Solver");
  addOption("interpolation_order",          OT_INTEGER,  3,  "Order of the interpolating polynomials");
  addOption("collocation_scheme",           OT_STRING,  "radau",  "Collocation scheme","radau|legendre");
  casadi_warning("CasADi::DirectCollocation is still experimental");
}

DirectCollocationInternal::~DirectCollocationInternal(){
}

void DirectCollocationInternal::init(){
  // Initialize the base classes
  OCPSolverInternal::init();
  
  // Free parameters currently not supported
  casadi_assert_message(np_==0, "Not implemented");

  // Interpolation order
  deg_ = getOption("interpolation_order");

  // All collocation time points
  std::vector<double> tau_root = collocationPoints(deg_,getOption("collocation_scheme"));

  // Size of the finite elements
  double h = tf_/nk_;

  // Coefficients of the collocation equation
  vector<vector<MX> > C(deg_+1,vector<MX>(deg_+1));

  // Coefficients of the collocation equation as DMatrix
  DMatrix C_num = DMatrix(deg_+1,deg_+1,0);

  // Coefficients of the continuity equation
  vector<MX> D(deg_+1);

  // Coefficients of the collocation equation as DMatrix
  DMatrix D_num = DMatrix(deg_+1,1,0);

  // Collocation point
  SXMatrix tau = ssym("tau");

  // For all collocation points
  for(int j=0; j<deg_+1; ++j){
    // Construct Lagrange polynomials to get the polynomial basis at the collocation point
    SXMatrix L = 1;
    for(int j2=0; j2<deg_+1; ++j2){
      if(j2 != j){
        L *= (tau-tau_root[j2])/(tau_root[j]-tau_root[j2]);
      }
    }

    SXFunction lfcn(tau,L);
    lfcn.init();

    // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    lfcn.setInput(1.0);
    lfcn.evaluate();
    D[j] = lfcn.output();
    D_num(j) = lfcn.output();

    // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    for(int j2=0; j2<deg_+1; ++j2){
      lfcn.setInput(tau_root[j2]);
      lfcn.setFwdSeed(1.0);
      lfcn.evaluate(1,0);
      C[j][j2] = lfcn.fwdSens();
      C_num(j,j2) = lfcn.fwdSens();
    }
  }

  C_num(std::vector<int>(1,0),ALL) = 0;
  C_num(0,0)   = 1;

  // All collocation time points
  vector<vector<double> > T(nk_);
  for(int k=0; k<nk_; ++k){
          T[k].resize(deg_+1);
          for(int j=0; j<=deg_; ++j){
                  T[k][j] = h*(k + tau_root[j]);
          }
  }

  // Total number of variables
  int nlp_nx = 0;
  nlp_nx += nk_*(deg_+1)*nx_;   // Collocated states
  nlp_nx += nk_*nu_;            // Parametrized controls
  nlp_nx += nx_;                       // Final state

  // NLP variable vector
  MX nlp_x = msym("x",nlp_nx);
  int offset = 0;

  // Get collocated states and parametrized control
  vector<vector<MX> > X(nk_+1);
  vector<MX> U(nk_);
  for(int k=0; k<nk_; ++k){
    // Collocated states
        X[k].resize(deg_+1);
    for(int j=0; j<=deg_; ++j){
        // Get the expression for the state vector
        X[k][j] = nlp_x[Slice(offset,offset+nx_)];
        offset += nx_;
    }

    // Parametrized controls
    U[k] = nlp_x[Slice(offset,offset+nu_)];
    offset += nu_;
  }

  // State at end time
  X[nk_].resize(1);
  X[nk_][0] = nlp_x[Slice(offset,offset+nx_)];
  offset += nx_;
  casadi_assert(offset==nlp_nx);

  // Constraint function for the NLP
  vector<MX> nlp_g;

  // Objective function
  MX nlp_j = 0;

  // For all finite elements
  for(int k=0; k<nk_; ++k){

    // For all collocation points
    for(int j=1; j<=deg_; ++j){

        // Get an expression for the state derivative at the collocation point
        MX xp_jk = 0;
        for(int r=0; r<=deg_; ++r){
            xp_jk += C[r][j]*X[k][r];
        }

        // Add collocation equations to the NLP
        MX fk = ffcn_.call(daeIn("x",X[k][j],"p",U[k]))[DAE_ODE];
        nlp_g.push_back(h*fk - xp_jk);
    }

    // Get an expression for the state at the end of the finite element
    MX xf_k = 0;
    for(int r=0; r<=deg_; ++r){
        xf_k += D[r]*X[k][r];
    }

    // Add continuity equation to NLP
    nlp_g.push_back(X[k+1][0] - xf_k);

    // Add path constraints
    if(nh_>0){
      MX pk = cfcn_.call(daeIn("x",X[k+1][0],"p",U[k])).at(0);
      nlp_g.push_back(pk);
    }

    // Add integral objective function term
        //    [Jk] = lfcn.call([X[k+1,0], U[k]])
        //    nlp_j += Jk
  }

  // Add end cost
  MX Jk = mfcn_.call(mayerIn("x",X[nk_][0])).at(0);
  nlp_j += Jk;

  // Objective function of the NLP
  nlp_ = MXFunction(nlpIn("x",nlp_x), nlpOut("f",nlp_j,"g",vertcat(nlp_g)));

  // Get the NLP creator function
  NLPSolverCreator nlp_solver_creator = getOption("nlp_solver");
  
  // Allocate an NLP solver
  nlp_solver_ = nlp_solver_creator(nlp_);
  
  // Pass options
  if(hasSetOption("nlp_solver_options")){
    const Dictionary& nlp_solver_options = getOption("nlp_solver_options");
    nlp_solver_.setOption(nlp_solver_options);
  }
  
  // Initialize the solver
  nlp_solver_.init();
}

void DirectCollocationInternal::getGuess(vector<double>& V_init) const{
  // OCP solution guess
  const Matrix<double> &p_init = input(OCP_P_INIT);
  const Matrix<double> &x_init = input(OCP_X_INIT);
  const Matrix<double> &u_init = input(OCP_U_INIT);
  
  // Running index
  int el=0;
  
  // Pass guess for parameters
  for(int i=0; i<np_; ++i){
    V_init[el++] = p_init.elem(i);
  }
  
  for(int k=0; k<nk_; ++k){
    // Pass guess for state
    for(int j=0; j<=deg_; ++j){
      for(int i=0; i<nx_; ++i){
         V_init[el++] = x_init.elem(i,k);
      }
    }
    
    // Pass guess for control
    for(int i=0; i<nu_; ++i){
      V_init[el++] = u_init.elem(i,k);
    }
  }

  // Pass guess for final state
  for(int i=0; i<nx_; ++i){
    V_init[el++] = x_init.elem(i,nk_);
  }
  
  casadi_assert(el==V_init.size());
}

void DirectCollocationInternal::getVariableBounds(vector<double>& V_min, vector<double>& V_max) const{
  // OCP variable bounds 
  const Matrix<double> &p_min = input(OCP_LBP);
  const Matrix<double> &p_max = input(OCP_UBP);
  const Matrix<double> &x_min = input(OCP_LBX);
  const Matrix<double> &x_max = input(OCP_UBX);
  const Matrix<double> &u_min = input(OCP_LBU);
  const Matrix<double> &u_max = input(OCP_UBU);

  // Running index
  int min_el=0, max_el=0;
  
  // Pass bounds on parameters
  for(int i=0; i<np_; ++i){
    V_min[min_el++] = p_min.elem(i);
    V_max[max_el++] = p_max.elem(i);
  }

  for(int k=0; k<nk_; ++k){
    
    // Pass bounds on state
    for(int i=0; i<nx_; ++i){
      V_min[min_el++] = x_min.elem(i,k);
      V_max[max_el++] = x_max.elem(i,k);
    }
    
    // Pass bounds on collocation points
    for(int j=0; j<deg_; ++j){
      for(int i=0; i<nx_; ++i){
        V_min[min_el++] = std::min(x_min.elem(i,k),x_min.elem(i,k+1));
        V_max[max_el++] = std::max(x_max.elem(i,k),x_max.elem(i,k+1));
      }
    }

    // Pass bounds on control
    for(int i=0; i<nu_; ++i){
      V_min[min_el++] = u_min.elem(i,k);
      V_max[max_el++] = u_max.elem(i,k);
    }
  }

  // Pass bounds on final state
  for(int i=0; i<nx_; ++i){
    V_min[min_el++] = x_min.elem(i,nk_);
    V_max[max_el++] = x_max.elem(i,nk_);
  }
  
  casadi_assert(min_el==V_min.size() && max_el==V_max.size());
}

void DirectCollocationInternal::getConstraintBounds(vector<double>& G_min, vector<double>& G_max) const{
  // OCP constraint bounds
  const Matrix<double> &h_min = input(OCP_LBH);
  const Matrix<double> &h_max = input(OCP_UBH);
  
  // Running index
  int min_el=0, max_el=0;
  
  for(int k=0; k<nk_; ++k){
        for(int j=0; j<=deg_; ++j){
      for(int i=0; i<nx_; ++i){
        G_min[min_el++] = 0.;
        G_max[max_el++] = 0.;
      }
    }
    
    for(int i=0; i<nh_; ++i){
      G_min[min_el++] = h_min.elem(i,k);
      G_max[max_el++] = h_max.elem(i,k);
    }
  }
  casadi_assert(min_el==G_min.size() && max_el==G_max.size());
}

void DirectCollocationInternal::setOptimalSolution( const vector<double> &V_opt ){
  // OCP solution
  Matrix<double> &p_opt = output(OCP_P_OPT);
  Matrix<double> &x_opt = output(OCP_X_OPT);
  Matrix<double> &u_opt = output(OCP_U_OPT);
  
  // Running index
  int el=0;

  // Pass optimized state
  for(int i=0; i<np_; ++i){
    p_opt(i) = V_opt[el++];
  }
    
  for(int k=0; k<nk_; ++k){
    
    // Pass optimized state
    for(int i=0; i<nx_; ++i){
      x_opt(i,k) = V_opt[el++];
    }
    
    // Skip collocation points
    el += deg_*nx_;

    // Pass optimized control
    for(int i=0; i<nu_; ++i){
      u_opt(i,k) = V_opt[el++];
    }
  }

  // Pass optimized terminal state
  for(int i=0; i<nx_; ++i){
    x_opt(i,nk_) = V_opt[el++];
  }
  casadi_assert(el==V_opt.size());
}

void DirectCollocationInternal::evaluate(int nfdir, int nadir){
  // get NLP variable bounds and initial guess
  getGuess(nlp_solver_.input(NLP_SOLVER_X0).data());
  getVariableBounds(nlp_solver_.input(NLP_SOLVER_LBX).data(),nlp_solver_.input(NLP_SOLVER_UBX).data());
       
  // get NLP constraint bounds
  getConstraintBounds(nlp_solver_.input(NLP_SOLVER_LBG).data(), nlp_solver_.input(NLP_SOLVER_UBG).data());
       
  //Solve the problem
  nlp_solver_.solve();
  
  // Save the optimal solution
  setOptimalSolution(nlp_solver_.output(NLP_SOLVER_X).data());

  // Save the optimal cost
  output(OCP_COST).set(nlp_solver_.output(NLP_SOLVER_F));
}


void DirectCollocationInternal::reportConstraints(std::ostream &stream) { 
  stream << "Reporting Collocation constraints" << endl;
 
  CasADi::reportConstraints(stream,output(OCP_X_OPT),input(OCP_LBX),input(OCP_UBX), "states");
  CasADi::reportConstraints(stream,output(OCP_U_OPT),input(OCP_LBU),input(OCP_UBU), "controls");
  CasADi::reportConstraints(stream,output(OCP_P_OPT),input(OCP_LBP),input(OCP_UBP), "parameters");
 
}

} // namespace CasADi
