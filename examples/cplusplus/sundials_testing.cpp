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

#include <interfaces/sundials/cvodes_integrator.hpp>
#include <interfaces/sundials/idas_integrator.hpp>
#include <interfaces/lapack/lapack_lu_dense.hpp>
#include <interfaces/csparse/csparse.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <symbolic/fx/simulator.hpp>
#include <symbolic/fx/c_function.hpp>
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/sx_function.hpp"
#include "symbolic/fx/jacobian.hpp"
#include <symbolic/casadi_calculus.hpp>

#include <fstream>
#include <iostream>

using namespace std;
using namespace CasADi;

// Use CVodes or IDAS
const bool implicit_integrator = false;

// use plain c instead of SX
const bool plain_c = false;

// test adjoint sensitivities
const bool with_asens = true;

// use exact jacobian
const bool exact_jacobian = plain_c ? false : true;

// Calculate the forward sensitivities using finite differences
const bool finite_difference_fsens = !exact_jacobian;

// Calculate initial condition (for IDAS only)
const bool calc_ic = true;

// Perturb x or u
const bool perturb_u = true;

// Use a user_defined linear solver
const bool user_defined_solver = false;

// Use sparse direct solver (CSparse)
const bool sparse_direct = true;

// Second order sensitivities by a symbolic-numeric approach
const bool second_order = true;

// The DAE residual in plain c (for IDAS)
void dae_res_c(double t, const double *x, const double* xdot, const double* p, double* ode, double* quad){
  // Get the arguments
  //double s = x[0];
  double v = x[1];
  double m = x[2];
  double u = p[0];
  double sdot = xdot[0], vdot = xdot[1], mdot = xdot[2];

  // Calculate the DAE residual
  ode[0] = sdot - v;
  ode[1] = vdot - (u-0.02*v*v)/m;
  ode[2] = mdot - (-0.01*u*u);
  
  double u_ref = 3-sin(t);
  double u_dev = u-u_dev;
  quad[0] = u_dev*u_dev;
}

// Wrap the function to allow creating an CasADi function
void dae_res_c_wrapper(CFunction &f, int nfwd, int nadj, void* user_data){
  casadi_assert(nfwd==0 && nadj==0);
  dae_res_c(f.input(DAE_T).front(), &f.input(DAE_X).front(), &f.input(DAE_XDOT).front(), &f.input(DAE_P).front(), &f.output(DAE_ODE).front(), &f.output(DAE_QUAD).front());
}

// Create an IDAS instance (fully implicit integrator)
Integrator create_Sundials(){
  // Time 
  SX t("t");

  // Differential states
  SX s("s"), v("v"), m("m");
  vector<SX> x(3); 
  x[0] = s;
  x[1] = v;
  x[2] = m;
  
  // State derivatives
  SX sdot("sdot"), vdot("vdot"), mdot("mdot");
  vector<SX> xdot(3); 
  xdot[0] = sdot;
  xdot[1] = vdot;
  xdot[2] = mdot;

  // Control
  SX u("u");
  
  // Reference trajectory
  SX u_ref = 3-sin(t);
  
  // Square deviation from the state trajectory
  SX u_dev = u-u_ref;
  u_dev *= u_dev;
  
  // Differential equation (fully implicit form)
  vector<SX> res(3);
  res[0] = v - sdot;
  res[1] = (u-0.02*v*v)/m - vdot;
  res[2] = -0.01*u*u - mdot;

  // Input/output of the DAE residual function
  vector<SXMatrix> ffcn_in = daeIn<SXMatrix>("x",x, "p",u, "t",t, "xdot",xdot);
  vector<SXMatrix> ffcn_out = daeOut<SXMatrix>("ode",res, "quad",u_dev);

  // DAE residual function
  FX ffcn = SXFunction(ffcn_in,ffcn_out);

  // Overwrite ffcn with a plain c function (avoid this!)
  if(plain_c){
    // Use DAE residual defined in a c-function
    ffcn = CFunction(dae_res_c_wrapper);
    
    // Specify the number of inputs and outputs
    ffcn.setNumInputs(DAE_NUM_IN);
    ffcn.setNumOutputs(DAE_NUM_OUT);
    
    // Specify dimensions of inputs and outputs
    ffcn.input(DAE_T)    = DMatrix(1,1,0);
    ffcn.input(DAE_X)    = DMatrix(3,1,0);
    ffcn.input(DAE_XDOT) = DMatrix(3,1,0);
    ffcn.input(DAE_P)    = DMatrix(1,1,0);
    ffcn.output(DAE_ODE) = DMatrix(3,1,0);
    ffcn.output(DAE_QUAD) = DMatrix(1,1,0);
  }
  
  if(implicit_integrator){
    // Create an IDAS instance
    IdasIntegrator integrator(ffcn);
    
    // Set IDAS specific options
    integrator.setOption("calc_ic",calc_ic);

    // Return the integrator
    return integrator;
  } else {
    // Create an CVodes instance
    CVodesIntegrator integrator(ffcn);

    // Return the integrator
    return integrator;
  }
}

int main(){
  // Time horizon
  double t0 = 0,  tf = 10;
  
  // Bounds on the control
  double /*u_lb = -0.5, u_ub = 1.3,*/ u_init = 1;

  // Initial conditions
  vector<double> x0(3);
  x0[0] = 0;
  x0[1] = 0;
  x0[2] = 1;
  
  // Integrator
  Integrator integrator = create_Sundials();
  
  // Attach user-defined linear solver
  if(user_defined_solver){
    if(sparse_direct){
      integrator.setOption("linear_solver_creator",CSparse::creator);
      // integrator.setOption("linear_solver_creator",SuperLU::creator);
    } else {
      integrator.setOption("linear_solver_creator",LapackLUDense::creator);
      integrator.setOption("linear_solver","user_defined");
    }
    // integrator.setOption("linear_solver","user_defined"); // FIXME: bug for second order
  }
  
  // Set common integrator options
  integrator.setOption("fsens_err_con",true);
  integrator.setOption("quad_err_con",true);
  integrator.setOption("abstol",1e-12);
  integrator.setOption("reltol",1e-12);
  integrator.setOption("fsens_abstol",1e-6);
  integrator.setOption("fsens_reltol",1e-6);
  integrator.setOption("asens_abstol",1e-6);
  integrator.setOption("asens_reltol",1e-6);
  integrator.setOption("exact_jacobian",exact_jacobian);
  integrator.setOption("finite_difference_fsens",finite_difference_fsens);
  integrator.setOption("max_num_steps",100000);
//  integrator.setOption("max_multistep_order",4);
  integrator.setOption("t0",t0);
  integrator.setOption("tf",tf);

  // Initialize the integrator
  integrator.init();
 
  // Set parameters
  integrator.setInput(u_init,INTEGRATOR_P);
  
  // Set inital state
  integrator.setInput(x0,INTEGRATOR_X0);
  
  // Integrate
  integrator.evaluate();

  // Save the result
  Matrix<double> res0_xf = integrator.output(INTEGRATOR_XF);
  Matrix<double> res0_qf = integrator.output(INTEGRATOR_QF);

  // Perturb in some direction
  if(perturb_u){
    double u_pert = u_init + 0.01;
    integrator.setInput(u_pert,INTEGRATOR_P);
  } else {
    vector<double> x_pert = x0;
    x_pert[1] += 0.01;
    integrator.setInput(x_pert,INTEGRATOR_X0);
  }
  
  // Integrate again
  integrator.evaluate();
  
  // Print statistics
  integrator.printStats();

  // Calculate finite difference approximation
  Matrix<double> fd_xf = (integrator.output(INTEGRATOR_XF) - res0_xf)/0.01;
  Matrix<double> fd_qf = (integrator.output(INTEGRATOR_QF) - res0_qf)/0.01;

  cout << "unperturbed                     " << res0_xf << "; " << res0_qf << endl;
  cout << "perturbed                       " << integrator.output(INTEGRATOR_XF) << "; " << integrator.output(INTEGRATOR_QF) << endl;
  cout << "finite_difference approximation " << fd_xf << "; " << fd_qf << endl;

  if(perturb_u){
    integrator.setFwdSeed(1.0,INTEGRATOR_P);
  } else {
    vector<double> x0_seed(x0.size(),0);
    x0_seed[1] = 1;
    integrator.setFwdSeed(x0_seed,INTEGRATOR_X0);
  }
    
  // Reset parameters
  integrator.setInput(u_init,INTEGRATOR_P);
  
  // Reset initial state
  integrator.setInput(x0,INTEGRATOR_X0);

  if(with_asens){
    // backward seeds
    vector<double> &bseed = integrator.adjSeed(INTEGRATOR_XF).data();
    fill(bseed.begin(),bseed.end(),0);
    bseed[1] = 1;

    // evaluate with forward and adjoint sensitivities
    integrator.evaluate(1,1);
  } else {
    // evaluate with only forward sensitivities
    integrator.evaluate(1,0);
  }
    
  Matrix<double> fsens_xf = integrator.fwdSens(INTEGRATOR_XF);
  Matrix<double> fsens_qf = integrator.fwdSens(INTEGRATOR_QF);
  cout << "forward sensitivities           " << fsens_xf << "; " << fsens_qf << endl;

  if(with_asens){
    cout << "adjoint sensitivities           ";
    cout << integrator.adjSens(INTEGRATOR_X0) << "; ";
    cout << integrator.adjSens(INTEGRATOR_P) << "; ";
    cout << endl;
  }
  
  if(second_order){
    // Preturb the forward seeds
    if(perturb_u){
      double u_seed = 1.001;
      integrator.setFwdSeed(u_seed,INTEGRATOR_P);
    } else {
      vector<double> x0_seed(x0.size(),0);
      x0_seed[1] = 1.001;
      integrator.setFwdSeed(x0_seed,INTEGRATOR_X0);
    }
    
    // evaluate again with forward sensitivities
    integrator.evaluate(1,0);

    vector<double> fsens_pret_xf = integrator.fwdSens(INTEGRATOR_XF).data();
    vector<double> fsens_pret_qf = integrator.fwdSens(INTEGRATOR_QF).data();
    cout << "forward sensitivities preturbed " << fsens_pret_xf << "; " << fsens_pret_qf << endl;

    vector<double> fd2_xf(fsens_xf.size());
    vector<double> fd2_qf(fsens_qf.size());
    for(int i=0; i<fd2_xf.size(); ++i)
      fd2_xf[i] = (fsens_pret_xf.at(i)-fsens_xf.at(i))/0.001;
    for(int i=0; i<fd2_qf.size(); ++i)
      fd2_qf[i] = (fsens_pret_qf.at(i)-fsens_qf.at(i))/0.001;
    
    cout << "finite differences, 2nd order   " << fd2_xf << "; " << fd2_qf << endl;
    
    // Generate the jacobian by creating a new integrator for the sensitivity equations by source transformation
    FX intjac  = integrator.jacobian(INTEGRATOR_P,INTEGRATOR_XF);
    FX intjac2 = Jacobian(integrator,INTEGRATOR_P,INTEGRATOR_XF);

    // Set options
    intjac.setOption("number_of_fwd_dir",0);
    intjac.setOption("number_of_adj_dir",1);
    
    // Initialize the integrator
    intjac.init();
    intjac2.init();

    // Set inputs
    intjac.setInput(u_init,INTEGRATOR_P);
    intjac.setInput(x0,INTEGRATOR_X0);
    
    intjac2.setInput(u_init,INTEGRATOR_P);
    intjac2.setInput(x0,INTEGRATOR_X0);
    
    // Set adjoint seed
    vector<double> jacseed(3*1,0);
    jacseed[0] = 1;
    intjac.setAdjSeed(jacseed);
    
    // Evaluate the Jacobian
    intjac.evaluate(0,0);
    intjac2.evaluate(0,0);

    
    cout << "jacobian                  " << intjac.output(0) << endl;
    cout << "jacobian2                 " << intjac2.output(0) << endl;
    
    // Get the results
/*    cout << "unperturbed via jacobian        " << intjac.output(1+INTEGRATOR_XF) << endl;*/
    cout << "second order (fwd-over-adj)     " ;
    cout << intjac.adjSens(INTEGRATOR_X0) << ", ";
    cout << intjac.adjSens(INTEGRATOR_P) << endl;

    // Save the unpreturbed value
    Matrix<double> unpret = intjac.output();
    Matrix<double> unpret2 = intjac2.output();
    
    // Perturb X0
    intjac.setInput(u_init+0.01,INTEGRATOR_P);
    intjac2.setInput(u_init+0.01,INTEGRATOR_P);

    intjac.evaluate();
    intjac2.evaluate();
    Matrix<double> pret = intjac.output();
    Matrix<double> pret2 = intjac2.output();
    
    cout << "unperturbed fwd sens            " << unpret << endl;
    cout << "perturbed fwd sens              " << pret << endl;
    cout << "finite diff. (augmented dae)    " << (pret-unpret)/0.01 << endl;
    
    cout << "unperturbed fwd sens 2           " << unpret2 << endl;
    cout << "perturbed fwd sens 2             " << pret2 << endl;
    cout << "finite diff. (augmented dae) 2   " << (pret2-unpret2)/0.01 << endl;
    
  }
  
  return 0;
}

