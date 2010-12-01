/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include <integrator/sundials_interface/cvodes_integrator.hpp>
#include <casadi/stl_vector_tools.hpp>
#include <casadi/expression_tools.hpp>
#include <fstream>
#include <ctime>
#include <iterator>
#include "convection_diffusion_parameters.hpp"
#include "convection_diffusion_model.hpp"

using namespace std;
using namespace OPTICON;


// Generate a two-state model of the storage
SXFunction two_state_model(){
  // Approximate the temperature by a model with 2 parameters
  SX x0("x0");
  SX q("q");

  // parameters
  SX T_l = (T_receiver+273)/(T_return+273); // inlet temperature during loading
  
  Matrix x_meas = linspace(0,1,nz_out); // spatial points
  Matrix T_pred; // predicted T
  for(int j=0; j<nz_out; ++j)
    T_pred << T_l + 0.5*(1-T_l)*(1-erf((x0-x_meas[j])/q));

  // Create the function and initialize it for automatic differentiation
  Matrix par;
  par << x0 << q;
  SXFunction F(par,T_pred);
  F->setOption("ad_order",1);
  F->init();
  return F;
}

// --- Create a function that solves an unconstrained least squares problem of the appropiate size
FX create_LS_solver(){
  // minimize f(x) = 0.5 * ||J*x - y||^2
  // Allocate J and y
  Matrix J("J",nz_out,2);
  Matrix y("y",nz_out);
  
  // Make a symbolic QR factorization of J
  Matrix Q_J, R_J;
  qr(J,Q_J,R_J);

  // Due to orthogonality of Q_J, f(x) is not equal to ||R_J*x - trans(Q_J)*y||^2
  // Solve by a forward substitution
  Matrix res_LSSolver = solve(R_J,trans(Q_J)*y);

  // Create the function
  vector<Matrix> arg_LSSolver(2); // arguments of the least square function
  arg_LSSolver[0] = J;
  arg_LSSolver[1] = y;
  return SXFunction(arg_LSSolver,res_LSSolver);
}

int main(){
  try{
       
  clock_t time0 = clock();
  cout << "program started " << endl;
     
  // create a model instance
  CDModel mod(NZ);

  // Create integrators for the initial value problems
  Integrator integrator_pos = mod.make_integrator_new(1);
  Integrator integrator_neg = mod.make_integrator_new(-1);
  Integrator integrator_zero = mod.make_integrator_new(0);  

  // Generate a simple model
  SXFunction F = two_state_model();

  // Create a solver for the least squares problem
  FX LSSolver = create_LS_solver();

  // Jacobian
  vector<double> jac(x_dim*nz_out);
  
  // Create a file for saving the results
  ofstream resfile;
  resfile.precision(12);
  resfile.open ("results_convection_diffusion.txt");

  vector<double> u_disc(n_u_disc);
  
  // Integrator corresponding to each control
  Integrator integrator[n_u_disc];

  for(int i=0; i<n_u_disc; ++i){
    u_disc[i] = (i-n_u_pos)/double(n_u_pos);
    if(i<n_u_pos){      
      // Negative u
      integrator[i] = integrator_neg;      
    } else if(i==n_u_pos){
      // Zero u
      integrator[n_u_pos] = integrator_zero;
    } else {
      // Positive u
      integrator[i] = integrator_pos;
    }
  }
  
  // Total number of spatial discretizations
  int n_x_disc_total = 1;
  for(int i=0; i<x_dim; ++i)
    n_x_disc_total *= n_x_disc[i];

  // Value for x for each spatial discretization
  vector< vector<double> > x_disc(x_dim);
  for(int i=0; i<x_dim; ++i){
      x_disc[i].resize(n_x_disc[i]);

      for(int j=0; j<n_x_disc[i]; ++j){
	x_disc[i][j] = x_min[i] + (x_max[i]-x_min[i])*j/double(n_x_disc[i]-1);
      }
    }
  {
  
  // Tabulate the new value of the state given a control action
  vector< vector< int > > transition;

  // Tabulate the corresponding value of the objective function
  vector< vector< double > > transition_cost;

  // The least square misfit for the corresponding transition
  vector< vector<double> > transition_misfit;

  // The discretization error (round-off) for the corresponding transition
  vector< vector<double> > transition_roundoff;

  // Current point in the state space
  vector<double> x_cur(x_dim);

  // Discretized 
  vector<int> x_int(x_dim);

  // Current temperature profile
  vector<double> T_prof_cur(nz_out);

  // New temperature profile
  vector<double> T_prof_new(nz_out);

  // deviation from the profile
  vector<double> T_prof_dev(nz_out);

  // for each point in state space
  transition.resize(n_x_disc_total);
  transition_cost.resize(n_x_disc_total);
  transition_misfit.resize(n_x_disc_total);
  transition_roundoff.resize(n_x_disc_total);


   clock_t time1 = clock();
   int progress = 0;
   
   cout << "starting tabulation after " << (time1-time0)/double(CLOCKS_PER_SEC) << " seconds." << endl;
  for(int iss=0; iss<n_x_disc_total; ++iss){
    if((100*iss)/n_x_disc_total > progress){
      progress = (100*iss)/n_x_disc_total;
      cout  << progress << " \%" << endl;
    }
  
    // Get the point in the state space
    int ii=iss;
    for(int j=x_dim-1; j>=0; --j){
      x_cur[j] = x_disc[j][ii%n_x_disc[j]];
      ii /= n_x_disc[j];
    }

    // Print
    if(verbose)	
      cout << "current state space point: " << x_cur << endl;

    // Evaluate the model function to obtain the temperature profile
    F->getInput() = x_cur;
    F->evaluate();
    T_prof_cur = F->getOutput();

    if(verbose)	
      cout << "initial temperature profile: " << T_prof_cur << endl;

    // for each control action
    transition[iss].resize(n_u_disc);
    transition_cost[iss].resize(n_u_disc);
    transition_misfit[iss].resize(n_u_disc);
    transition_roundoff[iss].resize(n_u_disc);

    for(int iu=0; iu<n_u_disc; ++iu){
      if(verbose)
	cout << "current control action: " << u_disc[iu] << endl;

      // Value of the objective function
      double fk = numeric_limits<double>::quiet_NaN();

      // Pass the initial condition to the integrator
      integrator[iu].setInput(T_prof_cur,2);

      // Pass the control action
      integrator[iu].setInput(&u_disc[iu],3);
      
      // Integrate
      integrator[iu].evaluate();

      // Get temperature profile
      integrator[iu].getOutput(T_prof_new);

      // New point in the state space
      vector<double> x_new = x_cur;
      
      // doesn't work well, give some values:
      x_new[0] = 0.5;
      x_new[1] = 0.5;


      if(verbose)	
	cout << "resulting temperature profile: " << T_prof_new << endl;
  
      // Tolerance
      double tol = 1e-6;

      // Find the new point by solving the least squares problem
      for(int iter=0; iter<30; ++iter){ // 20 iterations

	assert(x_new[0] == x_new[0]); // make sure it's not nan
	assert(x_new[1] == x_new[1]); // make sure it's not nan
  
        // Pass the parameter value to the function
        F.setInput(x_new);
  
        // Evaluate the function
        F.evaluate();

        // Get the difference between the fitted profile and the simulated one
        F.getOutput(T_prof_dev);
        for(int j=0; j<nz_out; ++j)
          T_prof_dev[j] -= T_prof_new[j];

        // Calculate the value of the objective function
        fk = 0;
        for(int j=0; j<nz_out; ++j)
          fk += T_prof_dev[j]*T_prof_dev[j];

	// Break if the function is small enough
	if(fk<tol) break;

        // Solve the subproblem: min 0.5*||jac*p + T_cur||^2

        // Right hand side of the LS Solver
        LSSolver.getInput(rhs,1);
        for(int j=0; j<rhs.size(); ++j)
          rhs[j] = -T_prof_dev[j];

	// Get jacobian of F
//	F->evaluateJac(&jac[0]);

	
	
	
	
	// Seed vector
        vector<double> &fseed = F->getInputSeed();
        for(int j=0; j<fseed.size(); ++j)
          fseed[j] = 0;

        // Calculate the jacobian, row by row
        for(int j=0; j<x_dim; ++j){
          // Give a seed in the s-th direction
          fseed[j] = 1;
      
          // Evaluate the directional derivatives
          F->evaluateFwd();

          // Save to jacobian
          for(int l=0; l<nz_out; ++l)
            jac[j + l*x_dim] = F->getOutputSeed()[l];

          // Remove the seed from the seed vector
          fseed[j] = 0;
        }

	// print

	// Jacobian of the LS Solver
        LSSolver.setInput(jac,0);

        // Solve the subproblem
        LSSolver.evaluate();

        // Get the result
        const vector<double> &p = LSSolver->getOutput();

        // Take a step
        for(int j=0; j<x_dim; ++j)
          x_new[j] += 0.5*p[j];

/*        cout << "argval[" << iss << "," << iu << "," << iter << "] = {" << x_new << ", obj = " << fk << endl;*/
	assert(x_new[0] == x_new[0]); // make sure it's not nan
	assert(x_new[1] == x_new[1]); // make sure it's not nan
	assert(fk == fk); // make sure it's not nan
      } // end of for(iter)
    

    if(verbose)	
      cout << "x_new = " << x_new << " obj = " << fk << endl;

    if(verbose)
      cout << "optimal temperature profile: " << F->getOutput() << endl;

    // Save the misfit
    transition_misfit[iss][iu] = fk;

    // Round off to nearest spatial point and save the round-off error
    double err = 0;
    for(int i=0; i<x_dim; ++i){
      x_int[i] = int((x_new[i]-x_min[i])/(x_max[i]-x_min[i])*(n_x_disc[i]-1));

      // Project to box (better: include bounds when solving least squares problem)
      if(x_int[i] < 0)             x_int[i] = 0;
      if(x_int[i] >= n_x_disc[i])  x_int[i] = n_x_disc[i]-1;

     double diff = x_disc[i][x_int[i]] - x_new[i];
     err += diff*diff;
    }

    // Save the round-off error
    transition_roundoff[iss][iu] = err;

    // Save the total cost contribution
    transition_cost[iss][iu] = 0.1*u_disc[iu]*u_disc[iu] + 10*err + 10*fk;

    // Save the transition
    int x_new_int = x_int[0];
    for(int i=1; i<x_dim; ++i){
      x_new_int = x_new_int*n_x_disc[i] + x_int[i];
    }
    transition[iss][iu] = x_new_int;

    // Simulate with the rounded x
    for(int i=0; i<x_dim; ++i){
      x_new[i] = x_disc[i][x_int[i]];
    }

    // Pass the parameter value to the function
    F->getInput() = x_new;
  
    // Evaluate the function
    F->evaluate();

    // Make sure temperature constraints are met
    if(F->getOutput()[0] < 2){
      transition_cost[iss][iu] = numeric_limits<double>::infinity();
    }

    if(verbose){
      cout << "discretized temperature profile: " << F->getOutput() << endl;
    }

    } // for iu
  } // for iss

    if(verbose)	
      cout << "x_disc = " << x_disc << endl;
  
   clock_t time2 = clock();
   cout << "transitions tabulated after " << (time2-time0)/double(CLOCKS_PER_SEC) << " seconds." << endl;

  resfile.close();
  resfile.open ("transition.txt");
  write_matlab(resfile,transition);
  resfile.close();
  resfile.open ("transition_misfit.txt");
  write_matlab(resfile,transition_misfit);
  resfile.close();
  resfile.open ("transition_roundoff.txt");
  write_matlab(resfile,transition_roundoff);
  resfile.close();
  resfile.open ("transition_cost.txt");
  write_matlab(resfile,transition_cost);
  resfile.close();
  resfile.open ("x_int.txt");
  write_matlab(resfile,x_int);
  resfile.close();

  
  cout << "saved to disk" << endl;
}

  return 0;
} catch (exception &e){
  cerr << e.what() << endl;  
  return 1;
} catch (const char * str){
  cerr << str << endl;
  return 1;
}

}
