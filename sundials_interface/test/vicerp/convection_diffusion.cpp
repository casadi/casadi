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

#include <fstream>
#include <ctime>
#include "convection_diffusion_model.hpp"

// Initial value problems corresponding to positive, negative and zero flow
// note: global variables since integrators contains a reference to an ocp, not a smart pointer, this will change
  OCP ivp_pos, ivp_neg, ivp_zero;

  // Control (flow to storage)
  SX u("u");

  // Maximal u
  SX u_max = 0.33; // fully loading/unloading takes 1.5 hours

  // Initial temperature profile
  Matrix T0("T0",NZ);

  // infinity
  double infty = numeric_limits<double>::infinity();

  // Number of times
  int nk_per_hour = 1;
  int NK = 24*nk_per_hour;

  // Time for a control
  SX DT = 1.0/nk_per_hour;

  // Discretize space
  const int n_dim = 2; // number of dimensions
  
  // Discretization for each dimension
  int n_x_disc[] = {20, 20};

  // Lower and upper bound on each variable
  double x_min[] = {-1, 0.01};
  double x_max[] = {1,  1};

  // Discretize the control
  const int n_u_pos = 10; // in the positive direction direction
//  const int n_u_pos = 4; // in the positive direction direction
  const int n_u_disc = 2*n_u_pos+1; // n_u_pos downwind, zero and n_u_pos upwind


// Assemble the initial value problems for one shooting interval
void calc_ivp(){

  // Common for all
  for(int i=0; i<3; ++i){
    // reference to the ivp
    OCP &ivp = i==0 ? ivp_pos : i==1 ? ivp_neg : ivp_zero;

    // Time variable
    Matrix t = ivp.getTime();
    
    // Specify the time horizon
    ivp.setInitialTime(0);
    ivp.setFinalTime(DT);
    
    // differential states
    ivp.addState(T,Tdot);

    // control
    ivp.addParameter(u);

    // Initial temperature in the storage
    ivp.addParameter(T0);
    ivp.addInitialCondition(T, T0);

    // Set output function (final T)
    ivp.addOutput(DT, T);
  }

  // Differential equation (positive u)
  ivp_pos.addEquation(Tdot, -u*u_max*dT_dz_pos + D*d2T_dz2);

  // Differential equation (negative u)
  ivp_neg.addEquation(Tdot, -u*u_max*dT_dz_neg + D*d2T_dz2);

  // Differential equation (zero u)
  ivp_zero.addEquation(Tdot, D*d2T_dz2);

}

// Generate a two-state model of the storage
SXFunction two_state_model(){
  // Approximate the temperature by a model with 2 parameters
  SX x0("x0");
  SX q("q");
  Matrix x_meas = linspace(0,1,nz_out); // spatial points
  Matrix T_pred; // predicted T
  for(int j=0; j<nz_out; ++j)
    T_pred << T_l + 0.5*(1-T_l)*(1-erf((x0-x_meas[j])/q));

  // Create the function and initialize it for automatic differentiation
  Matrix par;
  par << x0 << q;
  SXFunction F(par,T_pred);
  F->initAD(1);
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

  // Calculate the derivative approximations
  calc_Tdisc();

  // Create a file for saving the results
  ofstream resfile;
  resfile.open ("results_convection_diffusion.txt");

  // Vector for storing the temperature profile
  vector<double> t_prof;

  // Simulate the model
  simulate(resfile, t_prof);

  // Assemble the initial value problems
  calc_ivp();

  // Allocate integrators for the initial value problems
  CvodesIntegrator integrator_pos(ivp_pos);
  CvodesIntegrator integrator_neg(ivp_neg);
  CvodesIntegrator integrator_zero(ivp_zero);

  // Initialize
  for(int i=0; i<3; ++i){
    CvodesIntegrator& integrator = i==0 ? integrator_pos : i==1 ? integrator_neg : integrator_zero;

    // Set the linear solver
    integrator->setOption("linear_solver","band");

    // Upper and lower band-widths (only relevant when using a band linear solver)
    integrator->setOption("mupper",1);
    integrator->setOption("mlower",1);

    // Use exact jacobian
    integrator->setOption("exact_jacobian",false);

    // set tolerances 
    integrator->setOption("reltol",1e-6);
    integrator->setOption("abstol",1e-8);

    // Initialize the integrator
    integrator->init();
    cout << "initialized" << endl;

    // Set options
    integrator->setOption("use_events",  false);
    integrator->setOption("use_output",  true);
    integrator->setOption("use_fsens",   false);
    integrator->setOption("use_asens",   false);
    integrator->setOption("use_quad",    false);
  }

  // Generate a simple model
  SXFunction F = two_state_model();

  // Create a solver for the least squares problem
  FX LSSolver = create_LS_solver();

  vector<double> u_disc(n_u_disc);
  
  // Integrator corresponding to each control
  Integrator integrator[n_u_disc];

  // Initial value problem for each control
  OCP* ivp[n_u_disc];

  // The indices of the control
  vector<int> u_ind[n_u_disc];

  // The indices of the initial conditions
  vector<int> ic_ind[n_u_disc];

  // Negative u
  for(int i=0; i<n_u_pos; ++i){
    integrator[i] = integrator_neg;
    ivp[i] = &ivp_neg;
    u_disc[i] = i/double(n_u_pos) - 1;
    u_ind[i] = ivp_neg.getVarInds(u);
    ic_ind[i] = ivp_neg.getVarInds(T0);
  }

  // Zero u
  integrator[n_u_pos] = integrator_zero;
  ivp[n_u_pos] = &ivp_zero;
  u_disc[n_u_pos] = 0;
  u_ind[n_u_pos] = ivp_zero.getVarInds(u);
  ic_ind[n_u_pos] = ivp_zero.getVarInds(T0);

  // Positive u
  for(int i=0; i<n_u_pos; ++i){
    integrator[n_u_pos+1+i] = integrator_pos;
    ivp[n_u_pos+1+i] = &ivp_pos;
    u_disc[n_u_pos+1+i] = (i+1)/double(n_u_pos);
    u_ind[n_u_pos+1+i] = ivp_pos.getVarInds(u);
    ic_ind[n_u_pos+1+i] = ivp_pos.getVarInds(T0);
  }

  // Test the integrator
  vector<double> temp(NZ);
  for(int i=0; i<NZ; ++i)
    temp[i] = 1 + ((680+273)/double(120+273)-1) * (i<NZ/2);

  // Loop over all controls
  for(int i=0; i<n_u_disc; ++i){
    // Pass the temperature profile
    ivp[i]->setArg(ic_ind[i], &temp[0]);

    // Pass a control
    ivp[i]->setArg(u_ind[i], &u_disc[i]);

    // Evaluate to get the output
    integrator[i]->evaluate();
/*    cout << "res" << i << " " << integrator[i]->getOutput() << endl;*/
    resfile << "res" << i << " " << integrator[i]->getOutput() << endl;
    
  }


  // Total number of spatial discretizations
  int n_x_disc_total = 1;
  for(int i=0; i<n_dim; ++i)
    n_x_disc_total *= n_x_disc[i];

  // Value for x for each spatial discretization
  vector< vector<double> > x_disc(n_dim);
  for(int i=0; i<n_dim; ++i){
    x_disc[i].resize(n_x_disc[i]);

    for(int j=0; j<n_x_disc[i]; ++j){
      x_disc[i][j] = x_min[i] + (x_max[i]-x_min[i])*j/double(n_x_disc[i]-1);
    }
  }

  // Tabulate the new value of the state given a control action
  vector< vector< int > > transition;

  // Tabulate the corresponding value of the objective function
  vector< vector< double > > transition_cost;

  // The least square misfit for the corresponding transition
  vector< vector<double> > transition_misfit;

  // The discretization error (round-off) for the corresponding transition
  vector< vector<double> > transition_roundoff;

  // Current point in the state space
  vector<double> x_cur(n_dim);

  // Discretized 
  vector<int> x_int(n_dim);

  // Current temperature profile
  vector<double> T_prof_cur(nz_out);

  // deviation from the profile
  vector<double> T_prof_dev(nz_out);

  // for each point in state space
  transition.resize(n_x_disc_total);
  transition_cost.resize(n_x_disc_total);
  transition_misfit.resize(n_x_disc_total);
  transition_roundoff.resize(n_x_disc_total);


   clock_t time1 = clock();
   cout << "starting tabulation after " << (time1-time0)/double(CLOCKS_PER_SEC) << " seconds." << endl;
  for(int iss=0; iss<n_x_disc_total; ++iss){
  
    // Get the point in the state space
    int ii=iss;
    for(int j=n_dim-1; j>=0; --j){
      x_cur[j] = x_disc[j][ii%n_x_disc[j]];
      ii /= n_x_disc[j];
    }

    // Print
    if(print)	
      cout << "current state space point: " << x_cur << endl;

    // Evaluate the model function to obtain the temperature profile
    F->getInput() = x_cur;
    F->evaluate();
    T_prof_cur = F->getOutput();

    if(print)	
      cout << "initial temperature profile: " << T_prof_cur << endl;

    // for each control action
    transition[iss].resize(n_u_disc);
    transition_cost[iss].resize(n_u_disc);
    transition_misfit[iss].resize(n_u_disc);
    transition_roundoff[iss].resize(n_u_disc);

    for(int iu=0; iu<n_u_disc; ++iu){
      if(print)	
	cout << "current control action: " << u_disc[iu] << endl;

      // Value of the objective function
      double fk = numeric_limits<double>::quiet_NaN();

//       // Allocate memory for the result
//       transition[iss][iu].resize(n_dim);

      // Pass the initial condition to the integrator
      ivp[iu]->setArg(ic_ind[iu], &T_prof_cur[0]);

      // Pass the control action
      ivp[iu]->setArg(u_ind[iu], &u_disc[iu]);
      
      // Integrate
      try{
        integrator[iu]->evaluate();
      } catch(...){
         // The integration failed: set the value to be not a number
          for(int l=0; l<n_dim; ++l){
	    assert(0);
//            transition[iss][iu][l] = numeric_limits<double>::quiet_NaN();
	  }
      }

      // New point in the state space
      vector<double> x_new = x_cur;
      
      // doesn't work well, give some values:
      x_new[0] = 0.5;
      x_new[1] = 0.5;

      // Temperature profile
      const vector<double> &T_prof_new = integrator[iu]->getOutput();

      if(print)	
	cout << "resulting temperature profile: " << T_prof_new << endl;
  
      // Tolerance
      double tol = 1e-6;

      // Find the new point by solving the least squares problem
      for(int iter=0; iter<30; ++iter){ // 20 iterations

	assert(x_new[0] == x_new[0]); // make sure it's not nan
	assert(x_new[1] == x_new[1]); // make sure it's not nan

  
        // Pass the parameter value to the function
        F->getInput() = x_new;
  
        // Evaluate the function
        F->evaluate();

        // Get the difference between the fitted profile and the simulated one
        T_prof_dev = F->getOutput();
        for(int j=0; j<nz_out; ++j)
          T_prof_dev[j] -= T_prof_new[j];

//  	cout << "T_prof_dev = " << T_prof_dev << endl;

        // Calculate the value of the objective function
        fk = 0;
        for(int j=0; j<nz_out; ++j)
          fk += T_prof_dev[j]*T_prof_dev[j];

	// Break if the function is small enough
	if(fk<tol) break;

        // Solve the subproblem: min 0.5*||jac*p + T_cur||^2

        // Right hand side of the LS Solver
        vector<double> &rhs = LSSolver->getInput(1);
        for(int j=0; j<rhs.size(); ++j)
          rhs[j] = -T_prof_dev[j];

        // Seed vector
        vector<double> &fseed = F->getInputSeed();
        for(int j=0; j<fseed.size(); ++j)
          fseed[j] = 0;

        // Jacobian of the LS Solver
        vector<double> &jac = LSSolver->getInput(0);


// cout << endl<< "F->work[0] = " << F->work[0] << endl<< endl;
// cout << endl<< "F->work[1] = " << F->work[1] << endl<< endl;

        // Calculate the jacobian, row by row
        for(int j=0; j<n_dim; ++j){
          // Give a seed in the s-th direction
          fseed[j] = 1;
      
          // Evaluate the directional derivatives
          F->evaluateFwd();

          // Save to jacobian
          for(int l=0; l<nz_out; ++l)
            jac[j + l*n_dim] = F->getOutputSeed()[l];

          // Remove the seed from the seed vector
          fseed[j] = 0;
        }

	// print
//  	cout << "jacobian: " << jac << endl;

        // Solve the subproblem
        LSSolver->evaluate();

        // Get the result
        const vector<double> &p = LSSolver->getOutput();

        // Take a step
        for(int j=0; j<n_dim; ++j)
          x_new[j] += 0.5*p[j];



/*        cout << "argval[" << iss << "," << iu << "," << iter << "] = {" << x_new << ", obj = " << fk << endl;*/
	assert(x_new[0] == x_new[0]); // make sure it's not nan
	assert(x_new[1] == x_new[1]); // make sure it's not nan
	assert(fk == fk); // make sure it's not nan
      } // end of for(iter)
    

    if(print)	
      cout << "x_new = " << x_new << " obj = " << fk << endl;

    if(print)
      cout << "optimal temperature profile: " << F->getOutput() << endl;

    // Save the misfit
    transition_misfit[iss][iu] = fk;

    // Round off to nearest spatial point and save the round-off error
    double err = 0;
    for(int i=0; i<n_dim; ++i){
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
    for(int i=1; i<n_dim; ++i){
      x_new_int = x_new_int*n_x_disc[i] + x_int[i];
    }
    transition[iss][iu] = x_new_int;

    // Simulate with the rounded x
    for(int i=0; i<n_dim; ++i){
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

    if(print){
      cout << "discretized temperature profile: " << F->getOutput() << endl;
    }

    } // for iu
  } // for iss

    if(print)	
      cout << "x_disc = " << x_disc << endl;
  
   clock_t time2 = clock();
   cout << "transitions tabulated after " << (time2-time0)/double(CLOCKS_PER_SEC) << " seconds." << endl;

  resfile << "transition " << endl << transition << endl;
  resfile << "transition_misfit " << endl << transition_misfit << endl;
  resfile << "transition_roundoff " << endl << transition_roundoff << endl;
  resfile << "transition_cost " << endl << transition_cost << endl;

  cout << "starting dynamic programming algorithm" << endl;
  
  // Cost to arrive function
  vector< vector<double> > J(NK+1);

  // Give a starting point (half-full storage)
  J[0].resize(n_x_disc_total,infty);
  x_int[0] = n_x_disc[0]/2;
  x_int[1] = 0;
  int init_ind = x_int[0];
  for(int i=1; i<n_dim; ++i){
    init_ind = init_ind*n_x_disc[i] + x_int[i];
  }
  J[0][init_ind] = 0;
  
  // Optimal control trajectory table
  vector< vector<int> > u_table(NK);

  // state trajectory table
  vector< vector<int> > x_table(NK);

  // DP algorithm
  for(int k=0; k<NK; ++k){

    // Allocate space for the optimal control and cost to arrive
    u_table[k].resize(n_x_disc_total);
    x_table[k].resize(n_x_disc_total);
    J[k+1].resize(n_x_disc_total,infty);

    // Loop over the points in state space
    for(int i=0; i<n_x_disc_total; ++i){

      if(J[k][i]==0)
	cout << "transition[" << i << "] = " << transition[i] << endl;

      // Loop over all possible control actions
      for(int j=0; j<n_u_disc; ++j){

	// Evaluate the dynamic equation
	int i_new = transition[i][j];
      
	// Flow to heat exchanger
	double flow = 0;
	if(k<12*nk_per_hour)
	  flow = (k*(12*nk_per_hour-k))/double(12*nk_per_hour*12*nk_per_hour/4);
	flow -= u_disc[j];

	// Heat exchanger flow setpoint
	double flow_ref = 0;
	if(k>=nk_per_hour && k<nk_per_hour*11){
	  // plant turned on
	  flow_ref = 0.8;
	}

	// Cost at time k
	double cost_k = flow-flow_ref;
	cost_k *= cost_k;

	// Calculate the optimal cost
	double cost = J[k][i] + transition_cost[i][j] + cost_k;

	// Infinite cost if flow to heat exchanger is negative
	if(flow<0)
	  cost = cost + 100;
//	  cost = infty;

	// Save, if better
	if(J[k+1][i_new] > cost){
	  J[k+1][i_new] = cost;
	  u_table[k][i_new] = j;
	  x_table[k][i_new] = i;
	}
     } // for j
   } // for i
 } // for k 

  // Find the optimal cost
  int i_min = -1;
  double c_min = infty;
  for(int i=0; i<n_x_disc_total; ++i)
    if(J[NK][i] < c_min){
      c_min = J[NK][i];
      i_min = i;
    }

  cout << "i_min = " << i_min << ", J[NK][i_min] = " << J[NK][i_min] << endl;
  
  // Recover the optimal policy
  cout << "rollback" << endl;

  vector<double> u_opt(NK);  
  for(int k=NK-1; k>=0; --k){
    u_opt[k] = u_disc[u_table[k][i_min]];
    i_min = x_table[k][i_min];
  }

  // Print the optimal control
  cout << "optimal control: " << u_opt << endl;



  // Print the cost to arrive
//  cout << "J[NK] = " << J[NK] << endl;
  

   clock_t time3 = clock();
   cout << "program finished after " << (time3-time0)/double(CLOCKS_PER_SEC) << " seconds." << endl;



return 0;
//  vector<int>



#if 0

  // give the result back to the initial condition
  ivp_neg.setArg(ic_neg, &integrator_neg->getOutput()[0]);

  integrator_neg->evaluate();
  cout << "res2 " << integrator_neg->getOutput() << endl;

#endif



  // Fit two parameter curve to the data
  vector<double> x0_opt(nt_out);
  vector<double> q_opt(nt_out);
  vector<double> res_opt(nz_out*nt_out);

  // Residual vector
  double rk[nz_out];

  // Jacobian of residual vector
  double Jk[nz_out][2];

  // temporaries
  double argval[2];
  double temppar[2];
  double tempres[nz_out];

  // Fit a curve to the temperature profile
 for(int i=0; i<nt_out; ++i){
    // Make a guess for the parameter values
    argval[0] = 0.0; // x0
    argval[1] = 0.1; // q
//    cout << "argval0[" << i << "] = {" << argval[0] << "," << argval[1] << "}" << endl;

    // Temperature profile
    double *T_meas2 = &t_prof[i*nz_out];

    for(int k=0; k<20; ++k){
  
    // Pass the parameter value to the function
    F->getInput()[0] = argval[0];
    F->getInput()[1] = argval[1];
  
    // Evaluate the function
    F->evaluate();

    // Get the result
    F->getOutputDense(&rk[0]);

    // Subtact the "real" value
    for(int j=0; j<nz_out; ++j)
      rk[j] -= T_meas2[j];

    // Calculate the value of the objective function
    double fk = 0;
    for(int j=0; j<nz_out; ++j)
      fk += rk[j]*rk[j];

    // Seed vector
    double fseed[2] = {0}; // all zeros

    // Calculate the jacobian, row by row
    for(int j=0; j<2; ++j){
      // Clear all derivative values
      F->clear(1);

      // Give a seed in the j-th direction
      fseed[j] = 1;
      for(int r=0; r<2; ++r)
        F->getInputSeed()[r] = fseed[r];
      
      // Evaluate the directional derivatives
      F->evaluateFwd();

      // Get the result
      F->getOutputDense(&tempres[0],0,1);
  
      // Save to jacobian
      for(int l=0; l<nz_out; ++l)
         Jk[l][j] = tempres[l];

      // Remove the seed from the seed vector
      fseed[j] = 0;
    }

    // Solve the subproblem: min 0.5*||J_k*p + rk||^2
    // Pass the jacobian to the LSSolver
    vector<double> &ip = LSSolver->getInput();
    for(int r=0; r<ip.size(); ++r)
      ip[r] = Jk[0][r];

    // Pass the right hand side to the LSSolver
    for(int l=0; l<nz_out; ++l)
       LSSolver->getInput(1)[l] = -rk[l];
    
    // Solve the subproblem
    LSSolver->evaluate();

    // Get the result
    LSSolver->getOutputDense(&temppar[0]);

    // Take a step
    for(int j=0; j<2; ++j)
      argval[j] += 0.5*temppar[j];


//    cout << "argval[" << i << "," << k << "] = {" << argval[0] << "," << argval[1] << "}, obj = " << fk << endl;
}

  // Save the parameters

  argval[1] = i/double(nt_out-1);
  argval[0] = -0.1;

  x0_opt[i] = argval[0];
  q_opt[i] = argval[1];
 
  // Calculate the final residual
  vector<double> &ip = F->getInput();
  for(int r=0; r<ip.size(); ++r)
    ip[r] = argval[r];

  F->evaluate();
  F->getOutputDense(&res_opt[i*nz_out]);
}

  resfile << "x0_opt " << x0_opt << endl;
  resfile << "q_opt " << q_opt << endl;
  resfile << "res_opt " << res_opt << endl;

  // Close the results file
  resfile.close();




# if 0

    for(int j=0; j<nz_out; ++j)
      cout << Ti[j] << ", ";
  
    cout << endl << endl;



    // Create a function
    Matrix arg; arg << x0 << q;

    // Create the NLP
    NLP nlp(R,arg);

    // Bounds on u and initial condition
    vector<double> umin(2), umax(2), usol(2);
    umin[0] = 0;
    umax[0] = 1;
    usol[0] = 0.2;

    umin[1] = 0;
    umax[1] = 100;
    usol[1] = sqrt(2*1e-3*1);

    nlp.setVariableBounds(&umin[0],&umax[0]);
    nlp.setPrimal(&usol[0]);

    // Allocate an NLP solver
    IpoptSolver solver(nlp);

    // Set options
    solver.setOption("exact_hessian",false);
    solver.setOption("abstol",1e-10);

    // initialize the solver
    solver.init();

    // Solve the problem
    solver.solve();

    // Get the optimal cost
    double fopt = nlp.eval_f();
    cout << "optimal cost: " << fopt << endl;

    // Get the optimal solution
    nlp.getPrimal(&usol[0]);
    cout << "optimal solution: " << usol << endl;


  }

  // Calculate the corresponding curve for the two-state model
  {

  // Parameters of the two-state model

  SX x("x");
  SX X0("X0");
  SX T0("T0");

  SX Uc = T_u; // unloading
  SX Uh = T_l; // loading

  // Temperature profile
  SX U = Uc+(Uh-Uc)*(1-erf((x-X0)/sqrt(4*D*T0)))/2;
  cout << "U = " << U << endl;

  // Calculate the energy
  Matrix E = gauss_quadrature((U-Uc)/(Uh-Uc),x,0,1); // integrate
  cout << "E = " << E << endl;

  // Make a function
  Matrix arg;
  arg << X0 << T0;
  SXFunction Efcn(arg,E);

  // Evaluate at some point 
  double argval[] = {0.5, 10};
  Efcn.setArgument(&argval[0]);
  Efcn.evaluate();
  double Eres;
  Efcn.getOutputDense(&Eres);  
  cout << "Eres = " << Eres << endl;
  
  // Calculate the exergy
  Matrix A = E - gauss_quadrature(log(U/Uc)/log(Uh/Uc),x,0,1);

  // Make a function
  SXFunction Afcn(arg,A);

  // Evaluate  
  Afcn.setArgument(&argval[0]);
  Afcn.evaluate();
  double Ares;
  Afcn.getOutputDense(&Ares);  
  cout << "Ares = " << Ares << endl;

  


  }
  
#endif  

// 
// 
// Matrix T_l = 1; // inlet temperature during loading
//   Matrix T_u = (120+273)/double(680+273); // inlet temperature during unloading
//   Matrix D = 1e-3; // diffusion coefficient













  // some constants
  enum Gas       {       N2,       O2,       Ar,      C02,       Ne,      Kr};
 
  // Mass concentration in air
  double Ck[] =  { 0.755300, 0.231400, 0.012800, 0.000450, 0.000012, 0.00003};

  // Molar mass
  double Mk[] =  { 28.01348, 31.99880, 39.94800, 44.00980, 20.17970, 83.8000};  

  // 











  return 0;
} catch (const char * str){
  cerr << str << endl;
  return 1;
}

}
