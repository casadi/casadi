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

bool on_off = true;
int startup_time = 2;
bool init_free = true;

// cost functions
double ls_cost(double phi, double phi_ref){
  double diff = phi-phi_ref;
  return diff*diff;
}

double energy_cost(double phi, double phi_ref){
  // efficiency of steam cycle
  double eta = (0.80-phi)/(0.80-0.25);
  eta = 1 - eta*eta;

  // usefull efficiency
  return -phi*eta;
}

typedef double (*cost_function)(double, double);
cost_function calc_cost[] = {ls_cost,energy_cost};


  // Tables
  vector< vector< int > > transition;
  vector< vector< double > > transition_cost;
  vector< double > phi_ref(NK,0);
  vector< double > phi_r(NK,0);
  vector<double> u_disc(n_u_disc);
  vector<int> x_int(x_dim);
  vector< vector<double> > x_disc(x_dim);
  int n_x_disc_total;
  
  // Calculate tables
  void calc_tables(){

      // read tabulated data from disk
  ifstream indata;
  indata.precision(12);

  indata.open ("phir_meas.txt");
  read_matlab(indata,phi_r);
  indata.close(); 
      
  indata.open ("transition.txt");
  read_matlab(indata,transition);
  indata.close();
    
  indata.open ("transition_cost.txt");
  read_matlab(indata,transition_cost);
  indata.close();

#if 0
  vector< vector<double> > transition_misfit;
  indata.open ("transition_misfit.txt");
  read_matlab(indata,transition_misfit);
  indata.close();
  
  vector< vector<double> > transition_roundoff;
  indata.open ("transition_roundoff.txt");
  read_matlab(indata,transition_roundoff);
  indata.close();

  #endif
     
  indata.open("x_int.txt");
  read_matlab(indata,x_int);
  indata.close();

  // MISC
  for(int i=0; i<n_u_disc; ++i)
    u_disc[i] = (i-n_u_pos)/double(n_u_pos);

  // Total number of spatial discretizations
  n_x_disc_total = 1;
  for(int i=0; i<x_dim; ++i)
    n_x_disc_total *= n_x_disc[i];

  // Value for x for each spatial discretization
  for(int i=0; i<x_dim; ++i){
    x_disc[i].resize(n_x_disc[i]);

    for(int j=0; j<n_x_disc[i]; ++j){
      x_disc[i][j] = x_min[i] + (x_max[i]-x_min[i])*j/double(n_x_disc[i]-1);
    }
  }  
    
  // receiver flow and reference trajectory
  for(int k=0; k<NK; ++k){
 
    // k during one day
    int kk = k % (24*nk_per_hour);

//     // receiver mass flow
//     if(kk<12*nk_per_hour)
//       phi_r[k] = (kk*(12*nk_per_hour-kk))/double(12*nk_per_hour*12*nk_per_hour/4);

    // Heat exchanger flow setpoint
    if(kk>=nk_per_hour && kk<nk_per_hour*11){       // plant turned on
      phi_ref[k] = 0.8;
    }
  }
  }
  
  // Transition function
  int transition_fcn(int i, int j, int k){
    int time_on = i / n_x_disc_total;
    i %= n_x_disc_total;

    bool plant_on = j >= n_u_disc;
    j %= n_u_disc;

    int new_time_on = min((time_on+1)*plant_on,startup_time);
    return transition[i][j] + new_time_on*n_x_disc_total;
  }
  
  // Objective function contribution
  double stage_cost(int i, int j, int k){
      int time_on = i / n_x_disc_total;
      i %= n_x_disc_total;
   
      bool plant_on = j >= n_u_disc;
      j %= n_u_disc;
              
      bool generates_electricity = plant_on && time_on==startup_time;
      
      // flow to boiler
      double phi = phi_r[k]-u_disc[j];

      // efficiency of steam cycle
      double eta = (0.80-phi)/(0.80-0.25);
      eta = 1 - eta*eta;
           
      // Cost at time k
      double cost = -100 * generates_electricity * phi*eta; // usefull energy
      
      // penalize mass flow through storage
      cost += 0.1*u_disc[j]*u_disc[j];

// alternative objective function
      if(obj==LEAST_SQUARE){
        double diff = phi-phi_ref[k];
	cost = diff*diff;
      }
      
      // Add the cost resulting from the transition (zero or infinity)
      cost += transition_cost[i][j];

      // Infinite cost if flow to heat exchanger is negative
      if(phi<0)
	cost += 1000000; 	// cost += infty;

      // Artificial cost for turning the plant on or off
//      cost += change_on_off*10;
      
      // Infinite cost if tries to run plant with phi<0.25
      if((plant_on && phi<0.25) || phi>1)
	cost += infty;
      
      return cost;
  }

void get_u(vector<double>& u_opt, int j){
  if(on_off){
    u_opt.resize(2);    
    bool plant_on = j>=n_u_disc;
    j %= n_u_disc;
    u_opt[0] = u_disc[j];
    u_opt[1] = plant_on;    
  } else {
    u_opt.resize(1);
    u_opt[0] = u_disc[j];  
  }
}

void get_x(vector<double>& x_opt, int i){
  i %= n_x_disc_total; // discregard on-off  
  x_opt.resize(x_dim);    
  for(int j=x_dim-1; j>=0; --j){
      x_opt[j] = x_disc[j][i % n_x_disc[j]];
      i /= n_x_disc[j];
  }
}

int main(){
  try{
       
  clock_t time0 = clock();
  cout << "program started " << endl;
     
  // Calculate tables
  calc_tables();

  // Enumeration of the state, control and time
  int n_i, n_j, n_k;  
  if(on_off){ // allow the plant to turn on-off
    n_i = n_x_disc_total*(1+startup_time);
    n_j = n_u_disc*2;
    n_k = NK;
  } else {
    n_i = n_x_disc_total;
    n_j = n_u_disc;
    n_k = NK;
  }
  
  // Cost to arrive function
  vector< vector<double> > J(n_k+1);
  for(int k=0; k<n_k+1; ++k){
     J[k].resize(n_i);
      for(int i=0; i<n_i; ++i)
	J[k][i] = infty;
  }

if(init_free){
  // initial condition free
  for(int i=0; i<n_i; ++i)
    J[0][i] = 0;
  
} else {
  // Initial condition fixed
  x_int[0] = n_x_disc[0]/2;
  x_int[1] = 10;
  int init_ind = x_int[0];
  for(int i=1; i<x_dim; ++i){
    init_ind = init_ind*n_x_disc[i] + x_int[i];
  }
  J[0][init_ind] = 0;
}
  
  // Optimal control trajectory table
  vector< vector<int> > u_table(n_k);
  for(int k=0; k<n_k; ++k)
    u_table[k].resize(n_i);
    
  // state trajectory table
  vector< vector<int> > x_table(n_k);
  for(int k=0; k<n_k; ++k)
    x_table[k].resize(n_i);
      
  // Progress in percent
  int progress = 0;
    
  // DP algorithm
  for(int k=0; k<n_k; ++k){

    // Print progress
    if((100*k)/n_k>progress){
      progress = (100*k)/n_k;
      cout << progress << " \%" << endl;
    }

    // Loop over the points in state space
    for(int i=0; i<n_i; ++i){

      // Loop over all possible control actions
      for(int j=0; j<n_j; ++j){

	// Evaluate the dynamic equation
	int i_new = transition_fcn(i,j,k);

	// Evaluate the cost function
	double cost = J[k][i] + stage_cost(i,j,k);
	
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
  for(int i=0; i<n_i; ++i)
    if(J[n_k][i] < c_min){
      c_min = J[n_k][i];
      i_min = i;
    }

  cout << "terminal value for state = " << i_min << ", optimal cost = " << J[n_k][i_min] << endl;
  
  // Recover the optimal policy
  cout << "starting rollback" << endl;

  // optimal control trajectory
  vector<int> u_opt(n_k);

  // optimal state trajectory
  vector<int > x_opt(n_k);

  // rollback
  int i_opt = i_min;
  for(int k=n_k-1; k>=0; --k){
    x_opt[k] = i_opt;
    u_opt[k] = u_table[k][i_opt];
    i_opt = x_table[k][i_opt];
  }

  // Get the corresponding numerical values
  vector<vector<double> > u_opt_val(n_k);
  vector<vector<double> > x_opt_val(n_k);
  for(int k=0; k<n_k; ++k){
    get_u(u_opt_val[k],u_opt[k]);
    get_x(x_opt_val[k],x_opt[k]);
  }
      
   clock_t time3 = clock();
   cout << "program finished after " << (time3-time0)/double(CLOCKS_PER_SEC) << " seconds." << endl;

   // Create a file for saving the results
  ofstream resfile;

   // save to disk
  resfile.open ("optimal_u.txt");
  write_matlab(resfile,u_opt_val);
  resfile.close();

  resfile.open ("optimal_x.txt");
  write_matlab(resfile,x_opt_val);
  resfile.close();
   
  resfile.open ("phi_r.txt");
  write_matlab(resfile,phi_r);
  resfile.close();
   
  resfile.open ("phi_ref.txt");
  write_matlab(resfile,phi_ref);
  resfile.close();



  return 0;
} catch (exception &e){
  cerr << e.what() << endl;  
  return 1;
} catch (const char * str){
  cerr << str << endl;
  return 1;
}

}
