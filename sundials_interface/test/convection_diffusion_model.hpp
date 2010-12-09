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

#ifndef CONVECTION_DIFFUSION_MODEL_HPP
#define CONVECTION_DIFFUSION_MODEL_HPP

#include "convection_diffusion_parameters.hpp"
#include <casadi/expression_tools.hpp>
#include <integrator/sundials_interface/cvodes_integrator.hpp>
#include <integrator/sundials_interface/cvodes_interface.hpp>

using namespace std;
using namespace OPTICON;


class CDModel{
  public:
  CDModel(int NZ);
     
  // Just simulate a scenario
  void simulate(ostream& resfile, vector<double>& t_prof);
  void simulate_new(ostream& resfile, vector<double>& t_prof);

  CVodesInterface make_integrator_new(int dir);
  CvodesIntegrator make_integrator(int dir);
  
  // differential states
  Matrix T, Tdot;
  
  // Upwind finite difference approximation of dT_dz (for flow in positive and negative direction)
  Matrix dT_dz_pos, dT_dz_neg;

  // Central difference approximation of d2T_dz2
  Matrix d2T_dz2;
  
  // parameters
  SX T_l; // inlet temperature during loading
  SX T_u; // inlet temperature during unloading
  SX D;   // diffusion coefficient

  // Initial temperature in the storage
  Matrix T0;

  // Initial value problems corresponding to positive, negative and zero flow
  OCP_old ivp_pos, ivp_neg, ivp_zero;  

  // Control (flow to storage)
  SX u;

  // Maximal u
  SX u_max; // fully loading/unloading takes 1.5 hours

};

#endif // CONVECTION_DIFFUSION_MODEL_HPP

