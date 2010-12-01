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

#ifndef CONVECTION_DIFFUSION_HPP
#define CONVECTION_DIFFUSION_HPP

#include <integrator/sundials_interface/cvodes_integrator.hpp>
#include <casadi/stl_vector_tools.hpp>
#include <casadi/expression_tools.hpp>
#include <nlp_solver/ipopt_interface/ipopt_interface.hpp>

using namespace std;
using namespace OPTICON;

  // print?
  static bool print = false;

  // Variables
  static int NZ = 100;

  // parameters
  static SX T_l = (680+273)/double(120+273); // inlet temperature during loading
  static SX T_u = 1; // inlet temperature during unloading
  static SX D = 1e-3; // diffusion coefficient

  // discretization
  static SX L = 1;
  static SX dz = L/NZ;
  static SX dz2 = dz*dz;

  // time points
  static SX t_sunrise = 0;
  static SX t_startup = 1;
  static SX t_loading_start = 6-sqrt(7.2);
  static SX t_loading_end = 6+sqrt(7.2);
  static SX t_shutdown = 11;
  static SX t_sunset = 12;
  static SX t_next_sunrise = 24;

  // outputs discretization
  static int nt_out = 100; 
  static int nz_out = 100;

  // differential states
  static Matrix T("T",NZ), Tdot("Tdot",NZ);
  
  // Upwind finite difference approximation of dT_dz (for flow in positive and negative direction)
  static Matrix dT_dz_pos, dT_dz_neg;

  // Central difference approximation of d2T_dz2
  static Matrix d2T_dz2;

    // Calculate the derivative approximations
  void calc_Tdisc();

  // square
  double square(double x);

  void simulate(ostream& resfile, vector<double>& t_prof);

  
#endif // CONVECTION_DIFFUSION_HPP
