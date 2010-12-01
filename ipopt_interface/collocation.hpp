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

#ifndef COLLOCATION_HPP
#define COLLOCATION_HPP

#include <vector>
#include "ipopt_interface.hpp"
#include <casadi/fx/sx_function.hpp>
#include <acado_modelling/ocp_tools.hpp>

namespace IpoptInterface{
    using namespace OPTICON;


class Collocation{
  public:

    // default constructor
    Collocation();

    void optimize();
  
    // Degree of interpolating polynomial
    int K;
  
    // Number of finite elements
    int N;

    // Derivative of the state vector
    SXMatrix xdot;
    SXFunction ffcn;

    SXMatrix c;
    SXFunction cfcn;

    SXMatrix m;
    SXFunction mfcn;

        CollocationPoints cp;               

    // Final time
    double tf;
    
    SXMatrix x;
    vector<double> x_lb, x_ub, x_init;

    SXMatrix u;
    vector<double> u_lb, u_ub, u_init;
                
    SX t;
    
  // Collocated times
    vector<vector<double> > T;

    vector<double> x_opt;
};

} // namespace IpoptInterface

#endif // COLLOCATION_HPP
