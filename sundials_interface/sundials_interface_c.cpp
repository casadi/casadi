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

#include "sundials_interface_c.h"
#include "cvodes_integrator.hpp"
#include "idas_integrator.hpp"
#include "casadi/c_interface/fx_c.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>

using namespace CasADi::Sundials;
using namespace std;

extern "C"
int casadi_cvodes_integrator(fx_ref fcn, fx_ref ffcn){
 try{
    get_fx(fcn) = CVodesIntegrator(get_fx(ffcn));
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_idas_integrator(fx_ref fcn, fx_ref ffcn){
 try{
    get_fx(fcn) = IdasIntegrator(get_fx(ffcn));
    return 0;
  } catch(...){
    return 1;
  }
}

