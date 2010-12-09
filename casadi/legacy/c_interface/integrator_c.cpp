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

#include "integrator_c.hpp"

using namespace CasADi;
using namespace std;

// Get a reference
Integrator& get_integrator(fx_ref ref){
  if(ref==0) throw CasadiException("get_integrator failed: null pointer");
  Integrator *r = (Integrator*)(ref);
  return *r;
}

extern "C"
int casadi_integrator_integrate(fx_ref ref, double t_out){
  try{
    get_integrator(ref).integrate(t_out);
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_integrator_reset(fx_ref ref, int with_sens){
  try{
    get_integrator(ref).reset(bool(with_sens));
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }  
}

