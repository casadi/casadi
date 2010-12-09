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

#ifndef INTEGRATOR_C_H 
#define INTEGRATOR_C_H 

#ifdef __cplusplus
  extern "C" {
#endif 

#include "fx_c.h"

/* Print to cerr: redirect to print to a file */
int casadi_integrator_integrate(fx_ref ref, double t_out);

/* Reset the solver and bring the time back to t0 */
int casadi_integrator_reset(fx_ref ref, int with_sens);

#ifdef __cplusplus
  }
#endif 

#endif /* INTEGRATOR_C_H  */
