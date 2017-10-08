/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


/** \brief Demonstration of sparsity propagation in CasADi
 * NOTE: Example is mainly intended for developers of CasADi.
 * This example demonstrates how dependencies can be propagated through symbolic expressions,
 * a low-level feature useful when determining sparsity patterns. This functionality is mostly
 * used internally in CasADi when calculating the sparsity pattern of Jacobians and Hessians.
 *
 * \author Joel Andersson
 * \date 2012
 */

#include "casadi/casadi.hpp"

using namespace casadi;
using namespace std;


int main(){
  MX x = MX::sym("x");
  MX y = MX::sym("y");
  MX z = MX::sym("z");

  Function f = Function("f",{x,y,z},{x+y,x-y,z*x});

  int N = 133;
  
uout() << N << std::endl;
uout() << f.mapaccum("f",N,1,{{"base",4}}) << std::endl;

  return 0;
}
