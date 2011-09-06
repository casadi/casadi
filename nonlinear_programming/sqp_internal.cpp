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

#include "sqp_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include <ctime>

using namespace std;
namespace CasADi{

SQPInternal::SQPInternal(const FX& F, const FX& G, const FX& H, const FX& J){
  casadi_warning("The SQP method has not been implemented");
}


SQPInternal::~SQPInternal(){
}

void SQPInternal::init(){
  // Call the init method of the base class
  NLPSolverInternal::init();
}

void SQPInternal::evaluate(int nfdir, int nadir){
  casadi_assert(nfdir==0 && nadir==0);
}

} // namespace CasADi
