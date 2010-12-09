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

#include "sundials_integrator.hpp"

namespace OPTICON{

SundialsIntegrator::SundialsIntegrator(OCP_old &ocp_) : IntegratorNode_old(ocp_){
  addOption("mupper",                    OT_INTEGER,  Option()); // upper band-width of banded jacobians
  addOption("mlower",                    OT_INTEGER,  Option()); // lower band-width of banded jacobians

  addOption("linear_solver",             OT_STRING, "dense"); // "dense", "band" or "sparse"
  addOption("sparse_solver",             OT_STRING, "spgmr"); // "spgmr", "spbcg" or  "sptfqmr"
  addOption("pretype",                   OT_STRING, "none"); // "none", "left", "right", "both"

}

SundialsIntegrator::~SundialsIntegrator(){

}

void SundialsIntegrator::init(){
  IntegratorNode_old::init();
}

void SundialsIntegrator::sundialsAssert(int flag) const{
  if(flag<0) throw "Error in Sundials function";
}

} // namespace OPTICON

