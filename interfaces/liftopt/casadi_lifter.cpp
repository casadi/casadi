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

#include "casadi_lifter.hpp"
using namespace std;
namespace CasADi{
  namespace Interfaces{

void CasadiLifter::liftfun(double *v, int n, void *user_data){
  static_cast<CasadiLifter*>(user_data)->setNode(v,n);
}

long CasadiLifter::evalUserFcn ( liftopt::TLifterArgs<double>& args )
{
  MXFunction& F = interface_->fcn_;
  if(m_activeProblemFormulation ==  Problem_Residual ){
    F.setInput(args.u,LO_U);
    F.evaluate();
    F.getOutput(args.objRes,LO_OBJRES);
    F.getOutput(args.eq,LO_EQ);
  }
  if(m_activeProblemFormulation ==  Problem_LagGrad){
    F.setInput(args.u,LO_U);
    F.setInput(args.lambda,LO_LAMBDA);
    F.adjSeed(LO_OBJRES).setAll(0);
    F.adjSeed(LO_EQ).setAll(0);
    F.setAdjSeed(0.0,LO_OBJ);
    F.setAdjSeed(1.0,LO_LAGFCN);
    F.evaluate(0,1);
    F.getOutput(args.eq,LO_EQ);
    F.getOutput(args.objVal,LO_OBJ);
    F.getAdjSens(args.laggrad);
  }
  return 0;
}

  } // namespace Interfaces
} // namespace CasADi
