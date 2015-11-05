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


#include "core/casadi.hpp"
#include "interfaces/ooqp/ooqpsol.hpp"
#include <ctime>

using namespace std;
using namespace casadi;

int main() {
  DMatrix H(2,2,0);
  H(0,0) = 1.0;
  H(1,1) = 0.5;
  
  DMatrix A(1,2,0);
  A(0,0) = 1.0;
  A(0,1) = 1.0;
  
  DMatrix g(2,1,0);
  g(0) = 1.5;
  g(1) = 1.0;
  
  DMatrix lb(2,1,0);
  lb(0) = 0.5;
  lb(1) = -2.0;
  
  DMatrix ub(2,1,0);
  ub(0) = 5.0;
  ub(1) = 2.0;
  
  DMatrix lbA(1,1,0);
  lbA(0) = -1.0;
  
  DMatrix ubA(1,1,0);
  ubA(0) = 2.0;
  
  for(int rep=0; rep<2; ++rep){
    OOQpsol qpsol(H.sparsity(), A.sparsity());
    qpsol.init();
    qpsol.setInput(A,"a");
    qpsol.setInput(H,"h");
    qpsol.setInput(g,"g");
    qpsol.setInput(lb,"lbx");
    qpsol.setInput(ub,"ubx");
    qpsol.setInput(lbA,"lba");
    qpsol.setInput(ubA,"uba");
    qpsol.evaluate();
    qpsol.evaluate();
    qpsol.evaluate();
  }
    
  return 0;
}
