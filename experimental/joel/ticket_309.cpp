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

#include "casadi/casadi.hpp"
#include <ctime>

using namespace std;
using namespace CasADi;

int main() {
  int N = 2000;
  int n = 100;
  SXMatrix x = ssym("x",N);
  SXFunction f(x,pow(x,SXMatrix(2)));
  f.init();

  MX X = msym("X",N,n);
  MX total = 0;
  for(int i=0; i<n; ++i){
    MX Xcol = X(range(X.size1()),i);
    total += f.call(vector<MX>(1,Xcol))[0];
  }
  MXFunction F(X,total);
  F.setOption("verbose",true);
  F.init();
  
  SXFunction FF(F);
  FF.setOption("verbose",true);
  FF.init();
  
  std::cout << "Please have a coffee break" << std::endl;
  CRSSparsity Jsp = FF.jacSparsity();
  std::cout << "Hooray, finished before the universe ended" << std::endl;
  std::cout << "Jsp = " << Jsp << std::endl;

  return 0;
}
