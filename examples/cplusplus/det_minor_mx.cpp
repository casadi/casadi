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

#include "symbolic/mx/mx_tools.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/fx/sx_function.hpp"
#include "symbolic/fx/mx_function.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/fx.hpp"
#include <ctime>

using namespace std;
using namespace CasADi;

MX det(const MX& a);
MX getMinor(const MX &x, int i, int j);
MX cofactor(const MX &x, int i, int j);
MX adj(const MX& a);

vector<FX> detFuncs;
const int N_EXPAND_LIMIT = 9;

MX my_det(const MX& a){
  int n = a.size1(); // Dimension
  if(n==1) return a;

  // Return expression
  MX ret = 0;

  // We expand the matrix along the first column
  for(int i=0; i<n; ++i){
    // Sum up the cofactors
    ret += a(i,0)*cofactor(a,i,0);
  }
  return ret;
}

MX getMinor(const MX &x, int i, int j){
  int n = x.size1();

  // Trivial return if scalar
  if(n==1) return 1;

  // Remove row i and column j
  MX M(n-1,n-1);

   for(int i1=0; i1<n; ++i1)
       for(int j1=0; j1<n; ++j1){
           if(i1 == i || j1 == j)
              continue;

            int i2 = (i1<i)?i1:i1-1;
            int j2 = (j1<j)?j1:j1-1;
    
            M(i2,j2) = x(i1,j1);
       }

  // Create determinant function, if there isn't already one for the current size
  if(detFuncs[n-1].isNull()){
    // Argument
    MX A("A",(n-1)*(n-1));
    MX dA = my_det(reshape(A,n-1,n-1));
    MXFunction detFun(A,dA);
    detFun.init();
    
    // If small size, expand to SXMatrix
    if(n<N_EXPAND_LIMIT){
      detFuncs[n-1] = SXFunction(detFun);
      detFuncs[n-1].init();
    } else {
      detFuncs[n-1] = detFun;
    }
  }

  // Make a function call
  vector<MX> dM = detFuncs[n-1].call(flatten(M));
  return dM[0];
}

MX cofactor(const MX &x, int i, int j){

    // Calculate the i,j minor
    MX minor_ij = getMinor(x,i,j);

    // Calculate the cofactor
    int sign_i = 1-2*((i+j) % 2);

    return sign_i * minor_ij;
}

MX adj(const MX& a){
  int n = a.size1();
  casadi_assert_message(n == a.size2(),"adj: matrix must be square");

  // Cofactor matrix
  MX C(n,n);
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      C(i,j) = cofactor(a,i,j);
  
  return trans(C);
}



int main() {
  int n = 10;
  detFuncs.resize(n+1);
  
  MX A("A",n*n);
  MX dA = my_det(reshape(A,n,n));
  cout << "constructed"<< endl;
  
  MXFunction detFun(A,dA);
  cout << "constructor called"<< endl;
  detFun.init();
  cout << "initialized"<< endl;
  
  // Some random matrix
  DMatrix A0(n*n,1,0);
  for(int i=0; i<n*n; ++i)
    A0[i] = 0.1 + sqrt(double(i*1323 % 124));
  cout << A0 << endl;
  
  detFun.setInput(A0);
  cout << "set input"<< endl;
  
  clock_t time1 = clock();
  detFun.evaluate();
  clock_t time2 = clock();
  cout << "evaluated no derivatives: "<< (double(time2 - time1)/CLOCKS_PER_SEC*1000) << " ms" << endl;
  
  cout << detFun.output() << endl;
  
/*  double det_num = my_det(reshape(A0,n,n));
  cout << det_num << endl;
  cout << "evaluated numerically "<< endl;*/
  
  detFun.setAdjSeed(1.);

  time1 = clock();
  detFun.evaluate(0,1);
  time2 = clock();
  cout << "evaluated adj: "<< (double(time2 - time1)/CLOCKS_PER_SEC*1000) << " ms" << endl;
  
  cout << detFun.adjSens() << endl;
  
  return 0;
}

