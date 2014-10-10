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


#include "core/std_vector_tools.hpp"
#include "core/function/sx_function.hpp"
#include "core/function/mx_function.hpp"
#include "core/function/external_function.hpp"
#include "core/sx/sx_tools.hpp"
#include "core/mx/mx_tools.hpp"
#include <ctime>
#include <stdlib.h>

using namespace std;
using namespace casadi;

// Function to be called if the size reaches k
Function det_n_min; 
int n_min;

// Forward declaration
template<class M>
M getMinor(const M& x, int i, int j);

// Determinant
template<class M>
M det(const M& a){
  int n = a.size1(); // Dimension
  if(n==n_min) return det_n_min.call(a)[0]; // Quick return
  if(n==1) return a; // Quick return
  
  // Expand along the first column
  M ret = 0;
  for(int i=0; i<n; ++i){
    int sign = i % 2 == 1 ? 1 : -1;
    ret += a(i,0) * sign * getMinor(a,i,0);
  }
  return ret;
}

// Minor
template<class M>
M getMinor(const M& x, int i, int j){
  int n = x.size1(); // Dimension

  // Remove row i and column j from x
  M m(n-1,n-1); // (n-1)-by-(n-1) matrix
  for(int i1=0; i1<n; ++i1){
    for(int j1=0; j1<n; ++j1){
      if(!(i1 == i || j1 == j)){
        int i2 = i1<i ? i1 : i1-1;
        int j2 = j1<j ? j1 : j1-1;
        m(i2,j2) = x(i1,j1);
      }
    }
  }
  return det(m);
}

clock_t time1, time2;

DMatrix x0(int n){
  DMatrix A0(n,n,0);
  for(int i=0; i<n*n; ++i)
    A0[i] = 0.1 + sqrt(double(i*1323 % 124));
    
  return A0;
}

void test_mx(){
  int n=8;
  cout << "Testing pure MX only, n = " << n << endl;
  
  MX x("x",n,n);
  time1 = clock();
  MX detx = det(x);
  time2 = clock();
  double t = double(time2 - time1)/CLOCKS_PER_SEC;
  cout << "time to create graph: " << (t*1e3) << " ms, " << endl;
  
  time1 = clock();
  MXFunction f(x,detx);
  f.init();
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  int nops = f.countNodes();   // Number of elementary operations
  cout << "time to create function: " << (t*1e3) << " ms, " << (t*1e9/nops) << " ns per elementary operation" << endl;
  
  // Some "random" matrix as input
  f.setInput(x0(n));

  // Evaluate without derivatives
  int nrep = 100;
  time1 = clock();
  for(int i=0; i<nrep; ++i){
    f.evaluate();
  }
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  t /= nrep;
  cout << "evaluated without derivatives: " << (t*1e3) << " ms, " << (t*1e9/nops) << " ns per elementary operation" << endl;

  // Evaluate with adjoint mode AD
  f.setAdjSeed(1.0);
  time1 = clock();
  for(int i=0; i<nrep; ++i){
    f.evaluate(0,1);
  }
  time2 = clock();
  //cout << "adjoint sensitivities = " << f.adjSens() << endl;
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  t /= nrep;
  cout << "evaluated adjoint derivatives: " << (t*1e3) << " ms, " << (t*1e9/nops) << " ns per elementary operation" << endl;
}

void test_sx(){
  int n=8;
  cout << "Testing pure SX only, n = " << n << endl;
  
  SX x = ssym("x",n,n);
  time1 = clock();
  SX detx = det(x);
  time2 = clock();
  double t = double(time2 - time1)/CLOCKS_PER_SEC;
  cout << "time to create graph: " << (t*1e3) << " ms, " << endl;
  
  time1 = clock();
  SXFunction f(vec(x),detx);
  f.setOption("live_variables",true);
  f.setOption("topological_sorting","depth-first");
  f.init();
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  int nops = f.countNodes();   // Number of elementary operations
  cout << "number of operations " << nops << endl;
  cout << "time to create function: " << (t*1e3) << " ms, " << (t*1e9/nops) << " ns per elementary operation" << endl;

  // Create a new function without live variables (a bug prevents them from working together with symbolic calculations)
  time1 = clock();
  SXFunction f_nolive(vec(x),detx);
  f_nolive.init();
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  cout << "time to create function (no live): " << (t*1e3) << " ms, " << (t*1e9/nops) << " ns per elementary operation" << endl;

  // Differentiating symbolically
  time1 = clock();
  SX grad = f_nolive.grad();
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  cout << "time to create gradient: " << (t*1e3) << " ms, " << (t*1e9/nops) << " ns per elementary operation" << endl;

  // Differentiating symbolically
  time1 = clock();
  SXFunction g(vec(x),grad);
  g.setOption("live_variables",true);
  g.setOption("topological_sorting","depth-first");
  g.init();
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  int nops_grad = g.countNodes();   // Number of elementary operations
  cout << "number of operations in gradient function: " << nops_grad << endl;
  cout << "time to create gradient function: " << (t*1e3) << " ms, " << (t*1e9/nops_grad) << " ns per elementary operation" << endl;

  // Some "random" matrix as input
  f.setInput(vec(x0(n)));

  // Evaluate without derivatives
  int nrep = 1000;
  time1 = clock();
  for(int i=0; i<nrep; ++i){
    f.evaluate();
  }
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  t /= nrep;
  cout << "evaluated without derivatives: " << (t*1e3) << " ms, " << (t*1e9/nops) << " ns per elementary operation" << endl;

  // Evaluate with adjoint mode AD, OO
  f.setAdjSeed(1.0);
  time1 = clock();
  for(int i=0; i<nrep; ++i){
    f.evaluate(0,1);
  }
  time2 = clock();
  //cout << "adjoint sensitivities = " << f.adjSens() << endl;
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  t /= nrep;
  cout << "evaluated adjoint derivatives (OO): " << (t*1e3) << " ms, " << (t*1e9/nops) << " / " << (t*1e9/nops_grad) << " ns per elementary operation" << endl;

  // Evaluate with adjoint mode AD, SCT
  g.setInput(f.input());
  time1 = clock();
  for(int i=0; i<nrep; ++i){
    g.evaluate();
  }
  time2 = clock();
  //cout << "adjoint sensitivities = " << f.adjSens() << endl;
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  t /= nrep;
  cout << "evaluated adjoint derivatives (SCT): " << (t*1e3) << " ms, " << (t*1e9/nops) << " / " << (t*1e9/nops_grad) << " ns per elementary operation" << endl;

  
  cout << "res SCT = " << g.output().data() << endl;
  cout << "res OO  = " << f.adjSens().data() << endl;

  // Generate c-code
  g.generateCode("det_minor_generated.c");
  time1 = clock();
  int flag;
  flag = system("time gcc -fPIC -shared det_minor_generated.c -o det_minor_generated_no_opt.so");
  flag = system("time gcc -fPIC -shared -O3 det_minor_generated.c -o det_minor_generated_O3.so");
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  cout << "compilation took: " << (t*1e3) << " ms, " << endl;

  #if 0
  ExternalFunction e("./det_minor_generated_no_opt.so");
  e.init();

  // Evaluate with adjoint mode AD, code gen
  e.setInput(f.input());
  time1 = clock();
  for(int i=0; i<nrep; ++i){
    e.evaluate();
  }
  time2 = clock();
  //cout << "adjoint sensitivities = " << f.adjSens() << endl;
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  t /= nrep;
  cout << "evaluated adjoint derivatives (code gen): " << (t*1e3) << " ms, " << (t*1e9/nops) << " / " << (t*1e9/nops_grad) << " ns per elementary operation" << endl;
  
  // Evaluate with adjoint mode AD, code gen, optimization
  e = ExternalFunction("./det_minor_generated_no_opt.so");
  e.init();

  e.setInput(f.input());
  time1 = clock();
  for(int i=0; i<nrep; ++i){
    e.evaluate();
  }
  time2 = clock();
  //cout << "adjoint sensitivities = " << f.adjSens() << endl;
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  t /= nrep;
  cout << "evaluated adjoint derivatives (code gen O3): " << (t*1e3) << " ms, " << (t*1e9/nops) << " / " << (t*1e9/nops_grad) << " ns per elementary operation" << endl;
  
  #endif



}

void test_sx_mx(){
  
  int n=9;
  
  SX x = ssym("x",n,n);
  time1 = clock();
  SX detx = det(x);
  time2 = clock();
  double t = double(time2 - time1)/CLOCKS_PER_SEC;
  cout << "time to create graph: " << (t*1e3) << " ms, " << endl;
  
  time1 = clock();
  SXFunction f(x,detx);
  f.setOption("number_of_fwd_dir",1);
  f.setOption("live_variables",true);
  f.setOption("topological_sorting","depth-first");
  f.init();
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  int nops = f.countNodes();   // Number of elementary operations
  cout << "number of operations " << nops << endl;
  cout << "time to create function: " << (t*1e3) << " ms, " << (t*1e9/nops) << " ns per elementary operation" << endl;

  // Some "random" matrix as input
  f.setInput(x0(n));

  // Evaluate without derivatives
  int nrep = 10;
  time1 = clock();
  for(int i=0; i<nrep; ++i){
    f.evaluate();
  }
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  t /= nrep;
  cout << "evaluated without derivatives: " << (t*1e3) << " ms, " << (t*1e9/nops) << " ns per elementary operation" << endl;

  // Evaluate with adjoint mode AD, OO
  f.setAdjSeed(1.0);
  time1 = clock();
  for(int i=0; i<nrep; ++i){
    f.evaluate(0,1);
  }
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  t /= nrep;
  cout << "evaluated adjoint derivatives (OO): " << (t*1e3) << " ms, " << (t*1e9/nops) << " ns per elementary operation" << endl;

  // Use this function for "small" matrices
  det_n_min = f; 
  n_min = n;
  
  // Now lets try a larger matrix
  n = 14;
  cout << "Testing SX and MX only, n = " << n << endl;

  MX x_mx("x",n,n);
  time1 = clock();
  MX detx_mx = det(x_mx);
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  cout << "time to create graph: " << (t*1e3) << " ms, " << endl;
  
  time1 = clock();
  MXFunction f_mx(x_mx,detx_mx);
  f_mx.init();
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  nops = f_mx.countNodes();   // Number of elementary operations
  cout << "time to create function: " << (t*1e3) << " ms, " << (t*1e9/nops) << " ns per elementary operation" << endl;
  
  // Some "random" matrix as input
  f_mx.setInput(x0(n));

  // Evaluate without derivatives
  nrep = 1;
  time1 = clock();
  for(int i=0; i<nrep; ++i){
    f_mx.evaluate();
  }
  time2 = clock();
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  t /= nrep;
  cout << "evaluated without derivatives: " << (t*1e3) << " ms, " << (t*1e9/nops) << " ns per elementary operation" << endl;

  // Evaluate with adjoint mode AD
  f_mx.setAdjSeed(1.0);
  time1 = clock();
  for(int i=0; i<nrep; ++i){
    f_mx.evaluate(0,1);
  }
  time2 = clock();
  //cout << "adjoint sensitivities = " << f.adjSens() << endl;
  t = double(time2 - time1)/CLOCKS_PER_SEC;
  t /= nrep;
  cout << "evaluated adjoint derivatives: " << (t*1e3) << " ms, " << (t*1e9/nops) << " ns per elementary operation" << endl;
}

int main(){
//  test_mx();
//  test_sx();
  test_sx_mx();
  
/*  // Expand in scalar operations
  SXFunction f_sx(f);
  f_sx.init();
  f_sx.setInput(A0);
  
  time1 = clock();
  f_sx.evaluate(0,1);
  time2 = clock();
  double t_sx = double(time2 - time1)/CLOCKS_PER_SEC;
  cout << "evaluated adjoint derivatives (SX): " << (t_sx*1e3) << " ms, " << (t_sx*1e9/nops) << " ns per elementary operation" << endl;*/
  
  
  return 0;
}
