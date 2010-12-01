/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include "../sx/sx_matrix.hpp"
#include "../mx/mx.hpp"
#include "../stl_vector_tools.hpp"
#include "../fx/sx_function.hpp"
#include "../fx/mx_function.hpp"
#include "../fx/fx.hpp"
#include "../matrix.hpp"

using namespace std;
using namespace CasADi;

void test_build_tree(){

  MX x("x",3,3);
  MX y("y",3,3);  
  
  cout << "x = " << x << endl;
  cout << "y = " << y << endl;
  cout << "x+y = " << x+y << endl;
  cout << "prod(x,y) = " << prod(x,y) << endl;
  cout << "(prod(x,y)+y)(1,0) = " << (prod(x,y)+y)(1,0) << endl;
  cout << "norm_1(x+x) = " << norm_1(x+x) << endl;

  cout << "x[3] = " << x[3] << endl;
  // Evaluation
/*  MX f = (norm_1(x+x) + x[0]).function(x);
  MX f_y = f(y);*/
//   cout << "f = " << f << endl;
//   cout << "f_y = " << f_y << endl;

  // Transpose
  cout << "X' = " << trans(x) << endl;
  cout << "X'' = " << trans(trans(x)) << endl;

  // Construct a block matrix
  MX bmx = vertcat(horzcat(prod(x,x), y),horzcat(y+x, x));
  cout << "bmx = " << bmx << endl;
  

}


void test_runtime(){
  // Make a nonlinear function
  SX z("z");
  SX g = sin(z + 2*sqrt(z));
  SXFunction G(z,g);
  G.setOption("ad_order",1);
  G.init();
  
  // Create a matrix function
  MX X("X",3,3); // matrix variable
  MX Y("Y"); // 1-by-1 matrix variable
  MX F = prod(X,X);

    F.print(); // TODO: Mystical bug

  MX F1 = F;
  MX F2 = F1+F;
/*  F += F;*/
  F = F2;
  F1.print(); // TODO: Mystical bug
  F2.print(); // TODO: Mystical bug
  
  
  cout << "G = " << G << endl; // TODO: Mystical bug
  cout << Y << endl;
  MX GY = G(Y);
  GY.print(); cout << endl;// TODO: Mystical bug
  G.print();  cout << endl;// TODO: Mystical bug
  cout << "GY = " << GY << endl; // TODO: Mystical bug
  F *= G(Y); // function evaluation
  F += X;

  // Print the function
  cout << "F = " << F << endl; 

  // Create a matrix function
#ifdef WITH_CPP0X
  vector<MX> arg = {X,Y};  // C++0x style
#else
  vector<MX> arg(2);
  arg[0] = X;
  arg[1] = Y;
#endif
  MXFunction rt(arg,F);

  // Print the algorithm
  cout << "algorithm:" << endl;
  rt->printAlgorithm();

  // Give a value to x
  double xval[9] = {1,2,3,4,5,6,7,8,9};
  for(int i=0; i<9; ++i)
    rt.input().data()[i] = xval[i];

  // Give a value to y
  double yval = 100;
  rt.input().dataF()[0] = yval;
  
  // Print the values
  cout << "values:" << endl;
  rt->printValues();
  
  // Evaluate the algorithm
  rt.evaluate();

  // Print the values again
  cout << "values:" << endl;
  rt->printValues();

  // Give seed to the arguments
  double xseed1[9] = {1,0,0,0,0,0,0,0,0};
  for(int i=0; i<9; ++i)
    rt.input().dataF()[i] = xseed1[i];

  // Give a value to y
  double yseed1 = 0;
  rt.input(1).dataF()[0] = yseed1;

  // Print the seeds
  cout << "seeds (AD fwd):" << endl;
  rt->printValues(1);

  // Evaluate using AD forward
  rt.evaluate(1,0);

  // Print the seeds again
  cout << "dres (AD fwd):" << endl;
  rt->printValues(1);

  // Clear the seeds
  rt->clear(1);
  
  // Give a seed to the last element
  double fseed[] = {1,0,0,0,0,0,0,0,0};
  vector<double>& op = rt.output().dataF();
  for(int i=0; i<op.size(); ++i)
    op[i] = fseed[i];

  cout << "seeds (AD adj):" << endl;
  rt->printValues(1);

  // Do adjoint AD
  rt.evaluate(0,1);

  // Print the seeds again
  cout << "dres (AD adj):" << endl;
  rt->printValues(1);

  
    

}

int main(){

try{

//   SX xx("x");
//   SX yy = 2;
//   Matrix xxyy = {{xx,yy,3},{2,xx+2,4}};
//   Matrix xxyy2 = {3,22,3};
//   
//   cout << xxyy << endl;
//   cout << xxyy2 << endl;
//   return 0;

  
 test_build_tree();

 test_runtime();

  // Create a function
  SX x("x");
  SX f = x*x + sin(sqrt(x)); 

  SXFunction fcn(x,f);
  fcn.setOption("ad_order",1);
  fcn.init();

  // print algorithm
  fcn->printAlgorithm();

  // result
  double res;

  // Evaluate
  double x0 = 1.2;
  fcn.input().data()[0] = x0;
  fcn.evaluate();
  cout << "res = " << fcn.output().data() << endl;

  // Adjoint AD
  double df0 = 1;
  fcn.output().data()[0] = df0;
  fcn.evaluate(0,1);
  res = fcn.input().dataA()[0];
  cout << "dx (adj) = " << res << endl;

  // Forward AD
  double dx0 = 1;
  fcn->clear(1);
  fcn.input().dataF()[0] = dx0;
  fcn.evaluate(1,0); 
  cout << "dx (fwd) = " << fcn.input().dataF() << endl;

 
} catch (const char * str){
  std::cerr << str << std::endl;
  return 1;
}


}
