#include <casadi/stl_vector_tools.hpp>
#include <casadi/fx/fx_tools.hpp>
#include <casadi/mx/mx_tools.hpp>
#include <casadi/sx/sx_tools.hpp>
#include <casadi/matrix/matrix_tools.hpp>
#include <casadi/fx/jacobian.hpp>
#include <casadi/fx/mx_function.hpp>
#include <casadi/fx/sx_function.hpp>

using namespace CasADi;
using namespace std;

// Only symbolic variables and constants
void trivial(){
  cout << "symbolic variables and constants" << endl;
  MX X("X",10);
  MX V("V");
  vector<MX> f_in(2);
  f_in[0] = X;
  f_in[1] = V;
  vector<MX> f_out(2);
  f_out[0] = X;
  f_out[1] = MX::eye(3);
  MXFunction f(f_in,f_out);
  f.init();
  vector<MX> jacX = f.jac(0);
  vector<MX> jacV = f.jac(1);
  cout << "jacX = " << jacX << endl;
  cout << "jacV = " << jacV << endl;
}

void subtraction(){
  cout << "subtraction test" << endl;
  MX X("X",10);
  MX V("V",10);
  vector<MX> f_in(2);
  f_in[0] = X;
  f_in[1] = V;
  MXFunction f(f_in,X-V);
  f.init();
  vector<MX> jacX = f.jac(0);
  vector<MX> jacV = f.jac(1);
  cout << "jacX = " << jacX << endl;
  cout << "jacV = " << jacV << endl;
  
  MXFunction f2(f_in,V-X);
  f2.init();
  vector<MX> jacX2 = f2.jac(0);
  vector<MX> jacV2 = f2.jac(1);
  cout << "jacX2 = " << jacX2 << endl;
  cout << "jacV2 = " << jacV2 << endl;
  
}

void evaluation(){
  cout << "evaluation test" << endl;
  
  // Create an SXFunction
  SXMatrix x = symbolic("x",10);
  SXMatrix y = symbolic("y");
  SXMatrix f = y*(sin(x)+x);
  cout << "f = " << f << endl;
  vector<SXMatrix> xy(2); xy[0] = x;  xy[1] = y;
  SXFunction fcn(xy,f);
  fcn.init();
  
  // Create a trivial MX function
  MX X("X",10);
  MX Y("Y");
  vector<MX> XY(2); XY[0] = X;  XY[1] = Y;
  vector<MX> F = fcn.call(XY);
  MXFunction Fcn(XY,F);
  Fcn.init();
    
  // Get the symbolic Jacobian
  vector<MX> J = Fcn.jac();
  cout << J << endl;
  
  // Create a symbolic Jacobian function
  MXFunction Jac_sym(XY,J);
  Jac_sym.init();
  
  // Use the Jacobian function for comparison
  FX Jac_old = Jacobian(Fcn);
  Jac_old.init();
  
  // Finally, create the symbolic jacobian function
  FX Jac_sx = fcn.jacobian();
  Jac_sx.init();
  
  // Give arguments to x and y
  vector<double> x0(x.size());
  for(int i=0; i<x0.size(); ++i)
    x0[i] = 1./(i+1);
  double y0 = 10;
  
  // Evaluate the three functions and compare the result
  Jac_sx.setInput(x0,0);
  Jac_sx.setInput(y0,1);
  Jac_sx.evaluate();
  cout << "Using SXFunction directly, nnz = " << Jac_sx.output().size() << endl;
  cout << Jac_sx.output() << endl;
  
  Jac_sym.setInput(x0,0);
  Jac_sym.setInput(y0,1);
  Jac_sym.evaluate();
  cout << "Using MXFunction::jac function, nnz = " << Jac_sym.output().size() << endl;
  cout << Jac_sym.output() << endl;

  Jac_old.setInput(x0,0);
  Jac_old.setInput(y0,1);
  Jac_old.evaluate();
  cout << "Using Jacobian function, nnz = " << Jac_old.output().size() << endl;
  cout << Jac_old.output() << endl;

}

void mapping(){
  cout << "mapping " << endl;
  
  // Create a trivial MX function
  MX X("X",3);
  MX Y("Y",2);
  vector<MX> XY(2); XY[0] = X;  XY[1] = Y;
  MX F = vertcat(X,Y);
  MXFunction Fcn(XY,F);
  Fcn.init();
    
  // Get the symbolic Jacobian
  MX J = Fcn.jac()[0];
  cout << J << endl;
  DMatrix JJ(J.sparsity(),1);
  JJ.printDense();
  
}

int main(){

  // Only symbolic variables
  trivial();

  // Subtraction
  subtraction();
    
  // Function evaluation
  evaluation();
    
  // Nonzero mappings
  mapping();
  
  return 0;
}
