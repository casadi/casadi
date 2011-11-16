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
  
  // Jacobians
  vector<MX> jacX = f.jac(0);
  vector<MX> jacV = f.jac(1);
  cout << "jacX = " << jacX << endl;
  cout << "jacV = " << jacV << endl;
  
  // Gradients
  vector<MX> grad0 = f.grad(0);
  vector<MX> grad1 = f.grad(1);
  cout << "grad0 = " << grad0 << endl;
  cout << "grad1 = " << grad1 << endl;
  
  
  
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

  vector<MX> g = f.grad(0);
  cout << "g = " << g << endl;
  
  MXFunction f2(f_in,V-X);
  f2.init();
  vector<MX> jacX2 = f2.jac(0);
  vector<MX> jacV2 = f2.jac(1);
  cout << "jacX2 = " << jacX2 << endl;
  cout << "jacV2 = " << jacV2 << endl;
  
  vector<MX> g2 = f2.grad(0);
  cout << "g2 = " << g2 << endl;
  
}

void evaluation(){
  cout << "evaluation test" << endl;
  
  // Create an SXFunction
  SXMatrix x = ssym("x",10);
  SXMatrix y = ssym("y");
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
  
  // Create the symbolic jacobian function
  FX Jac_sx = fcn.jacobian();
  Jac_sx.init();

  // Calculate the gradient via adjoint mode AD, source code transformation
  MX G = Fcn.grad().at(0);

  // Create a symbolic Jacobian function using the gradient
  MXFunction Jac_adj(XY,trans(G));
  Jac_adj.init();
  
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

  Jac_adj.setInput(x0,0);
  Jac_adj.setInput(y0,1);
  Jac_adj.evaluate();
  cout << "Using MXFunction::grad function, nnz = " << Jac_adj.output().size() << endl;
  cout << Jac_adj.output() << endl;


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
  IMatrix JJ(J.sparsity(),1);
  JJ.printDense();
  
  // Get the symbolic Jacobian
  MX G = Fcn.grad()[0];
  cout << G << endl;
  IMatrix GG(G.sparsity(),1);
  GG.printDense();
}

void multiplication(){
  MX X("X",4);
  MX Y("Y",4);
  vector<MX> XY(2); XY[0] = X; XY[1] = Y;
  
  MX F = mul(reshape(X,2,2),reshape(Y,2,2));
  cout << "F = " << F << endl;
  
  MXFunction Fcn(XY,F);
  Fcn.init();
  
  // Generate Jacobian
  MX J1 = Fcn.jac()[0];
  cout << "J1 = " << J1 << endl;
  MXFunction Jfcn1(XY,J1);
  Jfcn1.init();
  SXFunction Jfcn2(Jfcn1);
  cout << "J (mx) = " << endl << Jfcn2.outputSX() << endl;
  cout << "nnz(J) = " << Jfcn2.outputSX().size() << endl;
  
  // Via adjoint mode
  MX G1 = Fcn.grad()[0];
  cout << "G1 = " << G1 << endl;
  MXFunction Gfcn1(XY,G1);
  Gfcn1.init();
  SXFunction Gfcn2(Gfcn1);
  cout << "G (mx) = " << endl << Gfcn2.outputSX() << endl;
  cout << "nnz(G) = " << Gfcn2.outputSX().size() << endl;
  
  SXMatrix x = ssym("x",4);
  SXMatrix y = ssym("y",4);
  SXMatrix f = mul(reshape(x,2,2),reshape(y,2,2));
  vector<SXMatrix> xy(2); xy[0] = x; xy[1] = y;
  SXFunction fcn(xy,f);
  fcn.init();
  SXMatrix Jfcn3 = fcn.jac();
  
  cout << "J (sx) = " << endl << Jfcn3 << endl;
  cout << "nnz(J) = " << Jfcn3.size() << endl;
  
/*  cout << */
  
  
/*  cout << "f = "*/
  
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

  // Matrix multiplication
  multiplication();

  return 0;
}
