#include "symbolic/casadi.hpp"
#include "nonlinear_programming/newton_implicit_solver.hpp"
#include "interfaces/csparse/csparse.hpp"


using namespace CasADi;
using namespace std;

template<typename M, typename F>
void test(const CRSSparsity& sp){
  M A = M::sym("A",sp);    
  M r = M::sym("r",3);
  M x = M::sym("x",3);
  vector<M> res_in;
  res_in.push_back(x);
  res_in.push_back(r);
  res_in.push_back(A);

  F res(res_in,mul(A,x)-r);
  NewtonImplicitSolver f(res);
  f.setOption("linear_solver",CSparse::creator);
  f.setOption("ad_mode","reverse");
  f.init();
  DMatrix::ones(f.jacSparsity(0,0)).printDense();
}

int main(){

  // DMatrix A = DMatrix::ones(3,3) + DMatrix::ones(sp_diag(3)); // dense
  // DMatrix A = DMatrix::ones(sp_diag(3)); // diagonal
  DMatrix A = 2*DMatrix::ones(sp_diag(3)); 
  A(1,0) = 1; 
  A(2,0) = 1; 
  A(2,1) = 1; // lower triangular
  
  // Permute
  vector<DMatrix> tt;
  tt.push_back(A(1,Slice()));
  tt.push_back(A(2,Slice()));
  tt.push_back(A(0,Slice()));
  A = vertcat(tt);

  cout << "SX" << endl;
  test<SXMatrix,SXFunction>(A.sparsity());

  cout << "MX" << endl;
  test<MX,MXFunction>(A.sparsity());

  return 0;
}
