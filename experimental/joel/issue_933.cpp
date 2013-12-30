#include "symbolic/casadi.hpp"
#include "nonlinear_programming/newton_implicit_solver.hpp"
#include "interfaces/csparse/csparse.hpp"


using namespace CasADi;
using namespace std;

void printBinary(bvec_t v){
  for(int k=0; k<bvec_size; ++k){
    if(k%4==0) cout << " ";
    if(v & (bvec_t(1)<<k)){
      cout << 1;
    } else {
      cout << 0;
    }
  }
  cout << endl;
}

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
  res.init();

  bvec_t* i0 = reinterpret_cast<bvec_t*>(res.input(0).ptr());
  bvec_t* i1 = reinterpret_cast<bvec_t*>(res.input(1).ptr());
  bvec_t* i2 = reinterpret_cast<bvec_t*>(res.input(2).ptr());
  bvec_t* o0 = reinterpret_cast<bvec_t*>(res.output(0).ptr());
  
  fill(o0,o0+3,0);
  o0[2] = 1;
  cout << "o0: " << endl;
  printBinary(o0[0]);
  printBinary(o0[1]);
  printBinary(o0[2]);

  for(int i=0; i<3; ++i){
    DMatrix& v = res.input(i);
    bvec_t* vv = reinterpret_cast<bvec_t*>(v.ptr());
    fill(vv,vv+v.size(),0);
  }
  
  res.spEvaluate(false);

  for(int i=0; i<3; ++i){
    DMatrix& v = res.input(i);
    bvec_t* vv = reinterpret_cast<bvec_t*>(v.ptr());
    cout << "input " << i << endl;
    for(int ii=0; ii<v.size(); ++ii){
      printBinary(vv[ii]);
    }
  }


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
