#include "../sx/sx_matrix.hpp"
#include "../expression_tools.hpp"
#include "../stl_vector_tools.hpp"
#include "../fx/sx_function.hpp"
#include <cassert>

// Get a pointer to the data for a ublas type
#define GET_DATA(x) ((x).data().begin())

using namespace std;
using namespace CasADi;

void test_evaluation(){

  // Declare variables
  SXMatrix x("x");
  SXMatrix y("y");
  SXMatrix z("z");


  SXMatrix f = 5+sin(x/2.2)+y*z;
  cout << "f = " << f << endl;

  SXMatrix xy;
  xy << x << y;
  SXMatrix jac_f = jacobian(f,xy);
  cout << "jac_f = " << jac_f << endl;
  
  

}

void test_making_smooth(){
  // Declare variables
  SXMatrix x("x");
  SXMatrix y("y");
  SXMatrix z("z");

  SXMatrix h = exp(sin(x) > y);
  std::cout << " h =  "<< h << std::endl;

  SXMatrix bvar, bexpr;
  makeSmooth(h,bvar, bexpr);
  std::cout << " h =  "<< h << std::endl;
  std::cout << "bvar = "<< bvar << std::endl;
  std::cout << "bexpr = "<< bexpr << std::endl;

}

int main(){

try{
    
  // This is how element access should work
  SXMatrix W(3,3);
  for(int i=0; i<9; ++i)
    W(i%3,i/3) = i+1;
  
  cout << W << endl;

  // Old way
  vector<int> Wi_old(3);
  Wi_old[0] = 1;
  Wi_old[1] = 1;
  Wi_old[2] = 2;
//  cout << W[Wi_old] << endl;
  
#ifdef WITH_CPP0X
  vector<int> Wi = {1,1,2};
  cout << W[Wi] << endl;

  //  cout << W[{1,2,3}] << endl; // still doesn't work -- bug or not yet implemented in gcc?
#endif
  
  SXMatrix W2 = W + W;
  cout << W << endl;
  cout << W2 << endl;
  SXMatrix W3 = binary(3,W,44);
  
  cout << W3 << endl;

//   return 0;

  
  test_evaluation();
  test_making_smooth();

  // Declare variables
  SX x("x");
  SX y("y");
  SX z("z");

  std::vector<int> uuu(2);
  uuu[0] = 1;
  uuu[1] = 22;
  
  std::vector<int> uuu2(2);
  uuu2[0] = 1000;
  uuu2[1] = 2000;

//   std::cout << (uuu + uuu2) << std::endl;


  std::vector<int> vvv(3);
  vvv[1] = 3;
  std::cout << toSXMatrix(vvv) << std::endl;
  
  // build an expression tree and print it
  SXMatrix n1 = x + 2;
  SXMatrix f = exp(y*z);
  std::cout << "f(x,y,z) = " << f << std::endl;

  // Make a vector of variables
  SXMatrix X(3,1);
  X(0) = x;
  X(1) = y;
  X(2) = z;
  std::cout << "X = " << X << std::endl;

  // differentiate an expression with respect to X
  SXMatrix grad_f = gradient(f,X);
  std::cout << "grad_f = " << grad_f << std::endl;

  // Vertical concatination
  SXMatrix m("m"),n("n"),o("o"),p("p"),q("q");
  SXMatrix mnopq;
  mnopq << m << n << o << p << q;
  std::cout << mnopq << std::endl; // appending

#ifdef HAVE_UBLAS
 // Evaluate at the point X0 using uBlas vectors
 ublas::vector<double> X0(3);
 X0[0] = 1.1;
 X0[1] = 6.1;
 X0[2] = 9.3;
 std::cout << "X0 = " << X0 << std::endl;
 
 // Evaluate the function numerically at X0
 ublas::vector<double> f0(1);

 f.evaluate(GET_DATA(X0));
 f.save(GET_DATA(f0));
 std::cout << "f(X0) = " << f0 << std::endl;

 // Also evaluate the gradient
 ublas::vector<double> grad_f0(3);
 grad_f.evaluate(GET_DATA(X0));
 grad_f.save(GET_DATA(grad_f0));
 std::cout << "grad_f(X0) = " << grad_f0 << std::endl;
#endif // HAVE_UBLAS

 // Evaluate the function symbolically
 SXMatrix Y(3,1);
 Y(0) = x;
 Y(1) = y;
 Y(2) = z;
 std::cout << "Y = " << Y << std::endl;

 // Matrices
 SXMatrix A(3,3);
 A(0,0) = x/x + x;
 A(2,2) = y*z;
 std::cout << "A = " << A << std::endl;
 std::cout << "X = " << X << std::endl;
 std::cout << "prod(A,X) = " << prod(A,X) << std::endl;
  
 std::cout << "A has " << nnz(A) << " non-zero elements" << std::endl;
 std::cout << A << std::endl;

 #ifdef HAVE_UBLAS
 // Evaluate to a ublas matrix
 A.init(Y);
 ublas::matrix<double> res(3,3);
 A.evaluate(GET_DATA(X0));
 A.save(GET_DATA(res));
#endif // HAVE_UBLAS

 // More matrices
  SXMatrix S(5,3);
  for(int i=0; i<5; ++i)
    for(int j=0; j<3; ++j)
      S(i,j) = 1;
  std::cout << S << std::endl;

  SXMatrix T(3,6);
  for(int i=0; i<3; ++i)
    for(int j=0; j<6; ++j)
      T(i,j) = 1;
  std::cout << T << std::endl;

  SXMatrix ST = prod(S,T); 
  std::cout << ST << std::endl;

  // Hessian
  SXMatrix g = 2*x*x-z*x+y;
  SXMatrix hess_g = hessian(g,Y);
  std::cout << hess_g << std::endl;
  std::cout << "the hessian has " << nnz(hess_g) << " (" << nnz_sym(hess_g) << ") non-zeros\n";

/*  SXMatrix grad_f = gradient(f,xx);
  cout << "grad_f(x,y,z) = " << grad_f << endl << endl;

  SXMatrix yy(3);
  yy[0] = x;
  yy[1] = y;
  yy[2] = z;

  SXMatrix gradf = f.gradient(yy);
  cout << "gradf(x,y,z) = " << gradf << endl << endl;
*/
// s.setArg(u);
 
//   SXMatrixMap x0;
//   x0[&x] = 1.2;
//   x0[&y] = 0.1;
//   x0[&z] = 11;

//  ublas::vector<double> res = f(x0);
//   cout << "f(x0) = " << f(x0) << endl << endl;

//  ublas::vector<double> vres = grad_f(x0);
/*  cout << "grad_f(x0) = " << grad_f(x0) << endl << endl;*/
  
//   SXMatrix a = 3.3;
//   TempSXMatrix b(a);
//   double c = b;
//   SXMatrix d = b;
// 
//   std::cout << a << std::endl;
// //  std::cout << b << std::endl;
//   std::cout << c << std::endl;
//   std::cout << d << std::endl;

  // SXMatrix of all ones
  SXMatrix aa = ones(2,9);

  // Identity matrix
  SXMatrix bb = eye(4);

/*  std::cout << aa.size() << std::endl;  */
  std::cout << aa << std::endl;  
  std::cout << bb << std::endl;
  std::cout << (aa.size1() == bb.size1() && aa.size2() == bb.size2()) << std::endl;

  // Check smoothness
  SXMatrix non_smooth = (x>2);
  std::cout << non_smooth << " is smooth? " << isSmooth(non_smooth) << std::endl;
  SXMatrix seven = 7;
  std::cout << x << " is smooth? " << isSmooth(x) << std::endl;
  std::cout << seven << " is smooth? " << isSmooth(seven) << std::endl;

  // Make a piececwise constant function in x
  SXMatrix tval = linspace(1,10,10);
  SXMatrix val = linspace(8,6,11);
  SXMatrix saw = pw_const(x,tval,val);
  std::cout << saw << std::endl;

  std::cout << "tval = " << tval << std::endl;
  std::cout << "val = " << val << std::endl;

  // Evaluate at some x
  SXMatrix t_test = linspace(0,10,40);
  
  SXFunction sawfcn(x,saw);
  
  SXMatrix saw_test = sawfcn.eval(t_test);

  std::cout << saw_test << std::endl;

  // Take the gradient
  SXMatrix grad_saw = gradient(saw,x);
  std::cout << grad_saw << std::endl;

  std::cout << "tan(atan(3)) = " << tan(atan(SXMatrix(3))) << std::endl;
  SXMatrix a2 = tan(atan(x));
  SXFunction a2fcn(x,a2);
  std::cout << "a2(3) = " << a2fcn.eval(3) << std::endl;

  // determinant
  SXMatrix Q("Q",2,2);
  std::cout << "Q = " << Q << ", det(Q) =  " << det(Q) << std::endl;

  SXMatrix R("R",3,3);
  std::cout << "R = " << R << ", det(R) =  " << det(R) << std::endl;

  // adjugate
  std::cout << "Q = " << Q << ", adj(Q) =  " << adj(Q) << std::endl;

#if 0
  // Test reading input
  SXMatrix inp1;
  std::cout << "Enter the first variable: ";
  std::cin >> inp1;
  if(inp1.isInteger())
    std::cout << "The variable is an integer" << std::endl;
  else if(inp1.isConstant())
    std::cout << "The variable is an real" << std::endl;
  else
    std::cout << "The variable is symbolic" << std::endl;

  SXMatrix inp2;
  std::cout << "Enter the second variable: ";
  std::cin >> inp2;

  if(inp1.isSymbolic() && inp2.isSymbolic() && inp1.getName().compare(inp2.getName())==0) inp2 = inp1;
  std::cout << "The sum is: " << inp1+inp2 << std::endl;
#endif

  SXMatrix u("u",4);
  std::cout << u << std::endl;

  SXMatrix u1_2 = u(1)*u(1);
  SXMatrix u3 = (u(0)*u(0) + u1_2 + u(2)*u(2)) - u1_2;


  SXMatrix u2 = prod(trans(u),u) + 3;
  std::cout << "u2 = " << u2 << std::endl;

  SXMatrix w,t;
  expand(u2,w,t);
  std::cout << "w = " << w << std::endl;
  std::cout << "t = " << t << std::endl;

  std::cout << "u3 = " << u3 << std::endl;
  simplify(u3);
  std::cout << "u3 (simplified) = " << u3 << std::endl;

  // Test substituting
  A = tan(x*x) - sqrt(x-y);
  std::cout << "A = " << A << std::endl;
  substitute(A,x,3*z);
  std::cout << "A (x->3*z) = " << A << std::endl;

  // Test making smooth
  SXMatrix nonsmooth = sin(3 - (x>=0));
  std::cout << "nonsmooth = " << nonsmooth << std::endl;
  SXMatrix bvar, bexpr;
  makeSmooth(nonsmooth,bvar,bexpr);

  std::cout << "smooth = " << nonsmooth << std::endl;
  std::cout << "binary variables = " << bvar << std::endl;
  std::cout << "binding expressions =  0 <= " << bvar - bexpr << " <= 1" << std::endl;

  // Test making banded matrices
//   SXMatrix AA(5,5);
//   AA(0,0) = 11;  AA(1,0) = 21;  AA(2,0) = 31;
//   AA(0,1) = 12;  AA(1,1) = 22;  AA(2,1) = 32;  AA(3,1) = 42;
//   AA(1,2) = 23;  AA(2,2) = 33;  AA(3,2) = 43;  AA(4,2) = 53;
//   AA(2,3) = 34;  AA(3,3) = 44;  AA(4,3) = 54;
//   AA(3,4) = 45;  AA(4,4) = 55;
//   std::cout << AA << std::endl;

  // Automatic diff
  SXMatrix G(2,2);
  G(0,0) = sin(x);
  G(1,1) = cos(x*x);
  G(0,1) = tan(x);
  std::cout << "G = " << G << std::endl;

  SXMatrix V(2,1);
  V(0) = exp(x);
  V(1) = x;
  std::cout << "V = " << V << std::endl;

  SXMatrix GV = prod(G,V);
  std::cout << "prod(G,V) = " << GV << std::endl;

  SXMatrix gradGV = gradient(GV,x);
  std::cout << "gradGV = " << gradGV << std::endl;

  // Integration using Gauss quadrature
  SXMatrix f1 = -x*x*exp(0.222*x*y);
  std::cout << "f1 = " << f1 << std::endl;
  SXMatrix f1_prim = gauss_quadrature(f1, x, 0, 1);
  std::cout << "f1_prim = " << f1_prim << std::endl;  

  // QR factorization using Modified Gram-Schmidt
  {
  double Adata_array[3][3] = {{2.3,  4.1, 9.2},{1.3,  0.1, 13.2},{-2.3, 0, -8.2}};
  vector<double> Adata(&Adata_array[0][0],&Adata_array[0][0]+9);
  SXMatrix A(Adata,3,3);
  
    
  // Factorize
  SXMatrix Q,R;
  qr(A,Q,R);
  cout << endl;
  cout << "QR decomposition of " << A << " is " << endl;
  cout << "Q = " << Q << endl;
  cout << "R = " << R << endl;
  cout << "prod(Q,R) = " << prod(Q,R) << endl;
  cout << "prod(trans(Q),Q) = " << prod(trans(Q),Q) << endl;
  cout << "prod(Q,trans(Q)) = " << prod(Q,trans(Q)) << endl;

  // Solve a linear system
  vector<double> bdata(3);
  bdata[0] = -4.4;
  bdata[1] = 1;
  bdata[2] = -22;
  SXMatrix b(bdata); // right-hand side

  cout << "prod(R\\trans(Q),x) = " <<  solve(R,prod(trans(Q),b)) << endl; // solve by QR factorization
  cout << "prod(inv(A),b) = " <<  prod(inv(A),b) << endl; // solve by calculating the inverse

  // Symbolic QR factorization
  SXMatrix Asym("A",10,3);

  // Factorize
  SXMatrix Qsym,Rsym;
  qr(Asym,Qsym,Rsym);
  cout << endl;
  cout << "QR decomposition of " << Asym << " is " << endl;
  cout << "Qsym = " << endl;
  SXFunction QsymFcn(SXMatrix(),Qsym);
  QsymFcn->printAlgorithm();
  cout << "Rsym = " << endl;
  SXFunction RsymFcn(SXMatrix(),Rsym);
  RsymFcn->printAlgorithm();

  cout << "spy(Qsym) = " << spy(Qsym) << endl;
  cout << "spy(Rsym) = " << spy(Rsym) << endl;

  // SXMatrix from blocks
  SXMatrix AAA("AAA"), BBB("BBB"), CCC("CCC"), DDD("DDD");
  SXMatrix TBT_data[2][2] = {{AAA,5},{CCC,DDD}};
  SXMatrix TBT = blockmatrix<2,2>(TBT_data);
  cout << "TBT = " << TBT << endl;

  // Calculate the inverse of a matrix in two different ways
  SXMatrix M("M",7,7);
  cout << "M = " << M << endl;
  
  // By Laplace's formula
  SXMatrix inv_M_laplace = inv(M);
  cout << "number of nodes, Laplace: " << numNodes(inv_M_laplace) << endl;

  // By QR factorization
  SXMatrix QM,RM;
  qr(M,QM,RM);
  SXMatrix inv_M_QR = solve(RM,trans(QM));
  cout << "number of nodes, QR: " << numNodes(inv_M_QR) << endl;
  
  }

  return 0;

} catch (const char * str){
  std::cerr << str << std::endl;
  return 1;
}


}
