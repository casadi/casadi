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

#include "sx_tools.hpp"
#include "../fx/sx_function_internal.hpp"
#include "binary_functions.hpp"

namespace CasADi{

void casadi_assert(bool cond, const std::string& msg="undefined error"){
  if(!cond){
    throw CasadiException(msg);
  }
}
  
SXMatrix gauss_quadrature(SXMatrix f, const SXMatrix &x, const SXMatrix &a, const SXMatrix &b, int order, const SXMatrix& w){
  casadi_assert(order == 5, "gauss_quadrature: order must be 5");
  casadi_assert(w.empty(),"gauss_quadrature: empty weights");

  // Change variables to [-1,1]
  if(!isEqual(a,-1) || !isEqual(b,1)){
    SXMatrix q1 = (b-a)/2;
    SXMatrix q2 = (b+a)/2;

    SXFunction fcn(x,f);

    return q1*gauss_quadrature(fcn.eval(q1*x+q2), x, -1, 1);
  }

  // Gauss points
  vector<double> xi;
  xi.push_back(-sqrt(5 + 2*sqrt(10.0/7))/3);
  xi.push_back(-sqrt(5 - 2*sqrt(10.0/7))/3);
  xi.push_back(0);
  xi.push_back(sqrt(5 - 2*sqrt(10.0/7))/3);
  xi.push_back(sqrt(5 + 2*sqrt(10.0/7))/3);

  // Gauss weights
  vector<double> wi;
  wi.push_back((322-13*sqrt(70.0))/900.0);
  wi.push_back((322+13*sqrt(70.0))/900.0);
  wi.push_back(128/225.0);
  wi.push_back((322+13*sqrt(70.0))/900.0);
  wi.push_back((322-13*sqrt(70.0))/900.0);
  
  // Evaluate at the Gauss points
  SXFunction fcn(x,f);
  vector<SX> f_val(5);
  for(int i=0; i<5; ++i)
    f_val[i] = fcn.eval(xi[i])[0];

  // Weighted sum
  SX sum;
  for(int i=0; i<5; ++i)
    sum += wi[i]*f_val[i];

  return sum;
}

SXMatrix linspace(const SXMatrix &a_, const SXMatrix &b_, int nsteps){
  if(!isScalar(a_) || !isScalar(b_)) throw "SXMatrix::linspace: a and b must be scalar";
  SX a = a_(0);
  SX b = b_(0);
  SXMatrix ret(nsteps,1);
  ret(0) = a;
  SX step = (b-a)/(nsteps-1);

  for(int i=1; i<nsteps-1; ++i)
    ret(i) = ret(i-1) + step;
  
  ret(nsteps-1) = b;
  return ret;
}

SX det(const SXMatrix& a){
  int n = a.size1();
  if(n != a.size2()) throw "det: matrix must be square";

  // Trivial return if scalar
  if(isScalar(a)) return a(0);

  // Return expression
  SX ret = 0;

  // We expand the matrix along the first column
  for(int i=0; i<n; ++i){

    // Sum up the cofactors
    ret += a(i,0)*cofactor(a,i,0);

  }
  return ret;
}

SX getMinor(const SXMatrix &x, int i, int j){
  int n = x.size1();
  if(n != x.size2()) throw "det: matrix must be square";

  // Trivial return if scalar
  if(n==1) return 1;

  // Remove row i and column j
  SXMatrix M(n-1,n-1);

   for(int i1=0; i1<n; ++i1)
       for(int j1=0; j1<n; ++j1){
	   if(i1 == i || j1 == j)
	      continue;

	    int i2 = (i1<i)?i1:i1-1;
	    int j2 = (j1<j)?j1:j1-1;
    
            M(i2,j2) = x(i1,j1);
       }
  return det(M);
}

SX cofactor(const SXMatrix &x, int i, int j){

    // Calculate the i,j minor
    SX minor_ij = getMinor(x,i,j);

    // Calculate the cofactor
    int sign_i = 1-2*((i+j) % 2);

    return sign_i * minor_ij;
}


SXMatrix adj(const SXMatrix& a){
  int n = a.size1();
  if(n != a.size2()) throw "adj: matrix must be square";

  // Cofactor matrix
  SXMatrix C(n,n);
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      C(i,j) = cofactor(a,i,j);
  
  return trans(C);
}

SXMatrix inv(const SXMatrix& a){
  // laplace formula
  return adj(a)/det(a);
}

SXMatrix pw_const(const SXMatrix &t, const SXMatrix &tval, const SXMatrix &val){
  // number of intervals
  int n = val.numel();

  if(!isScalar(t)) throw "SXMatrix::pw_const: t must be a scalar";
  if(tval.numel() != n-1 ) throw "SXMatrix::pw_const: dimensions do not match";

  SXMatrix ret = val(0);  
  for(int i=0; i<n-1; ++i){
    ret += (val(i+1)-val(i)) * (t>=tval(i));
  }

  return ret;
}

SXMatrix pw_lin(const SX &t, const SXMatrix &tval, const SXMatrix &val){
  // Number of points
  int N = tval.numel();
  casadi_assert(N>=2,"pw_lin: N>=2");

  // Gradient for each line segment
  SXMatrix g(N-1,1);
  for(int i=0; i<N-1; ++i){
    g(i) = (val(i+1)- val(i))/(tval(i+1)-tval(i));
  }

  // Line segments
  SXMatrix lseg(N-1,1);
  for(int i=0; i<N-1; ++i)
    lseg(i) = val(i) + g(i)*(t-tval(i)); 

  // interior time points
  SXMatrix tint;
  getSub(tint,tval,1,0,N-2); 

  // Return piecewise linear function
  return pw_const(t, tint, lseg);
}

SXMatrix if_else(const SXMatrix &cond, const SXMatrix &if_true, const SXMatrix &if_false){
  return if_false + (if_true-if_false)*cond;
}

SXMatrix ones(int n, int m){
  SXMatrix ret(n,m);
  for(int i=0; i<n; ++i)
    for(int j=0; j<m; ++j)
      ret(i,j) = 1;
  return ret;
}

SXMatrix zeros(int n, int m){
  return SXMatrix(n,m);
}

SXMatrix inf(int n, int m){
  SXMatrix ret(n,m);
  for(int i=0; i<n; ++i)
    for(int j=0; j<m; ++j)
      ret(i,j) = SX::inf;
  return ret;
}

SXMatrix eye(int n){
  SXMatrix ret(n,n);
  for(int i=0; i<n; ++i)
    ret(i,i) = 1;
  return ret;
}

SXMatrix heaviside(const SXMatrix& a){
  return (1+sign(a))/2;
}

SXMatrix sign(const SXMatrix& a){
  return ((a>0) + (a>=0))/2;
}


SXMatrix vec(const SXMatrix &expr){  
  // Quick return if empty or scalar
  if(expr.empty() || isScalar(expr)) return expr;
  SXMatrix ret(expr.numel(),1);
  int n = expr.size1();
  for(int i=0; i<expr.numel(); ++i)
    ret(i) = expr(i%n,i/n);
  return ret;
}

void getSub(SXMatrix &res, const SXMatrix &expr, int i, int j, int ni, int nj, int ki, int kj){
  casadi_assert(ki==1 && kj==1,"getSub: ki==1 && kj==1");
  casadi_assert(i+ni <= expr.size1() && j+nj <= expr.size2(),"getSub: i+ni <= expr.size1() && j+nj <= expr.size2()");
  res = SXMatrix(ni,nj);
  for(int r=0; r<ni; ++r)
    for(int c=0; c<nj; ++c)
      if(!expr(i+r,j+c)->isZero())
        res(r,c) = expr(i+r,j+c);
}

//void SXMatrix::set(const SXMatrix &expr, int i, int j, int ni, int nj, int ki, int kj){
void setSub(const SXMatrix &expr, SXMatrix &res, int i, int j){
  for(int r=0; r<expr.size1(); ++r)
    for(int c=0; c<expr.size2(); ++c)
      if(!expr(r,c)->isZero())
        res(i+r,j+c) = expr(r,c);
}

void getRow(SXMatrix &res, const SXMatrix &expr, int i, int ni, int ki){
  casadi_assert(i<expr.size1(),"getRow: i<expr.size1()");
  res = SXMatrix(ni,expr.size2());
  for(int ii=0; ii<ni; ++ii)
    for(int j=0; j<expr.size2(); ++j){
      SX temp = expr(i+ii,j);
      if(!temp->isZero()) res(ii,j) = temp;
      }
}

SXMatrix getRow( const SXMatrix &expr, int i, int ni, int ki){
  SXMatrix res(ni,expr.size2());
  getRow(res,expr,i,ni,ki);
  return res;
}

void getColumn(SXMatrix &res, const SXMatrix &expr, int j, int nj, int kj){
  casadi_assert(j<expr.size2(),"getColumn: j<expr.size2()");
  res = SXMatrix(expr.size1(),nj);
  for(int i=0; i<expr.size1(); ++i)
    for(int jj=0; jj<nj; ++jj){
      SX temp = expr(i,j+jj);
      if(!temp->isZero()) res(i,jj) = temp;
      }
}

SXMatrix getColumn( const SXMatrix &expr, int j, int nj, int kj){
  SXMatrix res(expr.size1(),nj);
  getColumn(res,expr,j,nj,kj);
  return res;
}

void setRow(const SXMatrix& expr, SXMatrix &res, int i, int ni, int ki){
  casadi_assert(i<res.size1(),"setRow: i<res.size1()");
  for(int j=0; j<res.size2(); ++j){
    if(!expr(0,j)->isZero())
      res(i,j) = expr(0,j);
    }
}

void setColumn(const SXMatrix& expr, SXMatrix &res, int j, int nj, int kj){
  casadi_assert(j<res.size2(),"setColumn: j<res.size2()");
  for(int i=0; i<res.size1(); ++i){
    if(!expr(i,0)->isZero())
      res(i,j) = expr(i,0);
    }
}

bool contains(const SXMatrix &list, const SX &e) {
  for (int i=0;i<nnz(list);i++) {
    if (list(i)->isEqual(e)) return true;
  }
  return false;
}


void simplify(SXMatrix &ex){
  // simplify all non-zero elements
  for(int el=0; el<ex.size(); ++el)
    simplify(ex[el]);
}

void compress(SXMatrix &ex, int level){

  throw "SXMatrix::compress: Not implemented";

  if(level>0)
    compress(ex,level-1);
}

void substitute(SXMatrix &ex, const SXMatrix &var, const SXMatrix &expr){
  if(var.empty()) return; // quick return if empty
  if(!isSymbolic(var)) throw "SXMatrix::substitute: the variable is not symbolic";
  if(var.size1() != expr.size1() || var.size2() != expr.size2()) throw "SXMatrix::substitute: the dimensions do not match";

  // evaluate with var == expr
  SXFunction fcn(var,ex);
  ex = fcn.eval(expr);
}


#if 0
void replaceDerivatives(SXMatrix &ex, const SXMatrix &var, const SXMatrix &dvar){
  // Initialize with an empty expression
  SXFunction fcn(ex);

  // Map from the var-node to the new der-node
  std::map<int, SX> dermap;
  for(int i=0; i<var.size(); ++i)
    dermap[fcn.treemap[var[i].get()]] = dvar[i];

  // Replacement map
  std::map<int, SX> replace;

  // Go through all nodes and check if any node is a derivative
  for(int i=0; i<fcn.algorithm.size(); ++i){
        if(fcn.algorithm[i].op == DER_NODE){

          // find the corresponding derivative
          std::map<int, SX>::iterator r = dermap.find(fcn.algorithm[i].ch0);
          casadi_assert(r != dermap.end());

          replace[i] = r->second;
        }
  }
  SXMatrix res;
  SXMatrix repres;
  fcn.eval_symbolic(SXMatrix(),res,replace,repres);
  ex = res;

  casadi_assert(0);

}
#endif


void makeSmooth(SXMatrix &ex, SXMatrix &bvar, SXMatrix &bexpr){
  // Initialize
  SXFunction fcn(SXMatrix(),ex);

  casadi_assert(bexpr.empty());

  // Nodes to be replaced
  std::map<int,SX> replace;

  // Go through all nodes and check if any node is non-smooth
  for(int i=0; i<fcn->algorithm.size(); ++i){

      // Check if we have a step node
      if(fcn->algorithm[i].op == STEP_NODE){

        // Get the index of the child
        int ch0 = fcn->algorithm[i].ch[0];

        // Binary variable corresponding to the the switch
        SXMatrix sw;

#if 0 
        // Find out if the switch has already been added
        for(int j=0; j<bexpr.size(); ++j)
          if(bexpr[j].isEqual(algorithm[i]->child0)){
            sw = bvar[j];
            break;
          }
#endif

        if(sw.empty()){ // the switch has not yet been added
          // Get an approriate name of the switch
          std::stringstream name;
          name << "sw_" << bvar.size1();
          sw = SX(name.str());
  
          // Add to list of switches
          bvar << sw;
//        bexpr << algorithm[i]->child0;
        }
        
        // Add to the substition map
        replace[i] = sw[0];
      }
  }
  SXMatrix res;
  fcn->eval(SXMatrix(),res,replace,bexpr);

  for(int i=0; i<bexpr.size(); ++i)
    bexpr[i] = bexpr[i]->dependent(0);

  ex = res;

#if 0
  // Make sure that the binding expression is smooth
  bexpr.init(SXMatrix());
  SXMatrix b;
  bexpr.eval_symbolic(SXMatrix(),b,replace,bexpr);
  bexpr = b;
#endif
}

void qr(const SXMatrix& A, SXMatrix& Q, SXMatrix &R){
  // The following algorithm is taken from J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.)
  int m = A.size1();
  int n = A.size2();
  casadi_assert(m>=n);

  // Transpose of A
  SXMatrix AT = trans(A);

  // Transposes of the output matrices
  SXMatrix QT, RT;

  // compute Q and R column by column
  for(int i=0; i<n; ++i){
    // Initialize qi to be the i-th column of A
    SXMatrix ai = getRow(AT,i);
    SXMatrix qi = ai;
    // The i-th column of R
    SXMatrix ri(1,n);
  
    // subtract the projection of of qi in the previous directions from ai
    for(int j=0; j<i; ++j){
      // Get the j-th column of Q
      SXMatrix qj = getRow(QT,j);

      ri(0,j) = prod(qj,trans(qi))[0]; // Modified Gram-Schmidt
      // ri[j] = inner_prod(qj,ai); // Classical Gram-Schmidt

      // Remove projection in direction j
      qi -= SX(ri(0,j)) * qj;
    }

    // Normalize qi
    ri(0,i) = norm_2(trans(qi))[0];
    qi /= SX(ri(0,i));

    // Update RT and QT
    QT << qi;
    RT << ri;

  }

  // Save to output
  Q = trans(QT);
  R = trans(RT);
}

SXMatrix spy(const SXMatrix& A){
  SXMatrix s(A.size1(),A.size2());
  for(int i=0; i<A.size1(); ++i)
    for(int j=0; j<A.size2(); ++j)
      if(!A(i,j)->isZero())
        s(i,j) = 1;
  return s;
}

bool isTril(const SXMatrix &A){
  // loop over rows
  for(int i=0; i<A.size1(); ++i){
    if(A.rowind(i) != A.rowind(i+1)){ // if there are any elements of the row
      // check column of the right-most element of the row
      int col = A.col(A.rowind(i+1)-1);

      // not lower triangular if col>i
      if(col>i) return false;
    }
  }
  // all rows ok
  return true;
}

bool isTriu(const SXMatrix &A){
  // loop over rows
  for(int i=0; i<A.size1(); ++i){
    if(A.rowind(i) != A.rowind(i+1)){ // if there are any elements of the row
      // check column of the left-most element of the row
      int col = A.col(A.rowind(i));

      // not lower triangular if col>i
      if(col<i) return false;
    }
  }
  // all rows ok
  return true;
}

SXMatrix solve(const SXMatrix& A, const SXMatrix& b){
  // check dimensions
  casadi_assert(A.size1() == b.size1());
  casadi_assert(A.size1() == A.size2()); // make sure that A is square  

  if(isTril(A)){
    // forward substiution if lower triangular
    SXMatrix x = b;
    for(int i=0; i<A.size1(); ++i){ // loop over rows
      for(int k=0; k<b.size2(); ++k){ // for every right hand side
        for(int j=0; j<i; ++j){
          x(i,k) -= A(i,j)*x(j,k);
        }
        x(i,k) /= A(i,i);
      }
    }
    return x;
  } else if(isTriu(A)){
    // backward substiution if upper triangular
    SXMatrix x = b;
    for(int i=A.size1()-1; i>=0; --i){ // loop over rows from the back
      for(int k=0; k<b.size2(); ++k){ // for every right hand side
        for(int j=A.size1()-1; j>i; --j){
          x(i,k) -= A(i,j)*x(j,k);
        }
        x(i,k) /= A(i,i);
      }
    }
    return x;
  } else {
    // Make a QR factorization
    SXMatrix Q,R;
    qr(A,Q,R);

    // Solve the factorized system
    return solve(R,prod(trans(Q),b));
  }
}

int numNodes(const SXMatrix& A){
  // Create a function
  SXFunction fcn(SXMatrix(),A);

  return fcn->tree.size();
}

bool dependsOn(const SXMatrix& ex, const SXMatrix &arg){
  SXMatrix g = gradient(vec(ex),arg);
  return !isZero(g);
}

bool isScalar(const SXMatrix& ex){
  return ex.size1()==1 && ex.size2()==1;
}

bool isVector(const SXMatrix& ex){
  return ex.size2()==1;
}

bool isSmooth(const SXMatrix& ex){
 // Make a function
 SXFunction temp(SXMatrix(),ex);
  
 // Run the function on the temporary variable
 return temp->isSmooth();
}

bool isInteger(const SXMatrix& ex){
  if(ex.empty())
    throw "isConstant(): SXMatrix is empty";
  else if(isScalar(ex))
    return ex(0)->isInteger();
  else
  {
    for(int k=0; k<ex.size(); ++k) // loop over non-zero elements
      if(!ex[k]->isInteger()) // if an element is not integer
        return false;
  }
  return true;
}

bool isConstant(const SXMatrix& ex){
  if(ex.empty())
    throw "SXMatrix::isConstant(): SXMatrix is empty";
  else
  {
    for(int k=0; k<ex.size(); ++k) // loop over non-zero elements
      if(!ex[k]->isConstant()) // if an element is not constant
        return false;
  }
  return true;
}

bool isDense(const SXMatrix& ex){
  return ex.size() == ex.numel();
}

bool isEmpty(const SXMatrix& ex){
  return ex.empty();
}

bool isSymbolic(const SXMatrix& ex){
  if(!isDense(ex)) return false;

  for(int k=0; k<ex.size(); ++k) // loop over non-zero elements
    if(!ex[k]->isSymbolic()) // if an element is not symbolic
      return false;
  
  return true;
}

bool isEqual(const SXMatrix& ex1,const SXMatrix &ex2){
  if ((nnz(ex1)!=0 || nnz(ex2)!=0) && (ex1.size1()!=ex2.size1() || ex1.size2()!=ex2.size2())) return false;
  SXMatrix difference = ex1 - ex2;  
  return isZero(difference);
}


SXMatrix gradient(const SXMatrix& ex, const SXMatrix &arg) {
  return trans(jacobian(ex,arg));
}
  
SXMatrix jacobian(const SXMatrix& ex, const SXMatrix &arg) {
  SXFunction temp(arg,ex); // make a runtime
  return temp->jac();
}

void hessian(const SXMatrix& ex, const SXMatrix &arg, SXMatrix &H, SXMatrix &g) {
  // this algorithm is _NOT_ linear time (but very easy to implement).. Change to higher order AD!
  g = gradient(ex,arg);  
  H = gradient(g,arg);
}

SXMatrix hessian(const SXMatrix& ex, const SXMatrix &arg) {
  SXMatrix H,g;
  hessian(ex,arg,H,g);
  return H;
}

double getValue(const SXMatrix& ex, int i, int j) {
  casadi_assert(i<ex.size1() && j<ex.size2());
  return ex(i,j)->getValue();
}

int getIntValue(const SXMatrix& ex, int i, int j) {
  casadi_assert(i<ex.size1() && j<ex.size2());
  return ex(i,j)->getIntValue();
}

void getValue(const SXMatrix& ex, double *res) {
  for(int i=0; i<ex.numel(); ++i)
    res[i] = ex(i)->getValue();
}

void getIntValue(const SXMatrix& ex, int *res) {
  for(int i=0; i<ex.numel(); ++i)
    res[i] = ex(i)->getIntValue();
}

const string& getName(const SXMatrix& ex) {
  if(!isScalar(ex)) throw "SXMatrix::getName: the expression must be scalar";
  return ex(0)->getName();
}

bool isZero(const SXMatrix& ex) {  

  // loop over (potentially) non-zero elements
  for(int el=0; el<ex.size(); ++el)
    if(!ex[el]->isZero())
      return false;
  
  return true;
}

int nnz(const SXMatrix& ex) {
  return ex.size();
}

int nnz_sym(const SXMatrix& ex) {
  int nz = 0; // number of non-zeros  
  for(int row=0; row<ex.size1(); ++row)
  {
    // Loop over the elements in the row
    for(int el=ex.rowind(row); el<ex.rowind(row+1); ++el){ // loop over the non-zero elements
      if(ex.col(el) > row) break; // break inner loop (only lower triangular part is used)
      nz++;
    }
  }
  return nz;
}

istream& operator>>(istream &stream, SXMatrix &expr){
    // try to read a double
    double realvalue;
    stream >> realvalue;
    if(!stream.fail()){
      // Check if the value is an integer
      int intvalue = (int)realvalue;
      if(intvalue == realvalue) expr = intvalue;
      else                      expr = realvalue;
      return stream;
    }
    stream.clear();

    // Here try to read a vector or matrix

    // else create a symbolic variable    
    string str;
    stream >> str;
    expr = SXMatrix(str);
    return stream;
}



// SXMatrix vertcat(const std::vector<SX>& comp){
//   return toSXMatrix(comp);
// }
// 
// SXMatrix vertcat(const std::vector<SXMatrix>& comp){
//   return toSXMatrix(comp);
// }
SXMatrix operator<=(const SXMatrix &a, const SXMatrix &b){
  if(!a.scalar() || !b.scalar()) throw "conditional operators only defined for scalars";
  return a.getElement() <= b.getElement();
}

SXMatrix operator>=(const SXMatrix &a, const SXMatrix &b){
  if(!a.scalar() || !b.scalar()) throw "conditional operators only defined for scalars";
  return a.getElement() >= b.getElement();
}

SXMatrix operator<(const SXMatrix &a, const SXMatrix &b){
  if(!a.scalar() || !b.scalar()) throw "conditional operators only defined for scalars";
  return a.getElement() < b.getElement();
}

SXMatrix operator>(const SXMatrix &a, const SXMatrix &b){
  if(!a.scalar() || !b.scalar()) throw "conditional operators only defined for scalars";
  return a.getElement() > b.getElement();
}

SXMatrix operator&&(const SXMatrix &a, const SXMatrix &b){
  if(!a.scalar() || !b.scalar()) throw "conditional operators only defined for scalars";
  return a.getElement() && b.getElement();
}

SXMatrix operator||(const SXMatrix &a, const SXMatrix &b){
  if(!a.scalar() || !b.scalar()) throw "conditional operators only defined for scalars";
  return a.getElement() || b.getElement();
}

SXMatrix operator==(const SXMatrix &a, const SXMatrix &b){
  if(!a.scalar() || !b.scalar()) throw "conditional operators only defined for scalars";
  return a.getElement() == b.getElement();
}

SXMatrix operator!=(const SXMatrix &a, const SXMatrix &b){
  if(!a.scalar() || !b.scalar()) throw "conditional operators only defined for scalars";
  return a.getElement() != b.getElement();
}

SXMatrix operator!(const SXMatrix &a){
  if(!a.scalar()) throw "conditional operators only defined for scalars";
  return !a.getElement();
}

SXMatrix norm_2(const SXMatrix& x){
  return sqrt(inner_prod(x,x));
}

SXMatrix inner_prod(const SXMatrix &x, const SXMatrix &y){
  if(!x.vector() || !y.vector()) throw CasadiException("inner_prod: arguments must be vectors");
  return prod(trans(x),y);
}

SXMatrix outer_prod(const SXMatrix &x, const SXMatrix &y){
  if(!x.vector() || !y.vector()) throw CasadiException("outer_prod: arguments must be vectors");
  return prod(x,trans(y));  
}
SXMatrix vertcat(const std::vector<SXMatrix> &v){
  SXMatrix ret;
  for(int i=0; i<v.size(); ++i)
    ret << v[i];
  return ret;
}

SXMatrix horzcat(const std::vector<SXMatrix> &v){
  SXMatrix ret;
  for(int i=0; i<v.size(); ++i)
    ret << trans(v[i]);
  return trans(ret);  
}

SXMatrix vertcat(const SXMatrix &x, const SXMatrix &y){
  vector<SXMatrix> xy(2);
  xy[0]=x;
  xy[1]=y;
  return vertcat(xy);
}

SXMatrix horzcat(const SXMatrix &x, const SXMatrix &y){
  vector<SXMatrix> xy(2);
  xy[0]=x;
  xy[1]=y;
  return horzcat(xy);
}

SXMatrix& operator<<(SXMatrix& expr, const SXMatrix& add){
  // Quick return if we are adding an empty expression
  if(add.empty()) return expr;

  // Likewise if expr is empty
  if(expr.empty()) return expr=add;

  // Check dimensions
  if(expr.size2() != add.size2()) throw "operator<<: dimensions do not match";

  // Resize the expression
  int oldn = expr.size1();
  int n    = expr.size1() + add.size1();  
  int m    = expr.size2();
  expr.resize(n,m);

  // Copy the lower expression to the end
  setSub(add, expr, oldn, 0);
  return expr;
}

  
std::vector< std::vector< std::vector< SX> > > create_symbolic(const std::string& name, int n, int m, int p);



void expand(const SXMatrix& ex2, SXMatrix &ww, SXMatrix& tt){
  casadi_assert(ex2.scalar());
  SX ex = ex2.getElement();
  
  // Terms, weights and indices of the nodes that are already expanded
  std::vector<std::vector<SXNode*> > terms;
  std::vector<std::vector<double> > weights;
  std::map<SXNode*,int> indices;

  // Stack of nodes that are not yet expanded
  std::stack<SXNode*> to_be_expanded;
  to_be_expanded.push(ex.get());

  while(!to_be_expanded.empty()){ // as long as there are nodes to be expanded

    // Check if the last element on the stack is already expanded
   if (indices.find(to_be_expanded.top()) != indices.end()){
      // Remove from stack
      to_be_expanded.pop();
      continue;
    }

    // Weights and terms
    std::vector<double> w; // weights
    std::vector<SXNode*> f; // terms

    if(to_be_expanded.top()->isConstant()){ // constant nodes are seen as multiples of one
      w.push_back(to_be_expanded.top()->getValue());
      f.push_back(SX::one.get());
    } else if(to_be_expanded.top()->isSymbolic()){ // symbolic nodes have weight one and itself as factor
      w.push_back(1);
      f.push_back(to_be_expanded.top());
    } else { // binary node

        casadi_assert(to_be_expanded.top()->isBinary()); // make sure that the node is binary

        // Check if addition, subtracton or multiplication
        BinarySXNode*     binnode = (BinarySXNode*)to_be_expanded.top();
        // If we have a binary node that we can factorize
        if(binnode->op == ADD_NODE || binnode->op == SUB_NODE || (binnode->op == MUL_NODE  && (binnode->child[0]->isConstant() || binnode->child[1]->isConstant()))){
          // Make sure that both children are factorized, if not - add to stack
          if (indices.find(binnode->child[0].get()) == indices.end()){
            to_be_expanded.push(binnode->child[0].get());
            continue;
          }
          if (indices.find(binnode->child[1].get()) == indices.end()){
             to_be_expanded.push(binnode->child[1].get());
             continue;
          }

          // Get indices of children
          int ind1 = indices[binnode->child[0].get()];
          int ind2 = indices[binnode->child[1].get()];
  
          // If multiplication
          if(binnode->op == MUL_NODE){
            double fac;
            if(binnode->child[0]->isConstant()){ // Multiplication where the first factor is a constant
              fac = binnode->child[0]->getValue();
              f = terms[ind2];
              w = weights[ind2];
            } else { // Multiplication where the second factor is a constant
              fac = binnode->child[1]->getValue();
              f = terms[ind1];
              w = weights[ind1];
            }
            for(int i=0; i<w.size(); ++i) w[i] *= fac;

          } else { // if addition or subtraction
            if(binnode->op == ADD_NODE){          // Addition: join both sums
              f = terms[ind1];      f.insert(f.end(), terms[ind2].begin(), terms[ind2].end());
              w = weights[ind1];    w.insert(w.end(), weights[ind2].begin(), weights[ind2].end());
            } else {      // Subtraction: join both sums with negative weights for second term
              f = terms[ind1];      f.insert(f.end(), terms[ind2].begin(), terms[ind2].end());
              w = weights[ind1];
              w.reserve(f.size());
              for(int i=0; i<weights[ind2].size(); ++i) w.push_back(-weights[ind2][i]);
            }
          // Eliminate multiple elements
          std::vector<double> w_new; w_new.reserve(w.size());   // weights
          std::vector<SXNode*> f_new;  f_new.reserve(f.size());   // terms
          std::map<SXNode*,int> f_ind; // index in f_new

          for(int i=0; i<w.size(); i++){
            // Try to locate the node
            std::map<SXNode*,int>::iterator it = f_ind.find(f[i]);
            if(it == f_ind.end()){ // if the term wasn't found
              w_new.push_back(w[i]);
              f_new.push_back(f[i]);
              f_ind[f[i]] = f_new.size()-1;
            } else { // if the term already exists
              w_new[it->second] += w[i]; // just add the weight
            }
          }
          w = w_new;
          f = f_new;
        }
      } else { // if we have a binary node that we cannot factorize
        // By default, 
        w.push_back(1);
        f.push_back(binnode);

      }
    }

    // Save factorization of the node
    weights.push_back(w);
    terms.push_back(f);
    indices[to_be_expanded.top()] = terms.size()-1;

    // Remove node from stack
    to_be_expanded.pop();
  }

  // Save expansion to output
  int thisind = indices[ex.get()];
  ww = toSXMatrix(weights[thisind]);

  vector<SX> termsv(terms[thisind].size());
  for(int i=0; i<termsv.size(); ++i)
    termsv[i] = SX(terms[thisind][i]);
  tt = toSXMatrix(termsv);
}

void simplify(SX& ex){
  // Start by expanding the node to a weighted sum
  SXMatrix terms, weights;
  expand(ex,weights,terms);

  // Make a scalar product to get the simplified expression
  SXMatrix s = prod(trans(weights),terms);
  ex = s.getElement();
}

SXMatrix trans(const SXMatrix& x){
  // quick return if empty or scalar
  if(x.empty() || x.scalar()) return x;

  // We do matrix transpose by the (linear time) "bucket sort" algorithm
  vector<vector<int> > buckets(x.size2()); // one bucket for each column
  
  // Create a vector with the rows for each non-zero element
  vector<int> row(x.size());
  
  // Loop over the rows of the original matrix
  for(int i=0; i<x.size1(); ++i)
  {
    // Loop over the elements in the row
    for(int el=x.rowind(i); el<x.rowind(i+1); ++el){ // loop over the non-zero elements
      int j=x.col(el);  // column
      
     // put the element into the right bucket
     buckets[j].push_back(el);
     
     // save the row index
     row[el] = i;
    }
  }

  SXMatrix ret(x.size2(),x.size1()); // create the return matrix

  // reserve space (to make the calculations quicker)
  ret.reserve(x.capacity());
  ret.col().reserve(x.col().size());
  ret.rowind().reserve(x.rowind().size());

  for(int j=0; j<x.size2(); ++j)   // loop over the columns
    for(int r=0; r<buckets[j].size(); ++r){ // loop over the bucket content
     int el =  buckets[j][r]; // the index of the non-zero element
     int i = row[el]; // the row of the element
     ret.getElementRef(j,i) = x[el]; // add the element
    }

    return ret;
}

SXMatrix prod(const SXMatrix &x, const SXMatrix &y){
  if(x.size2() != y.size1()) throw CasadiException("prod: dimension mismatch");

  SXMatrix ret(x.size1(),y.size2());
  SXMatrix b = trans(y); // take the transpose of the second matrix: linear time operation

  for(int i=0; i<x.size1(); ++i) // loop over the row of the resulting matrix)
    for(int j=0; j<b.size1(); ++j){ // loop over the column of the resulting matrix
      int el1 = x.rowind(i);
      int el2 = b.rowind(j);
      while(el1 < x.rowind(i+1) && el2 < b.rowind(j+1)){ // loop over non-zero elements
        int j1 = x.col(el1);
        int i2 = b.col(el2);      
        if(j1==i2){
          SX temp = x[el1++] * b[el2++];
          if(!temp->isZero())
            ret(i,j) += temp;
        } else if(j1<i2) {
          el1++;
        } else {
          el2++;
        }
      }
    }
return ret;
}

void fill(SXMatrix& mat, const SX& val){
  if(val->isZero())    mat.makeEmpty(mat.size1(),mat.size2());
  else                 mat.makeDense(mat.size1(),mat.size2(),val);
}

SXMatrix binary(int op, const SXMatrix &x, const SXMatrix &y){
  SXMatrix r;
  if(x.scalar())
    if(y.scalar())
      return sfcn[op](x(0),y(0));
    else
      return scalar_matrix(op,x(0),y);
  else if(y.scalar())
    return matrix_scalar(op,x,y(0));
  else
    return matrix_matrix(op,x,y);
}

SXMatrix unary(int op, const SXMatrix &x){
  if(x.scalar()){
    return sfcn[op](x(0),SX::nan);
  } else {
    SXMatrix r(x.size1(),x.size2());
    SX y; // dummy argument
    fill(r,sfcn[op](0,y));
    for(int i=0; i<r.size1(); ++i){ // loop over rows
      for(int el=x.rowind(i); el<x.rowind(i+1); ++el){
        int j = x.col(el);
        r(i,j) = sfcn[op](x[el],y);
      }
    }  
    return r;
  }
}

SXMatrix scalar_matrix(int op, const SX &x, const SXMatrix &y){
  SXMatrix r(y.size1(),y.size2());
  fill(r,sfcn[op](x,0));
  for(int i=0; i<r.size1(); ++i){ // loop over rows
    for(int el=y.rowind(i); el<y.rowind(i+1); ++el){
      int j = y.col(el);
      r(i,j) = sfcn[op](x,y[el]);
    }
  }
  return r;
}

SXMatrix matrix_scalar(int op, const SXMatrix &x, const SX &y){
  SXMatrix r(x.size1(),x.size2());
  fill(r,sfcn[op](0,y));
  for(int i=0; i<r.size1(); ++i){ // loop over rows
    for(int el=x.rowind(i); el<x.rowind(i+1); ++el){
      int j = x.col(el);
      r(i,j) = sfcn[op](x[el],y);
    }
  }
  return r;
}

SXMatrix matrix_matrix(int op, const SXMatrix &x, const SXMatrix &y){
if(x.size1() != y.size1() || x.size2() != y.size2()) throw CasadiException("matrix_matrix: dimension mismatch");
  SXMatrix r(x.size1(),x.size2());
  fill(r,sfcn[op](0,0));
  for(int i=0; i<r.size1(); ++i){ // loop over rows
    int el1 = x.rowind(i);
    int el2 = y.rowind(i);
    int k1 = x.rowind(i+1);
    int k2 = y.rowind(i+1);
    while(el1 < k1 || el2 < k2){
      int j1 = (el1 < k1) ? x.col(el1) : r.numel() ;
      int j2 = (el2 < k2) ? y.col(el2) : r.numel() ;
      
      if(j1==j2)
        r(i,j1) = sfcn[op](x[el1++],y[el2++]); 
      else if(j1>j2)
        r(i,j2) = sfcn[op](0,y[el2++]);
      else
        r(i,j1) = sfcn[op](x[el1++],0);
      }
    }
  return r;
}

void make_symbolic(SX& v, const std::string& name){
  v = SX(name);
}

vector<SX> create_symbolic(const string& name, int n){
  vector<SX> ret(n);
  make_symbolic(ret,name);
  return ret;
}

vector< vector<SX> > create_symbolic(const std::string& name, int n, int m){
  vector< vector<SX> > ret(n,vector<SX>(m));
  make_symbolic(ret,name);
  return ret;
}

vector< vector< vector<SX> > > create_symbolic(const std::string& name, int n, int m, int p){
  vector< vector< vector<SX> > > ret(n,vector< vector<SX> >(m, vector<SX>(p)));
  make_symbolic(ret,name);
  return ret;
}


} // namespace CasADi


namespace std{
using namespace CasADi;


SXMatrix sin(const SXMatrix& x){
  return unary(SIN_NODE,x);
}

SXMatrix cos(const SXMatrix& x){
  return unary(COS_NODE,x);
}

SXMatrix tan(const SXMatrix& x){
  return unary(TAN_NODE,x);
}

SXMatrix atan(const SXMatrix& x){
  return unary(ATAN_NODE,x);
}

SXMatrix asin(const SXMatrix& x){
  return unary(ASIN_NODE,x);
}

SXMatrix acos(const SXMatrix& x){
  return unary(ACOS_NODE,x);
}

SXMatrix exp(const SXMatrix& x){
  return unary(EXP_NODE,x);
}

SXMatrix log(const SXMatrix& x){
  return unary(LOG_NODE,x);
}

SXMatrix pow(const SXMatrix& x, const SXMatrix& n){
  return binary(POW_NODE,x,n);
}

SXMatrix sqrt(const SXMatrix& x){
  return unary(SQRT_NODE,x);
}

SXMatrix fmin(const SXMatrix& x, const SXMatrix& y){
  return binary(FMIN_NODE,x,y);
}

SXMatrix fmax(const SXMatrix& x, const SXMatrix& y){
  return binary(FMAX_NODE,x,y);
}

SXMatrix floor(const SXMatrix& x){
  return unary(FLOOR_NODE,x);
}

SXMatrix ceil(const SXMatrix& x){
  return unary(CEIL_NODE,x);
}

SXMatrix erf(const SXMatrix& x){
  return unary(ERF_NODE,x);
}


} // namespace std




