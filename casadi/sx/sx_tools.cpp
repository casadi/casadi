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
#include "../casadi_math.hpp"
#include "../matrix/matrix_tools.hpp"
using namespace std;

namespace CasADi{
  
Matrix<SX> gauss_quadrature(Matrix<SX> f, const Matrix<SX> &x, const Matrix<SX> &a, const Matrix<SX> &b, int order, const Matrix<SX>& w){
  casadi_assert_message(order == 5, "gauss_quadrature: order must be 5");
  casadi_assert_message(w.empty(),"gauss_quadrature: empty weights");

  // Change variables to [-1,1]
  if(!a(0)->isEqual(-1) || !b(0)->isEqual(1)){
    Matrix<SX> q1 = (b-a)/2;
    Matrix<SX> q2 = (b+a)/2;

    SXFunction fcn(x,f);

    return q1*gauss_quadrature(fcn.eval(q1*x+q2), x, -1, 1);
  }

  // Gauss points
  vector<double> xi;
  xi.push_back(-std::sqrt(5 + 2*std::sqrt(10.0/7))/3);
  xi.push_back(-std::sqrt(5 - 2*std::sqrt(10.0/7))/3);
  xi.push_back(0);
  xi.push_back(std::sqrt(5 - 2*std::sqrt(10.0/7))/3);
  xi.push_back(std::sqrt(5 + 2*std::sqrt(10.0/7))/3);

  // Gauss weights
  vector<double> wi;
  wi.push_back((322-13*std::sqrt(70.0))/900.0);
  wi.push_back((322+13*std::sqrt(70.0))/900.0);
  wi.push_back(128/225.0);
  wi.push_back((322+13*std::sqrt(70.0))/900.0);
  wi.push_back((322-13*std::sqrt(70.0))/900.0);
  
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

Matrix<SX> pw_const(const Matrix<SX> &t, const Matrix<SX> &tval, const Matrix<SX> &val){
  // number of intervals
  int n = val.numel();

  casadi_assert_message(isScalar(t),"t must be a scalar");
  casadi_assert_message(tval.numel() == n-1, "dimensions do not match");

  Matrix<SX> ret = val(0);  
  for(int i=0; i<n-1; ++i){
    ret += (val(i+1)-val(i)) * (t>=tval(i));
  }

  return ret;
}

Matrix<SX> pw_lin(const SX &t, const Matrix<SX> &tval, const Matrix<SX> &val){
  // Number of points
  int N = tval.numel();
  casadi_assert_message(N>=2,"pw_lin: N>=2");

  // Gradient for each line segment
  Matrix<SX> g(N-1,1);
  for(int i=0; i<N-1; ++i){
    g(i) = (val(i+1)- val(i))/(tval(i+1)-tval(i));
  }

  // Line segments
  Matrix<SX> lseg(N-1,1);
  for(int i=0; i<N-1; ++i)
    lseg(i) = val(i) + g(i)*(t-tval(i)); 

  // interior time points
  Matrix<SX> tint = tval(range(N-2),0);

  // Return piecewise linear function
  return pw_const(t, tint, lseg);
}

Matrix<SX> if_else(const Matrix<SX> &cond, const Matrix<SX> &if_true, const Matrix<SX> &if_false){
  return if_false + (if_true-if_false)*cond;
}

Matrix<SX> onesSX(int n, int m){
  return ones<SX>(n,m);
}

Matrix<SX> zerosSX(int n, int m){
  return zeros<SX>(n,m);
}

Matrix<SX> infSX(int n, int m){
  return Matrix<SX>(n,m,casadi_limits<SX>::inf);
}

Matrix<SX> eyeSX(int n){
  return eye<SX>(n);
}

Matrix<SX> heaviside(const Matrix<SX>& a){
  return (1+sign(a))/2;
}

Matrix<SX> sign(const Matrix<SX>& a){
  return ((a>0) + (a>=0))/2;
}

bool contains(const Matrix<SX> &list, const SX &e) {
  for (int i=0;i<nnz(list);i++) {
    if (list(i)->isEqual(e)) return true;
  }
  return false;
}


void simplify(Matrix<SX> &ex){
  // simplify all non-zero elements
  for(int el=0; el<ex.size(); ++el)
    simplify(ex[el]);
}

void compress(Matrix<SX> &ex, int level){

  throw CasadiException("Matrix<SX>::compress: Not implemented");

  if(level>0)
    compress(ex,level-1);
}

Matrix<SX> substitute(const Matrix<SX> &ex, const Matrix<SX> &var, const Matrix<SX> &expr){
  if(var.empty()) return ex; // quick return if empty
  casadi_assert_message(isSymbolic(var),"the variable is not symbolic");
  casadi_assert_message(var.size1() == expr.size1() && var.size2() == expr.size2(),"the dimensions do not match");

  // evaluate with var == expr
  SXFunction fcn(var,ex);
  return fcn.eval(expr);
}


#if 0
void replaceDerivatives(Matrix<SX> &ex, const Matrix<SX> &var, const Matrix<SX> &dvar){
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
        if(fcn.algorithm[i].op == DER){

          // find the corresponding derivative
          std::map<int, SX>::iterator r = dermap.find(fcn.algorithm[i].ch0);
          casadi_assert(r != dermap.end());

          replace[i] = r->second;
        }
  }
  Matrix<SX> res;
  Matrix<SX> repres;
  fcn.eval_symbolic(Matrix<SX>(),res,replace,repres);
  ex = res;

  casadi_assert(0);

}
#endif


void makeSmooth(Matrix<SX> &ex, Matrix<SX> &bvar, Matrix<SX> &bexpr){
  // Initialize
  SXFunction fcn(Matrix<SX>(),ex);

  casadi_assert(bexpr.empty());

  // Nodes to be replaced
  std::map<int,SX> replace;

  // Go through all nodes and check if any node is non-smooth
  for(int i=0; i<fcn->algorithm.size(); ++i){

      // Check if we have a step node
      if(fcn->algorithm[i].op == STEP){

        // Get the index of the child
        int ch0 = fcn->algorithm[i].ch[0];

        // Binary variable corresponding to the the switch
        Matrix<SX> sw;

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
  Matrix<SX> res;
  fcn->eval(Matrix<SX>(),res,replace,bexpr);

  for(int i=0; i<bexpr.size(); ++i)
    bexpr[i] = bexpr[i]->dep(0);

  ex = res;

#if 0
  // Make sure that the binding expression is smooth
  bexpr.init(Matrix<SX>());
  Matrix<SX> b;
  bexpr.eval_symbolic(Matrix<SX>(),b,replace,bexpr);
  bexpr = b;
#endif
}

Matrix<SX> spy(const Matrix<SX>& A){
  Matrix<SX> s(A.size1(),A.size2());
  for(int i=0; i<A.size1(); ++i)
    for(int j=0; j<A.size2(); ++j)
      if(!A(i,j)->isZero())
        s(i,j) = 1;
  return s;
}

int numNodes(const Matrix<SX>& A){
  // Create a function
  SXFunction fcn(Matrix<SX>(),A);

  return fcn->tree.size();
}

bool dependsOn(const Matrix<SX>& ex, const Matrix<SX> &arg){
  Matrix<SX> g = gradient(vec(ex),arg);
  return !isZero(g);
}


bool isSmooth(const Matrix<SX>& ex){
 // Make a function
 SXFunction temp(Matrix<SX>(),ex);
  
 // Run the function on the temporary variable
 return temp->isSmooth();
}


bool isSymbolic(const Matrix<SX>& ex){
  if(!isDense(ex)) return false;

  for(int k=0; k<ex.size(); ++k) // loop over non-zero elements
    if(!ex[k]->isSymbolic()) // if an element is not symbolic
      return false;
  
  return true;
}

Matrix<SX> gradient(const Matrix<SX>& ex, const Matrix<SX> &arg) {
  return trans(jacobian(ex,arg));
}
  
Matrix<SX> jacobian(const Matrix<SX>& ex, const Matrix<SX> &arg) {
  SXFunction temp(arg,ex); // make a runtime
  return temp->jac();
}

void hessian(const Matrix<SX>& ex, const Matrix<SX> &arg, Matrix<SX> &H, Matrix<SX> &g) {
  // this algorithm is _NOT_ linear time (but very easy to implement).. Change to higher order AD!
  g = gradient(ex,arg);  
  H = gradient(g,arg);
}

Matrix<SX> hessian(const Matrix<SX>& ex, const Matrix<SX> &arg) {
  Matrix<SX> H,g;
  hessian(ex,arg,H,g);
  return H;
}

double getValue(const Matrix<SX>& ex, int i, int j) {
  casadi_assert(i<ex.size1() && j<ex.size2());
  return ex(i,j)->getValue();
}

int getIntValue(const Matrix<SX>& ex, int i, int j) {
  casadi_assert(i<ex.size1() && j<ex.size2());
  return ex(i,j)->getIntValue();
}

void getValue(const Matrix<SX>& ex, double *res) {
  for(int i=0; i<ex.numel(); ++i)
    res[i] = ex(i)->getValue();
}

void getIntValue(const Matrix<SX>& ex, int *res) {
  for(int i=0; i<ex.numel(); ++i)
    res[i] = ex(i)->getIntValue();
}

const string& getName(const Matrix<SX>& ex) {
  casadi_assert_message(isScalar(ex),"the expression must be scalar");
  return ex(0)->getName();
}

istream& operator>>(istream &stream, Matrix<SX> &expr){
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
    expr = symbolic(str);
    return stream;
}

Matrix<SX> operator<=(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.getElement() <= b.getElement();
}

Matrix<SX> operator>=(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.getElement() >= b.getElement();
}

Matrix<SX> operator<(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.getElement() < b.getElement();
}

Matrix<SX> operator>(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.getElement() > b.getElement();
}

Matrix<SX> operator&&(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.getElement() && b.getElement();
}

Matrix<SX> operator||(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.getElement() || b.getElement();
}

Matrix<SX> operator==(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.getElement() == b.getElement();
}

Matrix<SX> operator!=(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.getElement() != b.getElement();
}

Matrix<SX> operator!(const Matrix<SX> &a){
  casadi_assert_message(a.scalar(), "conditional operators only defined for scalars");
  return !a.getElement();
}

Matrix<SX>& operator<<(Matrix<SX>& expr, const Matrix<SX>& add){
  append(expr,add);
  return expr;
}

void expand(const Matrix<SX>& ex2, Matrix<SX> &ww, Matrix<SX>& tt){
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
      f.push_back(casadi_limits<SX>::one.get());
    } else if(to_be_expanded.top()->isSymbolic()){ // symbolic nodes have weight one and itself as factor
      w.push_back(1);
      f.push_back(to_be_expanded.top());
    } else { // binary node

        casadi_assert(to_be_expanded.top()->hasDep()); // make sure that the node is binary

        // Check if addition, subtracton or multiplication
        BinarySXNode*     binnode = (BinarySXNode*)to_be_expanded.top();
        // If we have a binary node that we can factorize
        if(binnode->op == ADD || binnode->op == SUB || (binnode->op == MUL  && (binnode->child[0]->isConstant() || binnode->child[1]->isConstant()))){
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
          if(binnode->op == MUL){
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
            if(binnode->op == ADD){          // Addition: join both sums
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
  ww = Matrix<SX>(weights[thisind]);

  vector<SX> termsv(terms[thisind].size());
  for(int i=0; i<termsv.size(); ++i)
    termsv[i] = SX(terms[thisind][i]);
  tt = Matrix<SX>(termsv);
}

void simplify(SX& ex){
  // Start by expanding the node to a weighted sum
  Matrix<SX> terms, weights;
  expand(ex,weights,terms);

  // Make a scalar product to get the simplified expression
  Matrix<SX> s = prod(trans(weights),terms);
  ex = s.getElement();
}

void fill(Matrix<SX>& mat, const SX& val){
  if(val->isZero())    mat.makeEmpty(mat.size1(),mat.size2());
  else                 mat.makeDense(mat.size1(),mat.size2(),val);
}

// Matrix<SX> binary(int op, const Matrix<SX> &x, const Matrix<SX> &y){
//   Matrix<SX> r;
//   dynamic_cast<Matrix<SX>&>(r).binary(sfcn[op],x,y);
//   return r;
// }
// 
// Matrix<SX> scalar_matrix(int op, const SX &x, const Matrix<SX> &y){
//   Matrix<SX> r;
//   dynamic_cast<Matrix<SX>&>(r).scalar_matrix(sfcn[op],x,y);
//   return r;
// }
// 
// Matrix<SX> matrix_scalar(int op, const Matrix<SX> &x, const SX &y){
//   Matrix<SX> r;
//   dynamic_cast<Matrix<SX>&>(r).matrix_scalar(sfcn[op],x,y);
//   return r;
// }
// 
// Matrix<SX> matrix_matrix(int op, const Matrix<SX> &x, const Matrix<SX> &y){
//   Matrix<SX> r;
//   dynamic_cast<Matrix<SX>&>(r).matrix_matrix(sfcn[op],x,y);
//   return r;
// }

void make_symbolic(SX& v, const std::string& name){
  v = SX(name);
}

Matrix<SX> symbolic(const std::string& name, int n, int m){
  // Create a dense n-by-m matrix
  Matrix<SX> ret(n,m,0);
  
  // Fill with expressions
  for(int i=0; i<n; ++i){
    for(int j=0; j<m; ++j){
      stringstream ss;
      ss << name << "_" << i << "_" << j;
      ret[j+i*m] = SX(ss.str());
    }
  }
      
  // return
  return ret;
}

Matrix<SX> symbolic(const std::string& name, const CRSSparsity& sp){
  // Create a matrix
  Matrix<SX> ret(sp);
  
  // Fill with expressions
  for(int i=0; i<ret.size(); ++i){
    stringstream ss;
    ss << name << "_" << i;
    ret[i] = SX(ss.str());
  }
  
  return ret;
}

std::vector<Matrix<SX> > symbolic(const std::string& name, int n, int m, int p){
  std::vector<Matrix<SX> > ret(p);
  for(int k=0; k<p; ++k){
    stringstream ss;
    ss << name << k;
    ret[k] = symbolic(ss.str(),n,m);
  }
  return ret;
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


