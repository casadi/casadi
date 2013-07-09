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
#include "../casadi_exception.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../stl_vector_tools.hpp"
using namespace std;

namespace CasADi{
  
SXMatrix gauss_quadrature(SXMatrix f, const SXMatrix &x, const SXMatrix &a, const SXMatrix &b, int order, const SXMatrix& w){
  casadi_assert_message(order == 5, "gauss_quadrature: order must be 5");
  casadi_assert_message(w.empty(),"gauss_quadrature: empty weights");

  // Change variables to [-1,1]
  if(!a.toScalar().isEqual(-1) || !b.toScalar().isEqual(1)){
    SXMatrix q1 = (b-a)/2;
    SXMatrix q2 = (b+a)/2;

    SXFunction fcn(x,f);
    fcn.init();

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
    f_val[i] = fcn.eval(xi[i]).toScalar();

  // Weighted sum
  SX sum;
  for(int i=0; i<5; ++i)
    sum += wi[i]*f_val[i];

  return sum;
}

SXMatrix pw_const(const SXMatrix &t, const SXMatrix &tval, const SXMatrix &val){
  // number of intervals
  int n = val.numel();

  casadi_assert_message(isScalar(t),"t must be a scalar");
  casadi_assert_message(tval.numel() == n-1, "dimensions do not match");

  SXMatrix ret = val.at(0);  
  for(int i=0; i<n-1; ++i){
    ret += (val(i+1)-val(i)) * (t>=tval(i));
  }

  return ret;
}

SXMatrix pw_lin(const SX &t, const SXMatrix &tval, const SXMatrix &val){
  // Number of points
  int N = tval.numel();
  casadi_assert_message(N>=2,"pw_lin: N>=2");

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
  SXMatrix tint = tval(range(N-2),0);

  // Return piecewise linear function
  return pw_const(t, tint, lseg);
}

SXMatrix if_else(const SXMatrix &cond, const SXMatrix &if_true, const SXMatrix &if_false){
  return if_else_zero(cond,if_true) + if_else_zero(!cond,if_false);
}

SXMatrix heaviside(const SXMatrix& a){
  return (1+sign(a))/2;
}

SXMatrix ramp(const SXMatrix& a){
  return a*heaviside(a);
}

SXMatrix rectangle(const SXMatrix& a){
  return 0.5*(sign(a+0.5)-sign(a-0.5));
}

SXMatrix triangle(const SXMatrix& a){
  return rectangle(a.toScalar()/2)*(1-abs(a.toScalar()));
}

void simplify(SXMatrix &ex){
  // simplify all non-zero elements
  for(int el=0; el<ex.size(); ++el)
    simplify(ex.at(el));
}

void compress(SXMatrix &ex, int level){

  throw CasadiException("SXMatrix::compress: Not implemented");

  if(level>0)
    compress(ex,level-1);
}

std::vector<SXMatrix> substitute(const std::vector<SXMatrix> &ex, const std::vector<SXMatrix> &v, const std::vector<SXMatrix> &vdef){
  // Assert consistent dimensions
  casadi_assert(v.size()==vdef.size());

  // Quick return if all equal
  bool all_equal = true;
  for(int k=0; k<v.size(); ++k){
    if(!isEqual(v[k],vdef[k])){
      all_equal = false;
      break;
    }
  }
  if(all_equal) return ex;
  
  // Check sparsities
  for(int k=0; k<v.size(); ++k){
    if(v[k].sparsity()!=vdef[k].sparsity()) {
      if (vdef[k].scalar() && vdef[k].size()==1) { // Expand vdef to sparsity of v if vdef is scalar
        std::vector<SXMatrix> vdef_mod = vdef;
        vdef_mod[k] = SXMatrix(v[k].sparsity(),vdef[k].at(0));
        return substitute(ex,v,vdef_mod);
      } else {
        casadi_error("subsitute(ex,v,vdef): sparsities of v and vdef must match. Got v: " << v[k].dimString() << " and " << "vdef: " << vdef[k].dimString() << ".");
      }
    }
  }
  
    
  // Otherwise, evaluate symbolically     
  SXFunction F(v,ex);
  F.init();
  return F.eval(vdef);
}

SXMatrix substitute(const SXMatrix &ex, const SXMatrix &v, const SXMatrix &vdef){
  return substitute(vector<SXMatrix>(1,ex),vector<SXMatrix>(1,v),vector<SXMatrix>(1,vdef)).front();
}

Matrix<double> evalf(const SXMatrix &ex, const SXMatrix &v, const Matrix<double> &vdef) {
  SXFunction fcn(v,ex);
  fcn.init();
  fcn.input(0).set(vdef);
  fcn.evaluate();
  return fcn.output();
}

Matrix<double> evalf(const SXMatrix &ex) {
  SXFunction fcn(std::vector< SXMatrix >(0),ex);
  fcn.init();
  fcn.evaluate();
  return fcn.output();
}

void substituteInPlace(const SXMatrix &v, SXMatrix &vdef, bool reverse){
  // Empty vector
  vector<SXMatrix> ex;
  substituteInPlace(v,vdef,ex,reverse);
}

void substituteInPlace(const SXMatrix &v, SXMatrix &vdef, std::vector<SXMatrix>& ex, bool reverse){
  casadi_assert_message(isSymbolic(v),"the variable is not symbolic");
  casadi_assert_message(v.sparsity() == vdef.sparsity(),"the sparsity patterns of the expression and its defining expression do not match");
  if(v.empty()) return; // quick return if nothing to replace

  // Function inputs
  std::vector<SXMatrix> f_in;
  if(!reverse) f_in.push_back(v);

  // Function outputs
  std::vector<SXMatrix> f_out;
  f_out.push_back(vdef);
  f_out.insert(f_out.end(),ex.begin(),ex.end());
    
  // Write the mapping function
  SXFunction f(f_in,f_out);
  f.init();
  
  // Get references to the internal data structures
  const vector<ScalarAtomic>& algorithm = f.algorithm();
  vector<SX> work(f.getWorkSize());
  
  // Iterator to the binary operations
  vector<SX>::const_iterator b_it=f->operations_.begin();
  
  // Iterator to stack of constants
  vector<SX>::const_iterator c_it = f->constants_.begin();

  // Iterator to free variables
  vector<SX>::const_iterator p_it = f->free_vars_.begin();
  
  // Evaluate the algorithm
  for(vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it){
    switch(it->op){
      case OP_INPUT:
        // reverse is false, substitute out
        work[it->i0] = vdef.at(it->i2);  
        break;
      case OP_OUTPUT:
        if(it->i0==0){
          vdef.at(it->i2) = work[it->i1];
          if(reverse){
            // Use the new variable henceforth, substitute in
            work[it->i1] = v.at(it->i2);
          }
        } else {
          // Auxillary output
          ex[it->i0-1].at(it->i2) = work[it->i1];   
        }
        break;
      case OP_CONST:      work[it->i0] = *c_it++; break;
      case OP_PARAMETER:  work[it->i0] = *p_it++; break;
      default:
      {
        switch(it->op){
          CASADI_MATH_FUN_BUILTIN(work[it->i1],work[it->i2],work[it->i0])
        }
        
        // Avoid creating duplicates
        const int depth = 2; // NOTE: a higher depth could possibly give more savings
        work[it->i0].assignIfDuplicate(*b_it++,depth);
      }
    }
  }
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
        if(fcn.algorithm[i].op == DER){

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

#if 0
void makeSmooth(SXMatrix &ex, SXMatrix &bvar, SXMatrix &bexpr){
  // Initialize
  SXFunction fcn(SXMatrix(),ex);

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
    bexpr[i] = bexpr[i]->dep(0);

  ex = res;

#if 0
  // Make sure that the binding expression is smooth
  bexpr.init(SXMatrix());
  SXMatrix b;
  bexpr.eval_symbolic(SXMatrix(),b,replace,bexpr);
  bexpr = b;
#endif
}
#endif

SXMatrix spy(const SXMatrix& A){
  SXMatrix s(A.size1(),A.size2());
  for(int i=0; i<A.size1(); ++i)
    for(int j=0; j<A.size2(); ++j)
      if(!A(i,j).toScalar()->isZero())
        s(i,j) = 1;
  return s;
}

bool dependsOn(const SXMatrix& ex, const SXMatrix &arg){
  if(ex.size()==0) return false;

  SXFunction temp(arg,ex);
  temp.init();
  CRSSparsity Jsp = temp.jacSparsity();
  return Jsp.size()!=0;
}


bool isSmooth(const SXMatrix& ex){
 // Make a function
 SXFunction temp(SXMatrix(),ex);
 temp.init();
  
 // Run the function on the temporary variable
 return temp->isSmooth();
}


bool isSymbolic(const SXMatrix& ex){
  if(!isDense(ex)) return false;
  
  return isSymbolicSparse(ex);
}

bool isSymbolicSparse(const SXMatrix& ex) {
  for(int k=0; k<ex.size(); ++k) // loop over non-zero elements
    if(!ex.at(k)->isSymbolic()) // if an element is not symbolic
      return false;
  
  return true;
}

SXMatrix gradient(const SXMatrix& ex, const SXMatrix &arg) {
  SXFunction temp(arg,ex); // make a runtime
  temp.init();
  return temp.grad();
}
  
SXMatrix jacobian(const SXMatrix& ex, const SXMatrix &arg) {
  SXFunction temp(arg,ex); // make a runtime
  temp.init();
  return temp.jac();
}

void hessian(const SXMatrix& ex, const SXMatrix &arg, SXMatrix &H, SXMatrix &g) {
  g = gradient(ex,arg);  

  SXFunction temp(arg,g); // make a runtime
  temp.init();
  H = temp.jac(0,0,false,true);
}

SXMatrix hessian(const SXMatrix& ex, const SXMatrix &arg) {
  SXMatrix H,g;
  hessian(ex,arg,H,g);
  return H;
}

double getValue(const SXMatrix& ex, int i, int j) {
  casadi_assert(i<ex.size1() && j<ex.size2());
  return ex(i,j).toScalar().getValue();
}

int getIntValue(const SXMatrix& ex, int i, int j) {
  casadi_assert(i<ex.size1() && j<ex.size2());
  return ex(i,j).toScalar().getIntValue();
}

void getValue(const SXMatrix& ex, double *res) {
  for(int i=0; i<ex.numel(); ++i)
    res[i] = ex(i).toScalar()->getValue();
}

void getIntValue(const SXMatrix& ex, int *res) {
  for(int i=0; i<ex.numel(); ++i)
    res[i] = ex(i).toScalar().getIntValue();
}

const string& getName(const SXMatrix& ex) {
  casadi_assert_message(isScalar(ex),"the expression must be scalar");
  return ex.at(0)->getName();
}

void expand(const SXMatrix& ex2, SXMatrix &ww, SXMatrix& tt){
  casadi_assert(ex2.scalar());
  SX ex = ex2.toScalar();
  
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
        SXNode* node = to_be_expanded.top();
        // If we have a binary node that we can factorize
        if(node->getOp() == OP_ADD || node->getOp() == OP_SUB || (node->getOp() == OP_MUL  && (node->dep(0)->isConstant() || node->dep(1)->isConstant()))){
          // Make sure that both children are factorized, if not - add to stack
          if (indices.find(node->dep(0).get()) == indices.end()){
            to_be_expanded.push(node->dep(0).get());
            continue;
          }
          if (indices.find(node->dep(1).get()) == indices.end()){
             to_be_expanded.push(node->dep(1).get());
             continue;
          }

          // Get indices of children
          int ind1 = indices[node->dep(0).get()];
          int ind2 = indices[node->dep(1).get()];
  
          // If multiplication
          if(node->getOp() == OP_MUL){
            double fac;
            if(node->dep(0)->isConstant()){ // Multiplication where the first factor is a constant
              fac = node->dep(0)->getValue();
              f = terms[ind2];
              w = weights[ind2];
            } else { // Multiplication where the second factor is a constant
              fac = node->dep(1)->getValue();
              f = terms[ind1];
              w = weights[ind1];
            }
            for(int i=0; i<w.size(); ++i) w[i] *= fac;

          } else { // if addition or subtraction
            if(node->getOp() == OP_ADD){          // Addition: join both sums
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
        f.push_back(node);

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
  ww = SXMatrix(weights[thisind]);

  vector<SX> termsv(terms[thisind].size());
  for(int i=0; i<termsv.size(); ++i)
    termsv[i] = SX::create(terms[thisind][i]);
  tt = SXMatrix(termsv);
}

void simplify(SX& ex){
  // Start by expanding the node to a weighted sum
  SXMatrix terms, weights;
  expand(ex,weights,terms);

  // Make a scalar product to get the simplified expression
  SXMatrix s = mul(trans(weights),terms);
  ex = s.toScalar();
}

void fill(SXMatrix& mat, const SX& val){
  if(val->isZero())    mat.makeEmpty(mat.size1(),mat.size2());
  else                 mat.makeDense(mat.size1(),mat.size2(),val);
}

// SXMatrix binary(int op, const SXMatrix &x, const SXMatrix &y){
//   SXMatrix r;
//   dynamic_cast<SXMatrix&>(r).binary(sfcn[op],x,y);
//   return r;
// }
// 
// SXMatrix scalar_matrix(int op, const SX &x, const SXMatrix &y){
//   SXMatrix r;
//   dynamic_cast<SXMatrix&>(r).scalar_matrix(sfcn[op],x,y);
//   return r;
// }
// 
// SXMatrix matrix_scalar(int op, const SXMatrix &x, const SX &y){
//   SXMatrix r;
//   dynamic_cast<SXMatrix&>(r).matrix_scalar(sfcn[op],x,y);
//   return r;
// }
// 
// SXMatrix matrix_matrix(int op, const SXMatrix &x, const SXMatrix &y){
//   SXMatrix r;
//   dynamic_cast<SXMatrix&>(r).matrix_matrix(sfcn[op],x,y);
//   return r;
// }

SXMatrix ssym(const std::string& name, int n, int m){
  return ssym(name,sp_dense(n,m));
}

SXMatrix ssym(const std::string& name, const std::pair<int,int> & nm) {
  return ssym(name,nm.first,nm.second);
}

SXMatrix ssym(const std::string& name, const CRSSparsity& sp){
  // Create a dense n-by-m matrix
  vector<SX> retv;
  
  // Check if individial names have been provided
  if(name[0]=='['){

    // Make a copy of the string and modify it as to remove the special characters
    string modname = name;
    for(string::iterator it=modname.begin(); it!=modname.end(); ++it){
      switch(*it){
        case '(': case ')': case '[': case ']': case '{': case '}': case ',': case ';': *it = ' ';
      }
    }
    
    istringstream iss(modname);
    string varname;
    
    // Loop over elements
    while(!iss.fail()){
      // Read the name
      iss >> varname;
      
      // Append to the return vector
      if(!iss.fail())
        retv.push_back(SX(varname));
    }
  } else if(sp.scalar(true)){
    retv.push_back(SX(name));
  } else {
    // Scalar
    std::stringstream ss;
    for(int k=0; k<sp.size(); ++k){
      ss.str("");
      ss << name << "_" << k;
      retv.push_back(SX(ss.str()));
    }
  }

  // Determine dimensions automatically if empty
  if(sp.scalar(true)){
    return SXMatrix(retv);
  } else {
    return SXMatrix(sp,retv);
  }
}

std::vector<SXMatrix> ssym(const std::string& name, const CRSSparsity& sp, int p){
  std::vector<SXMatrix> ret(p);
  stringstream ss;
  for(int k=0; k<p; ++k){
    ss.str("");
    ss << name << "_" << k;
    ret[k] = ssym(ss.str(),sp);
  }
  return ret;
}

std::vector<std::vector<SXMatrix> > ssym(const std::string& name, const CRSSparsity& sp, int p, int r){
  std::vector<std::vector<SXMatrix> > ret(r);
  for(int k=0; k<r; ++k){
    stringstream ss;
    ss << name << "_" << k;
    ret[k] = ssym(ss.str(),sp,p);
  }
  return ret;
}

std::vector<SXMatrix> ssym(const std::string& name, int n, int m, int p){
  return  ssym(name,sp_dense(n,m),p);
}

std::vector<std::vector<SXMatrix> > ssym(const std::string& name, int n, int m, int p, int r){
  return ssym(name,sp_dense(n,m),p,r);
}

SXMatrix taylor(const SXMatrix& ex,const SXMatrix& x, const SXMatrix& a, int order) {
  casadi_assert(x.scalar() && a.scalar());
  if (ex.size()!=ex.numel())
   throw CasadiException("taylor: not implemented for sparse matrices");
  SXMatrix ff = vec(ex);
  
  SXMatrix result = substitute(ff,x,a);
  double nf=1; 
  SXMatrix dx = (x-a);
  SXMatrix dxa = (x-a);
  for (int i=1;i<=order;i++) {
    ff = jacobian(ff,x);
    nf*=i;
    result+=1/nf * substitute(ff,x,a) * dxa;
    dxa*=dx;
  }
  return trans(reshape(result,ex.size2(),ex.size1()));
}

SXMatrix mtaylor(const SXMatrix& ex,const SXMatrix& x, const SXMatrix& around,int order) {
  return mtaylor(ex,x,around,order,std::vector<int>(x.size(),1));
}

/// \cond
SXMatrix mtaylor_recursive(const SXMatrix& ex,const SXMatrix& x, const SXMatrix& a,int order,const std::vector<int>&order_contributions, const SX & current_dx=casadi_limits<SX>::one, double current_denom=1, int current_order=1) {
  SXMatrix result = substitute(ex,x,a)*current_dx/current_denom;
  for (int i=0;i<x.size();i++) {
    if (order_contributions[i]<=order) {
      result += mtaylor_recursive(
                  jacobian(ex,x.at(i)),
                  x,a,
                  order-order_contributions[i],
                  order_contributions,
                  current_dx*(x.at(i)-a.at(i)),
                  current_denom*current_order,current_order+1);
    }
  }
  return result;
}
/// \endcond

SXMatrix mtaylor(const SXMatrix& ex,const SXMatrix& x, const SXMatrix& a,int order,const std::vector<int>&order_contributions) {
  casadi_assert_message(ex.size()==ex.numel() && x.size()==x.numel(),"mtaylor: not implemented for sparse matrices");

  casadi_assert_message(x.size()==order_contributions.size(),
    "mtaylor: number of non-zero elements in x (" <<  x.size() << ") must match size of order_contributions (" << order_contributions.size() << ")"
  );

  return trans(reshape(mtaylor_recursive(vec(ex),x,a,order,order_contributions),ex.size2(),ex.size1()));
}

int countNodes(const SXMatrix& A){
  SXFunction f(SXMatrix(),A);
  f.init();
  return f.countNodes();
}


std::string getOperatorRepresentation(const SX& x, const std::vector<std::string>& args) {
  if (!x.hasDep()) throw CasadiException("getOperatorRepresentation: SX must be binary operator");
  if (args.size() == 0 || (casadi_math<double>::ndeps(x.getOp())==2 && args.size() < 2)) throw CasadiException("getOperatorRepresentation: not enough arguments supplied");
  std::stringstream s;
  casadi_math<double>::print(x.getOp(),s,args[0],args[1]);
  return s.str();
}

SXMatrix ssym(const Matrix<double>& x){
  return SXMatrix(x);
}

void makeSemiExplicit(const SXMatrix& f, const SXMatrix& x, SXMatrix& fe, SXMatrix& fi, SXMatrix& xe, SXMatrix& xi){
  casadi_assert(f.dense());
  casadi_assert(x.dense());
  
  // Create the implicit function
  SXFunction fcn(x,f);
  fcn.init();
  
  // Get the sparsity pattern of the Jacobian (no need to actually form the Jacobian)
  CRSSparsity Jsp = fcn.jacSparsity();
  
  // Free the function
  fcn = SXFunction();
  
  // Make a BLT sorting of the Jacobian (a Dulmage-Mendelsohn decomposition)
  std::vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
  Jsp.dulmageMendelsohn(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);
  
  // Make sure that the Jacobian is full rank
  casadi_assert(coarse_rowblock[0]==0);
  casadi_assert(coarse_rowblock[1]==0);
  casadi_assert(coarse_rowblock[2]==0);
  casadi_assert(coarse_rowblock[3]==coarse_rowblock[4]);

  casadi_assert(coarse_colblock[0]==0);
  casadi_assert(coarse_colblock[1]==0);
  casadi_assert(coarse_colblock[2]==coarse_colblock[3]);
  casadi_assert(coarse_colblock[3]==coarse_colblock[4]);

  // Permuted equations
  vector<SX> fp(f.size());
  for(int i=0; i<fp.size(); ++i){
    fp[i] = f.elem(rowperm[i]);
  }
  
  // Permuted variables
  vector<SX> xp(x.size());
  for(int i=0; i<xp.size(); ++i){
    xp[i]= x.elem(colperm[i]);
  }
  
  // Number of blocks
  int nb = rowblock.size()-1;

  // Block equations
  vector<SX> fb;

  // Block variables
  vector<SX> xb;

  // Block variables that enter linearly and nonlinearily respectively
  vector<SX> xb_lin, xb_nonlin;
  
  // The separated variables and equations
  vector<SX> fev, fiv, xev, xiv;
  
  // Loop over blocks
  for(int b=0; b<nb; ++b){
    
    // Get the local equations
    fb.clear();
    for(int i=rowblock[b]; i<rowblock[b+1]; ++i){
      fb.push_back(fp[i]);
    }
    
    // Get the local variables
    xb.clear();
    for(int i=colblock[b]; i<colblock[b+1]; ++i){
      xb.push_back(xp[i]);
    }

    // We shall find out which variables enter nonlinearily in the equations, for this we need a function that will depend on all the variables
    SXFunction fcnb_all(xb,inner_prod(SXMatrix(fb),ssym("dum1",fb.size())));
    fcnb_all.init();
    
    // Take the gradient of this function to find out which variables enter in the function (should be all)
    SXMatrix fcnb_dep = fcnb_all.grad();
    
    // Make sure that this expression is dense (otherwise, some variables would not enter)
    casadi_assert(fcnb_dep.dense());
    
    // Multiply this expression with a new dummy vector and take the jacobian to find out which variables enter nonlinearily
    SXFunction fcnb_nonlin(xb,inner_prod(fcnb_dep,ssym("dum2",fcnb_dep.size())));
    fcnb_nonlin.init();
    CRSSparsity sp_nonlin = fcnb_nonlin.jacSparsity();
    
    // Get the subsets of variables that appear nonlinearily
    vector<bool> nonlin(sp_nonlin.size2(),false);
    for(int el=0; el<sp_nonlin.size(); ++el){
      nonlin[sp_nonlin.col(el)] = true;
    }
/*    cout << "nonlin = " << nonlin << endl;*/
    
    // Separate variables
    xb_lin.clear();
    xb_nonlin.clear();
    for(int i=0; i<nonlin.size(); ++i){
      if(nonlin[i])
        xb_nonlin.push_back(xb[i]);
      else
        xb_lin.push_back(xb[i]);
    }
    
    // If there are only nonlinear variables
    if(xb_lin.empty()){
      // Substitute the already determined variables
      fb = substitute(SXMatrix(fb),SXMatrix(xev),SXMatrix(fev)).data();
      
      // Add to the implicit variables and equations
      fiv.insert(fiv.end(),fb.begin(),fb.end());
      xiv.insert(xiv.end(),xb.begin(),xb.end());
    } else {
      // Write the equations as a function of the linear variables
      SXFunction fcnb(xb_lin,fb);
      fcnb.init();
            
      // Write the equation in matrix form
      SXMatrix Jb = fcnb.jac();
      SXMatrix rb = -fcnb.eval(SXMatrix(xb_lin.size(),1,0));
      
      // Simple solve if there are no nonlinear variables
      if(xb_nonlin.empty()){
        
        // Check if 1-by-1 block
        if(Jb.numel()==1){
          // Simple division if Jb scalar
          rb /= Jb;
        } else {
          // Solve system of equations
          rb = solve(Jb,rb);
        }
        
        // Substitute the already determined variables
        rb = substitute(rb,SXMatrix(xev),SXMatrix(fev));
        
        // Add to the explicit variables and equations
        fev.insert(fev.end(),rb.begin(),rb.end());
        xev.insert(xev.end(),xb.begin(),xb.end());
        
      } else { // There are both linear and nonlinear variables
        
        // Make a Dulmage-Mendelsohn decomposition
        std::vector<int> rowpermb, colpermb, rowblockb, colblockb, coarse_rowblockb, coarse_colblockb;
        Jb.sparsity().dulmageMendelsohn(rowpermb, colpermb, rowblockb, colblockb, coarse_rowblockb, coarse_colblockb);
        
        Matrix<int>(Jb.sparsity(),1).printDense();
        Jb.printDense();
        
        


        

        cout << rowpermb << endl;
        cout << colpermb << endl;
        cout << rowblockb << endl;
        cout << colblockb << endl;
        cout << coarse_rowblockb << endl;
        cout << coarse_colblockb << endl;

        casadi_warning("tearing not implemented");
        
        
        // Substitute the already determined variables
        fb = substitute(SXMatrix(fb),SXMatrix(xev),SXMatrix(fev)).data();
        
        // Add to the implicit variables and equations
        fiv.insert(fiv.end(),fb.begin(),fb.end());
        xiv.insert(xiv.end(),xb.begin(),xb.end());
        
      }
    }
  }
  
  fi = SXMatrix(fiv);
  fe = SXMatrix(fev);
  xi = SXMatrix(xiv);
  xe = SXMatrix(xev);
}

SXMatrix getFree(const SXMatrix& ex){
  SXFunction f(vector<SXMatrix>(),ex);
  f.init();
  return f.getFree();
}

SXMatrix jacobianTimesVector(const SXMatrix &ex, const SXMatrix &arg, const SXMatrix &v, bool transpose_jacobian){
  SXFunction f(arg,ex);
  f.init();
  
  // Dimension of v
  int v1 = v.size1(), v2 = v.size2();
  
  // Make sure well-posed
  casadi_assert(v2 >= 1);
  casadi_assert(ex.size2()==1);
  casadi_assert(arg.size2()==1);
  if(transpose_jacobian){
    casadi_assert(v1==ex.size1());
  } else {
    casadi_assert(v1==arg.size1());
  }
  
  // Number of sensitivities
  int nfsens = transpose_jacobian ? 0 : v2;
  int nasens = transpose_jacobian ? v2 : 0;
  
  // Assemble arguments and directional derivatives
  vector<SXMatrix> argv = f.inputExpr();
  vector<SXMatrix> resv = f.outputExpr();
  vector<vector<SXMatrix> > fseed(nfsens,argv), fsens(nfsens,resv), aseed(nasens,resv), asens(nasens,argv);
  for(int dir=0; dir<v2; ++dir){
    if(transpose_jacobian){
      aseed[dir][0].set(v(Slice(0,v1),dir));
    } else {
      fseed[dir][0].set(v(Slice(0,v1),dir));
    }
  }
  
  // Evaluate with directional derivatives, output is the same as the funciton inputs
  f.evalSX(argv,resv,fseed,fsens,aseed,asens);
  
  // Get the results
  vector<SXMatrix> dirder(v2);
  for(int dir=0; dir<v2; ++dir){
    if(transpose_jacobian){
      dirder[dir] = asens[dir][0];
    } else {
      dirder[dir] = fsens[dir][0];
    }
  }
  return horzcat(dirder);
}

void extractShared(std::vector<SX>& ex, std::vector<SX>& v, std::vector<SX>& vdef, const std::string& v_prefix, const std::string& v_suffix){
  
  // Sort the expression
  SXFunction f(vector<SXMatrix>(),vector<SXMatrix>(1,ex));
  f.init();

  // Get references to the internal data structures
  const vector<ScalarAtomic>& algorithm = f.algorithm();
  vector<SX> work(f.getWorkSize());
  vector<SX> work2 = work;
  
  // Iterator to the binary operations
  vector<SX>::const_iterator b_it=f->operations_.begin();
  
  // Iterator to stack of constants
  vector<SX>::const_iterator c_it = f->constants_.begin();

  // Iterator to free variables
  vector<SX>::const_iterator p_it = f->free_vars_.begin();

  // Count how many times an expression has been used
  vector<int> usecount(work.size(),0);
  
  // New variables and definitions
  v.clear();
  vdef.clear();
  
  // Evaluate the algorithm
  for(vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it){
    // Increase usage counters
    switch(it->op){
      case OP_CONST:
      case OP_PARAMETER:
        break;
      CASADI_MATH_BINARY_BUILTIN // Binary operation
        if(usecount[it->i2]==0){
          usecount[it->i2]=1;
        } else if(usecount[it->i2]==1){
          // Get a suitable name
          vdef.push_back(work[it->i2]);
          usecount[it->i2]=-1; // Extracted, do not extract again
        }
        // fall-through
      case OP_OUTPUT: 
      default: // Unary operation, binary operation or output
        if(usecount[it->i1]==0){
          usecount[it->i1]=1;
        } else if(usecount[it->i1]==1){
          vdef.push_back(work[it->i1]);
          usecount[it->i1]=-1; // Extracted, do not extract again
        }
    }
    
    // Perform the operation
    switch(it->op){
      case OP_OUTPUT: 
        break;
      case OP_CONST:
      case OP_PARAMETER:
        usecount[it->i0] = -1; // Never extract since it is a primitive type
        break;
      default:
        work[it->i0] = *b_it++; 
        usecount[it->i0] = 0; // Not (yet) extracted
        break;
    }
  }
  
  // Create intermediate variables
  stringstream v_name;
  for(int i=0; i<vdef.size(); ++i){
    v_name.str(string());
    v_name << v_prefix << i << v_suffix;
    v.push_back(SX(v_name.str()));
  }
  
  // Mark the above expressions
  for(int i=0; i<vdef.size(); ++i){
    vdef[i].setTemp(i+1);
  }
  
  // Reset iterator
  b_it=f->operations_.begin();
  
  // Evaluate the algorithm
  for(vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it){
    switch(it->op){
      case OP_OUTPUT:     ex[it->i2] = work[it->i1];      break;
      case OP_CONST:      work2[it->i0] = work[it->i0] = *c_it++; break;
      case OP_PARAMETER:  work2[it->i0] = work[it->i0] = *p_it++; break;
      default:
      {
        switch(it->op){
          CASADI_MATH_FUN_BUILTIN(work[it->i1],work[it->i2],work[it->i0])
        }
        work2[it->i0] = *b_it++; 
        
        // Replace with intermediate variables
        int ind = work2[it->i0].getTemp()-1;
        if(ind>=0){
          vdef.at(ind) = work[it->i0];
          work[it->i0] = v.at(ind);
        }
      }
    }
  }

  // Unmark the expressions
  for(int i=0; i<vdef.size(); ++i){
    vdef[i].setTemp(0);
  }
}

void printCompact(const SXMatrix& ex, std::ostream &stream){
  // Extract shared subexpressions from ex
  vector<SX> v,vdef;
  SXMatrix ex_extracted = ex;
  extractShared(ex_extracted.data(),v,vdef,"@","");
  
  // Print the expression without shared subexpressions
  ex_extracted.print(stream);
  
  // Print the shared subexpressions
  if(!v.empty()){
    stream << endl << "where:" << endl;
    for(int i=0; i<v.size(); ++i){
      stream << v[i] << " := " << vdef[i] << endl;
    }
  }
}

  void substituteInPlace(const std::vector<SXMatrix>& v, std::vector<SXMatrix>& vdef, std::vector<SXMatrix>& ex, bool reverse){
    casadi_assert(v.size()==vdef.size());
    
    // Quick return if empty or single expression
    if(v.empty()){
      return;
     } else if(v.size()==1){
      substituteInPlace(v.front(),vdef.front(),ex,reverse);
      return;
    }

    // Count number of scalar variables
    int n =0;
    for(int i=0; i<v.size(); ++i){
      casadi_assert_message(v[i].sparsity() == vdef[i].sparsity(),"the sparsity patterns of the expression and its defining expression do not match");
      n += v[i].size();
    }

    // Gather all variables
    SXMatrix v_all(n,1,0);
    SXMatrix vdef_all(n,1,0);
    vector<SX>::iterator it_v = v_all.begin();
    vector<SX>::iterator it_vdef = vdef_all.begin();
    for(int i=0; i<v.size(); ++i){
      int nv = v[i].size();
      copy(v[i].begin(),v[i].end(),it_v);
      copy(vdef[i].begin(),vdef[i].end(),it_vdef);
      it_v += nv;  it_vdef += nv;
    }

    // Substitute
    substituteInPlace(v_all,vdef_all,ex,reverse);

    // Collect the result
    it_vdef = vdef_all.begin();
    for(int i=0; i<v.size(); ++i){
      int nv = v[i].size();
       copy(it_vdef,it_vdef+nv,vdef[i].begin());
       it_vdef += nv;
    }
  }
  
  bool isRegular(const SX& ex) {
    if (ex.isConstant()) {
      return !(ex.isNan() || ex.isInf() || ex.isMinusInf());
    } else {
      casadi_error("Cannot check regularity for symbolic SX");
    }
  }

  bool isRegular(const SXMatrix& ex) {
    // First pass: ignore symbolics
    for (int i=0;i<ex.size();++i) {
      const SX& x = ex.at(i);
      if (x.isConstant()) {
        if (x.isNan() || x.isInf() || x.isMinusInf()) return false;
      }
    }
    // Second pass: don't ignore symbolics
    for (int i=0;i<ex.size();++i) {
      if (!isRegular(ex.at(i))) return false;
    }
    return true;
  }
 
  std::vector<SX> getSymbols(const SXMatrix& e) {
    SXFunction f(std::vector<SXMatrix>(),e);
    f.init();
    return f.getFree();
  }
  
  SXMatrix poly_coeff(const SXMatrix& ex, const SXMatrix&x) {
    casadi_assert(ex.scalar());
    casadi_assert(x.scalar());
    casadi_assert(isSymbolic(x));
    
    SXMatrix ret;
    
    SXFunction f(x,ex);
    f.init();
    int mult = 1;
    bool success = false;
    for (int i=0;i<1000;++i) {
      ret.append(f.eval(casadi_limits<SX>::zero)/mult);
      SXMatrix j = f.jac();
      if (j.size()==0) {
        success = true;
        break;
      }
      f = SXFunction(x,j);
      f.init();
      mult*=i+1;
    }
    
    if (!success) casadi_error("poly: suplied expression does not appear to be polynomial.");

    std::reverse( ret.data().begin(), ret.data().end() );
    
    return ret;
    

  }
  
  SXMatrix poly_roots(const SXMatrix& p) {
    casadi_assert_message(p.size2()==1,"poly_root(): supplied paramter must be column vector but got " << p.dimString() << ".");
    casadi_assert(p.dense());
    if (p.size1()==2) { // a*x + b
      SX a = p.at(0);
      SX b = p.at(1);
      return -b/a;
    } else if (p.size1()==3) { // a*x^2 + b*x + c
      SX a = p.at(0);
      SX b = p.at(1);
      SX c = p.at(2);
      SX ds = sqrt(b*b-4*a*c);
      SX bm = -b;
      SX a2 = 2*a;
      SXMatrix ret;
      ret.append((bm-ds)/a2);
      ret.append((bm+ds)/a2);
      return ret;
    /**} else if (p.size1()==4) {
      SX a = p.at(0);
      SX b = p.at(1);
      SX c = p.at(2);
      SX d = p.at(3);
      SX ac = a*c;
      SX bb = b*b;
      SX D0 = bb-3*ac;
      SX D1 = (2*bb-9*ac)*b+27*a*a*d;
      SX D = D1*D1-4*pow(D0,3); // Assumed strictly positive
      SX C = pow((D1+sqrt(D))/2,1/3.0);
      SX DC = D0/C;
      SXMatrix ret = SXMatrix::ones(3)*b;
      double u = 1; double u_i = 1;
      ret(0)+= u*C + u_i*DC;
      u = -1/2; u_i = -0.5;
      ret(1)+= u*C + u_i*DC;
      u = -1/2; u_i = -0.5;
      ret(2)+= u*C + u_i*DC;
      return -1/ret/3/a;*/
    } else if (p.at(p.size()-1).isEqual(0)) {
      SXMatrix ret = poly_roots(p(range(p.size()-1)));
      ret.append(0);
      return ret;
    } else {
      casadi_error("poly_root(): can only solve cases for first or second order polynomial. Got order " << p.size1()-1 << ".");
    }
    
  }
  
  SXMatrix eig_symbolic(const SXMatrix& m) {
    casadi_assert_message(m.size1()==m.size2(),"eig(): supplied matrix must be square");
    
    SXMatrix ret;
    
    /// Bring m in block diagonal form, calculating eigenvalues of each block seperately
    std::vector<int> offset;
    std::vector<int> index;
    int nb = m.sparsity().stronglyConnectedComponents(offset,index);
    
    SXMatrix m_perm = m(offset,offset);
    
    SXMatrix l = ssym("l");
    
    for (int k=0;k<nb;++k) {
      std::vector<int> r = range(index.at(k),index.at(k+1));
      // det(lambda*I-m) = 0
      std::cout << r << std::endl;
      ret.append(poly_roots(poly_coeff(det(SXMatrix::eye(r.size())*l-m_perm(r,r)),l)));
    }
		
    return ret;
  }
  
} // namespace CasADi


