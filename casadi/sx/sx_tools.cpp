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
#include "../stl_vector_tools.hpp"
using namespace std;

namespace CasADi{
  
Matrix<SX> gauss_quadrature(Matrix<SX> f, const Matrix<SX> &x, const Matrix<SX> &a, const Matrix<SX> &b, int order, const Matrix<SX>& w){
  casadi_assert_message(order == 5, "gauss_quadrature: order must be 5");
  casadi_assert_message(w.empty(),"gauss_quadrature: empty weights");

  // Change variables to [-1,1]
  if(!a.toScalar().isEqual(-1) || !b.toScalar().isEqual(1)){
    Matrix<SX> q1 = (b-a)/2;
    Matrix<SX> q2 = (b+a)/2;

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

Matrix<SX> heaviside(const Matrix<SX>& a){
  return (1+sign(a))/2;
}

Matrix<SX> ramp(const Matrix<SX>& a){
  return a*heaviside(a);
}

Matrix<SX> rectangle(const Matrix<SX>& a){
  return 0.5*(sign(a+0.5)-sign(a-0.5));
}

Matrix<SX> triangle(const Matrix<SX>& a){
  return rectangle(a.toScalar()/2)*(1-abs(a.toScalar()));
}

bool contains(const Matrix<SX> &list, const SX &e) {
  for (int i=0;i<nnz(list);i++) {
    if (list(i).toScalar()->isEqual(e)) return true;
  }
  return false;
}


void simplify(Matrix<SX> &ex){
  // simplify all non-zero elements
  for(int el=0; el<ex.size(); ++el)
    simplify(ex.at(el));
}

void compress(Matrix<SX> &ex, int level){

  throw CasadiException("Matrix<SX>::compress: Not implemented");

  if(level>0)
    compress(ex,level-1);
}

std::vector<Matrix<SX> > substitute(const std::vector<Matrix<SX> > &ex, const Matrix<SX> &var, const Matrix<SX> &expr){
  SXFunction fcn(var,ex);
  fcn.init();
  return fcn.eval(vector<Matrix<SX> >(1,expr));
}

Matrix<SX> substitute(const Matrix<SX> &ex, const Matrix<SX> &var, const Matrix<SX> &expr){
  if(var.empty()) return ex; // quick return if empty
  casadi_assert_message(isSymbolic(var),"the variable is not symbolic");
  casadi_assert_message(var.size1() == expr.size1() && var.size2() == expr.size2(),"the dimensions do not match");

  // evaluate with var == expr
  SXFunction fcn(var,ex);
  fcn.init();
  return fcn.eval(expr);
}

Matrix<SX> substituteInPlace(const Matrix<SX> &var, const Matrix<SX> &def, bool eliminate_constants){
  casadi_assert_message(isSymbolic(var),"the variable is not symbolic");
  casadi_assert_message(var.size1() == def.size1() && var.size2() == def.size2(),"the dimensions do not match");
  casadi_assert_message(def.dense(),"Expression must be dense");
  if(var.empty()) return def; // quick return if nothing to replace

  // Write the mapping function
  SXFunction f(var,def);
  f.setOption("topological_sorting","depth-first");
  f.init();

  // Get references to the internal data structures
  const std::vector<int>& input_ind = f->input_ind_.front();
  std::vector<int>& output_ind = f->output_ind_.front();
  std::vector<SXAlgEl>& algorithm = f->algorithm_;

  // Create the filter
  vector<int> filter = range(f->work_.size());
  casadi_assert(input_ind.size()==output_ind.size());
  for(int k=0; k<input_ind.size(); ++k){
    filter[input_ind[k]] = output_ind[k];
  }
  
  // Start by filtering out the replacement variables from the output_ind vector
  for(vector<int>::iterator it=output_ind.begin(); it!=output_ind.end(); ++it){
    *it = filter[*it];
  }
  
  // Now filter out the variables from the algorithm
  for(vector<SXAlgEl>::iterator it=algorithm.begin(); it!=algorithm.end(); ++it){
    it->ch[0] = filter[it->ch[0]];
    it->ch[1] = filter[it->ch[1]];
  }
  
  // Replace expression
  std::vector<Matrix<SX> > outputv = f->outputv_;
  f->evaluateSX(f->inputv_, outputv, eliminate_constants);
  
  // Replace the result
  return outputv.front();
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

#if 0
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
#endif

Matrix<SX> spy(const Matrix<SX>& A){
  Matrix<SX> s(A.size1(),A.size2());
  for(int i=0; i<A.size1(); ++i)
    for(int j=0; j<A.size2(); ++j)
      if(!A(i,j).toScalar()->isZero())
        s(i,j) = 1;
  return s;
}

bool dependsOn(const Matrix<SX>& ex, const Matrix<SX> &arg){
  if(ex.size()==0) return false;

  SXFunction temp(arg,ex);
  temp.init();
  CRSSparsity Jsp = temp.jacSparsity();
  return Jsp.size()!=0;
}


bool isSmooth(const Matrix<SX>& ex){
 // Make a function
 SXFunction temp(Matrix<SX>(),ex);
 temp.init();
  
 // Run the function on the temporary variable
 return temp->isSmooth();
}


bool isSymbolic(const Matrix<SX>& ex){
  if(!isDense(ex)) return false;
  
  return isSymbolicSparse(ex);
}

bool isSymbolicSparse(const Matrix<SX>& ex) {
  for(int k=0; k<ex.size(); ++k) // loop over non-zero elements
    if(!ex.at(k)->isSymbolic()) // if an element is not symbolic
      return false;
  
  return true;
}

Matrix<SX> gradient(const Matrix<SX>& ex, const Matrix<SX> &arg) {
  return trans(jacobian(ex,arg));
}
  
Matrix<SX> jacobian(const Matrix<SX>& ex, const Matrix<SX> &arg) {
  SXFunction temp(arg,ex); // make a runtime
  temp.init();
  return temp.jac();
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
  return ex(i,j).toScalar().getValue();
}

int getIntValue(const Matrix<SX>& ex, int i, int j) {
  casadi_assert(i<ex.size1() && j<ex.size2());
  return ex(i,j).toScalar().getIntValue();
}

void getValue(const Matrix<SX>& ex, double *res) {
  for(int i=0; i<ex.numel(); ++i)
    res[i] = ex(i).toScalar()->getValue();
}

void getIntValue(const Matrix<SX>& ex, int *res) {
  for(int i=0; i<ex.numel(); ++i)
    res[i] = ex(i).toScalar().getIntValue();
}

const string& getName(const Matrix<SX>& ex) {
  casadi_assert_message(isScalar(ex),"the expression must be scalar");
  return ex(0).toScalar()->getName();
}

Matrix<SX> operator<=(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.toScalar() <= b.toScalar();
}

Matrix<SX> operator>=(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.toScalar() >= b.toScalar();
}

Matrix<SX> operator<(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.toScalar() < b.toScalar();
}

Matrix<SX> operator>(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.toScalar() > b.toScalar();
}

Matrix<SX> operator&&(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.toScalar() && b.toScalar();
}

Matrix<SX> operator||(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.toScalar() || b.toScalar();
}

Matrix<SX> operator==(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.toScalar() == b.toScalar();
}

Matrix<SX> operator!=(const Matrix<SX> &a, const Matrix<SX> &b){
  casadi_assert_message(a.scalar() && b.scalar(), "conditional operators only defined for scalars");
  return a.toScalar() != b.toScalar();
}

Matrix<SX> operator!(const Matrix<SX> &a){
  casadi_assert_message(a.scalar(), "conditional operators only defined for scalars");
  return !a.toScalar();
}

void expand(const Matrix<SX>& ex2, Matrix<SX> &ww, Matrix<SX>& tt){
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
    termsv[i] = SX::create(terms[thisind][i]);
  tt = Matrix<SX>(termsv);
}

void simplify(SX& ex){
  // Start by expanding the node to a weighted sum
  Matrix<SX> terms, weights;
  expand(ex,weights,terms);

  // Make a scalar product to get the simplified expression
  Matrix<SX> s = mul(trans(weights),terms);
  ex = s.toScalar();
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

Matrix<SX> ssym(const std::string& name){
  // Check if individial names have been provided
  if(name[0]=='['){
    // Create a return object (empty)
    Matrix<SX> ret;

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
        ret.append(SX(varname));
    }
    
    return ret;
  } else {
    // Scalar
    return SX(name);
  }
}


Matrix<SX> ssym(const std::string& name, int n, int m){
  // Create a dense n-by-m matrix
  Matrix<SX> ret(n,m,0);

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
    for(int i=0; i<n; ++i){
      for(int j=0; j<m; ++j){
        // Read the name
        iss >> varname;
        
        // Make sure that the reading was successful
        casadi_assert_message(!iss.fail(), "Not enough variable names in string");
        
        // Create the variable
        ret.data()[j+i*m] = SX(varname);
      }
    }
    
    // Try to read one more element
    iss >> varname;
    
    // Make sure that the reading failed
    casadi_assert_message(iss.fail(), "Too many variable names in string");
  } else {
  
    // Fill with expressions
    for(int i=0; i<n; ++i){
      for(int j=0; j<m; ++j){
        stringstream ss;
        ss << name;
        if(n>1)
          ss << "_" << i;
        if(m>1)
          ss << "_" << j;
        
        ret.data()[j+i*m] = SX(ss.str());
      }
    }
  }
      
  // return
  return ret;
}

Matrix<SX> ssym(const std::string& name, const std::pair<int,int> & nm) {
  return ssym(name,nm.first,nm.second);
}

Matrix<SX> ssym(const std::string& name, const CRSSparsity& sp){
  // Create a matrix
  Matrix<SX> ret(sp);
  
  // Fill with expressions
  for(int i=0; i<ret.size(); ++i){
    stringstream ss;
    ss << name << "_" << i;
    ret.data()[i] = SX(ss.str());
  }
  
  return ret;
}

std::vector<Matrix<SX> > ssym(const std::string& name, int n, int m, int p){
  std::vector<Matrix<SX> > ret(p);
  for(int k=0; k<p; ++k){
    stringstream ss;
    ss << name << "_" << k;
    ret[k] = ssym(ss.str(),n,m);
  }
  return ret;
}

std::vector<std::vector<Matrix<SX> > > ssym(const std::string& name, int n, int m, int p, int r){
  std::vector<std::vector<Matrix<SX> > > ret(r);
  for(int k=0; k<r; ++k){
    stringstream ss;
    ss << name << "_" << k;
    ret[k] = ssym(ss.str(),n,m,p);
  }
  return ret;
}

Matrix<SX> taylor(const Matrix<SX>& ex,const SX& x, const SX& a, int order) {
  if (ex.size()!=ex.numel())
   throw CasadiException("taylor: not implemented for sparse matrices");
  Matrix<SX> ff = vec(ex);
  
  Matrix<SX> result = substitute(ff,x,a);
  double nf=1; 
  SX dx = (x-a);
  SX dxa = (x-a);
  for (int i=1;i<=order;i++) {
    ff = jacobian(ff,x);
    nf*=i;
    result+=1/nf * substitute(ff,x,a) * dxa;
    dxa*=dx;
  }
  return trans(reshape(result,ex.size2(),ex.size1()));
}

Matrix<SX> mtaylor(const Matrix<SX>& ex,const Matrix<SX>& x, const Matrix<SX>& around,int order) {
  return mtaylor(ex,x,around,order,std::vector<int>(x.size(),1));
}

/// \cond
Matrix<SX> mtaylor_recursive(const Matrix<SX>& ex,const Matrix<SX>& x, const Matrix<SX>& a,int order,const std::vector<int>&order_contributions, const SX & current_dx=casadi_limits<SX>::one, double current_denom=1, int current_order=1) {
  Matrix<SX> result = substitute(ex,x,a)*current_dx/current_denom;
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

Matrix<SX> mtaylor(const Matrix<SX>& ex,const Matrix<SX>& x, const Matrix<SX>& a,int order,const std::vector<int>&order_contributions) {
  casadi_assert_message(ex.size()==ex.numel() && x.size()==x.numel(),"mtaylor: not implemented for sparse matrices");

  casadi_assert_message(x.size()==order_contributions.size(),
    "mtaylor: number of non-zero elements in x (" <<  x.size() << ") must match size of order_contributions (" << order_contributions.size() << ")"
  );

  return trans(reshape(mtaylor_recursive(vec(ex),x,a,order,order_contributions),ex.size2(),ex.size1()));
}

Matrix<SX> evaluateConstants(const Matrix<SX>& ex){
  // An empty input argument
  std::vector<Matrix<SX> > input_s;

  SXFunction fcn(input_s,ex); // Note: no input argument necessary
  fcn.init();
  
  // Evaluate symbolically, eliminating constants
  std::vector<Matrix<SX> > output_s = fcn->outputv_;
  fcn->evaluateSX(input_s,output_s,true);
  return output_s.front();
}

int countNodes(const Matrix<SX>& A){
  SXFunction f(SXMatrix(),A);
  f.init();
  return f.countNodes();
}


std::string getOperatorRepresentation(const SX& x, const std::vector<std::string>& args) {
  if (!x.isBinary()) throw CasadiException("getOperatorRepresentation: SX must be binary operator");
  if (args.size() == 0 || casadi_math<double>::ndeps(x.getOp())==2 && args.size() < 2) throw CasadiException("getOperatorRepresentation: not enough arguments supplied");
  std::stringstream s;
  casadi_math<double>::print(x.getOp(),s,args[0],args[1]);
  return s.str();
}

Matrix<SX> ssym(const Matrix<double>& x){
  return Matrix<SX>(x);
}

void makeSemiExplicit(const Matrix<SX>& f, const Matrix<SX>& x, Matrix<SX>& fe, Matrix<SX>& fi, Matrix<SX>& xe, Matrix<SX>& xi){
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

SX lift(const SX& x){
  return SX::create(new BinarySXNode(LIFT,x));
}

Matrix<SX> lift(const Matrix<SX>& x){
  // Return value
  Matrix<SX> ret = x;
  vector<SX>& retd = ret.data();
  
  // Lift all binary nodes
  for(int k=0; k<retd.size(); ++k){
    if(retd[k].isBinary()){
      retd[k] = lift(retd[k]);
    }
  }
  
  // Return the result
  return ret;
}

} // namespace CasADi


