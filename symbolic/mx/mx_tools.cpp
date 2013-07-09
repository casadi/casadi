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

#include "mx_tools.hpp"
#include "../fx/mx_function.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../stl_vector_tools.hpp"
#include "../fx/mx_function_internal.hpp"

using namespace std;

namespace CasADi{

  MX vertcat(const vector<MX>& comp){
    return MXNode::getVertcat(comp);
  }

  std::vector<MX> vertsplit(const MX& x, const std::vector<int>& offset){
    // Consistency check
    casadi_assert(offset.size()>=1);
    casadi_assert(offset.front()==0);
    casadi_assert(offset.back()<=x.size1());
    casadi_assert(isMonotone(offset));
    
    // Trivial return if possible
    if(offset.size()==1 || (offset.size()==2 && offset.back()==x.size1())){
      return vector<MX>(1,x);
    } else {
      return x->getVertsplit(offset);
    }
  }

  MX horzcat(const vector<MX>& comp){
    vector<MX> v(comp.size());
    for(int i=0; i<v.size(); ++i)
      v[i] = trans(comp[i]);
    return trans(vertcat(v));
  }

  MX vertcat(const MX& a, const MX& b){
    vector<MX> ab;
    ab.push_back(a);
    ab.push_back(b);
    return vertcat(ab);
  }

  MX horzcat(const MX& a, const MX& b){
    vector<MX> ab;
    ab.push_back(a);
    ab.push_back(b);
    return horzcat(ab);
  }

  MX veccat(const vector<MX>& comp) {
        MX (&f)(const MX&) = vec;
    return vertcat(applymap(f,comp));
  }

  MX vecNZcat(const vector<MX>& comp) {
    MX (&f)(const MX&) = vecNZ;
    return vertcat(applymap(vecNZ,comp));
  }

  MX norm_2(const MX &x){
    return x->getNorm2();
  }

  MX norm_F(const MX &x){
    return x->getNormF();
  }

  MX norm_1(const MX &x){
    return x->getNorm1();
  }

  MX norm_inf(const MX &x){
    return x->getNormInf();
  }

  MX mul(const MX &x, const MX &y){
    return x.mul(y);
  }

  MX mul(const std::vector< MX > &args){
    casadi_assert_message(args.size()>=1,"mul(std::vector< MX > &args): supplied list must not be empty.");
    if (args.size()==1) return args[0];
    MX ret = args[0].mul(args[1]);
    for (int i=2;i<args.size();++i) {
      ret = ret.mul(args[i]);
    }
    return ret;
  }


  bool isZero(const MX& ex){
    if(ex.size()==0){
      return true;
    } else {
      return ex->isZero();
    }
  }

  bool isOne(const MX& ex){
    return !ex.isNull() && ex->isOne();
  }

  bool isMinusOne(const MX& ex){
    return !ex.isNull() && ex->isValue(-1);
  }

  bool isIdentity(const MX& ex){
    return !ex.isNull() && ex->isIdentity();
  }

  MX inner_prod(const MX &x, const MX &y){
    return x.inner_prod(y);
  }

  MX outer_prod(const MX &x, const MX &y){
    return x.outer_prod(y);
  }

  void simplify(MX& ex){
    if(!ex.isNull()){
      ex->simplifyMe(ex);
    }
  }

  bool isTranspose(const MX& ex){
    if(ex.isNull()){
      return false;
    } else {
      return ex->getOp()==OP_TRANSPOSE;
    }
  }

  MX trans(const MX &x){
    // Quick return if null or scalar
    if(x.isNull() || x.numel()==1)
      return x;
    else
      return x->getTranspose();
  }

  MX reshape(const MX &x, const std::vector<int> sz){
    if(sz.size() != 2)
      throw CasadiException("MX::reshape: not two dimensions");
    return reshape(x,sz[0],sz[1]);
  }

  MX reshape(const MX &x, int n, int m){
    if(n==x.size1() && m==x.size2())
      return x;
    else
      return reshape(x,x.sparsity().reshape(n,m));
  }

  MX reshape(const MX &x, const CRSSparsity& sp){
    // quick return if already the right shape
    if(sp==x.sparsity())
      return x;
  
    // make sure that the number of zeros agree
    casadi_assert(x.size()==sp.size());
  
    // Create a reshape node
    return x->getReshape(sp);
  }

  MX vec(const MX &x) {
    if(x.size2()==1){
      return x;
    } else {
      return reshape(trans(x),x.numel(),1);
    }
  }

  MX flatten(const MX& x) {
    if(x.size2()==1){
      return x;
    } else {
      return reshape(x,x.numel(),1);
    }
  }

  MX vecNZ(const MX& x) {
    if(x.dense()){
      return vec(x);
    } else {
      return trans(x)->getGetNonzeros(sp_dense(x.size()),range(x.size()));
    }
  }

  MX if_else(const MX &cond, const MX &if_true, const MX &if_false){
    return if_else_zero(cond,if_true) + if_else_zero(!cond,if_false);
  }

  MX unite(const MX& A, const MX& B){
    // Join the sparsity patterns
    std::vector<unsigned char> mapping;
    CRSSparsity sp = A.sparsity().patternUnion(B.sparsity(),mapping);
  
    // Split up the mapping
    std::vector<int> nzA,nzB;
  
    // Copy sparsity
    for(int k=0; k<mapping.size(); ++k){
      if(mapping[k]==1){
        nzA.push_back(k);
      } else if(mapping[k]==2){
        nzB.push_back(k);
      } else {
        throw CasadiException("Pattern intersection not empty");
      }
    }

    // Create mapping
    MX ret = MX::zeros(sp);
    ret = A->getSetNonzeros(ret,nzA);
    ret = B->getSetNonzeros(ret,nzB);
    return ret;
  }

  bool isSymbolic(const MX& ex){
    if (ex.isNull())
      return false;
    return ex->getOp()==OP_PARAMETER;
  }

  bool isSymbolicSparse(const MX& ex){
    if(ex.isNull()){
      return false;
    } else if(ex.getOp()==OP_VERTCAT){
      // Check if the expression is a vertcat where all components are symbolic primitives
      for(int d=0; d<ex->ndep(); ++d){
        if(!ex->dep(d).isSymbolic()){
          return false;
        }
      }
      return true;
    } else {
      return isSymbolic(ex);
    }
  }

  MX trace(const MX& A){
    casadi_assert_message(A.size1() == A.size2(), "trace: must be square");
    MX res(0);
    for (int i=0; i < A.size1(); i ++) {
      res+=A(i,i);
    }
    return res;
  }

  MX repmat(const MX &A, int n, int m){
    // Quick return if possible
    if(n==1 &&  m==1)
      return A;
  
    // First concatenate horizontally
    MX row = horzcat(std::vector<MX >(m, A));
  
    // Then vertically
    return vertcat(std::vector<MX >(n, row));
  }

  /**
     MX clip(const MX& A, const CRSSparsity& sp) {
     // Join the sparsity patterns
     std::vector<int> mapping;
     CRSSparsity sp = A.sparsity().patternIntersection(sp,mapping);
  
     // Split up the mapping
     std::vector<int> nzA,nzB;
  
     // Copy sparsity
     for(int k=0; k<mapping.size(); ++k){
     if(mapping[k]<0){
     nzA.push_back(k);
     } else if(mapping[k]>0){
     nzB.push_back(k);
     } else {
     throw CasadiException("Pattern intersection not empty");
     }
     }
  
     // Create mapping
     MX ret;
     ret.assignNode(new Mapping(sp));
     ret->assign(A,range(nzA.size()),nzA);
     ret->assign(B,range(nzB.size()),nzB);
     return ret;
  
     }
  */

  MX densify(const MX& x){
    MX ret = x;
    makeDense(ret);
    return ret;
  }

  void makeDense(MX& x){
    // Quick return if already dense
    if(x.dense()) return;
  
    // Densify
    x = x.setSparse(sp_dense(x.size1(),x.size2()));
  }

  MX createParent(std::vector<MX> &deps) {
    // First check if arguments are symbolic
    for (int k=0;k<deps.size();k++) {
      if (!isSymbolic(deps[k])) throw CasadiException("createParent: the argumenst must be pure symbolic");
    }
  
    // Collect the sizes of the depenencies
    std::vector<int> index(deps.size()+1,0);
    for (int k=0;k<deps.size();k++) {
      index[k+1] =  index[k] + deps[k].size();
    }
  
    // Create the parent
    MX P("P",index[deps.size()],1);
  
    // Make the arguments dependent on the parent
    for (int k=0;k<deps.size();k++) {
      deps[k] = reshape(P(range(index[k],index[k+1])),deps[k].sparsity());
    }
  
    return P;
  }

  std::pair<MX, std::vector<MX> > createParent(const std::vector<CRSSparsity> &deps) {
    // Collect the sizes of the depenencies
    std::vector<int> index(deps.size()+1,0);
    for (int k=0;k<deps.size();k++) {
      index[k+1] =  index[k] + deps[k].size();
    }
  
    // Create the parent
    MX P("P",index[deps.size()],1);
  
    std::vector<MX> ret(deps.size());
  
    // Make the arguments dependent on the parent
    for (int k=0;k<deps.size();k++) {
      ret[k] =  reshape(P(range(index[k],index[k+1])),deps[k]);
    }
  
    return std::pair< MX, std::vector<MX> > (P,ret);
  }

  std::pair<MX, std::vector<MX> > createParent(const std::vector<MX> &deps) {
    std::vector<MX> ret(deps);
    MX P = createParent(ret);
    return std::pair< MX, std::vector<MX> > (P,ret);
  }

  MX diag(const MX& x){
    // Nonzero mapping
    std::vector<int> mapping;
  
    // Get the sparsity
    CRSSparsity sp = x.sparsity().diag(mapping);
  
    // Create a reference to the nonzeros
    return x->getGetNonzeros(sp,mapping);
  }
  
  MX blkdiag(const std::vector<MX> &A) {
    // This implementation does not pretend to be efficient
    int row=0;
    int col=0;
    for (int i=0;i<A.size();++i) {
      row+=A[i].size1();
      col+=A[i].size2();
    }
    
    MX ret = MX(row,col);
    
    row = 0;
    col = 0;
    
    for (int i=0;i<A.size();++i) {
      ret(range(row,row+A[i].size1()),range(col,col+A[i].size2())) = A[i];
      row+=A[i].size1();
      col+=A[i].size2();
    }
    
    return ret;
  }
  
  MX blkdiag(const MX &A, const MX& B) {
    std::vector<MX> ret;
    ret.push_back(A);
    ret.push_back(B);
    return blkdiag(ret);
  }

  int countNodes(const MX& A){
    MXFunction f(vector<MX>(),A);
    f.init();
    return f.countNodes();
  }

  MX sumRows(const MX &x) {
    return mul(MX::ones(1,x.size1()),x);
  }

  MX sumCols(const MX &x) {
    return mul(x,MX::ones(x.size2(),1));
  }

  MX sumAll(const MX &x) {
    return sumCols(sumRows(x));
  }


  MX polyval(const MX& p, const MX& x){
    casadi_assert_message(isDense(p),"polynomial coefficients vector must be a vector");
    casadi_assert_message(isVector(p) && p.size()>0,"polynomial coefficients must be a vector");
    MX ret = p[0];
    for(int i=1; i<p.size(); ++i){
      ret = ret*x + p[i];
    }
    return ret;
  }

  bool isVector(const MX& ex){
    return ex.size2()==1;
  }

  bool isDense(const MX& ex){
    return ex.size() == ex.numel();
  }

  MX msym(const std::string& name, int n, int m){
    return MX(name,n,m);
  }

  MX msym(const std::string& name, const std::pair<int,int> & nm) {
    return MX(name,nm.first,nm.second);
  }

  MX msym(const Matrix<double>& x){
    return MX(x);
  }

  MX msym(const std::string& name, const CRSSparsity& sp) {
    return MX(name,sp);
  }

  bool isEqual(const MX& ex1,const MX &ex2){
    if ((ex1.size()!=0 || ex2.size()!=0) && (ex1.size1()!=ex2.size1() || ex1.size2()!=ex2.size2())) return false;
    MX difference = ex1 - ex2;  
    return isZero(difference);
  }

  std::string getOperatorRepresentation(const MX& x, const std::vector<std::string>& args) {
    std::stringstream s;
    const MXNode* node = dynamic_cast<const MXNode*>(x.get());
    node->printPart(s,0);
    if (x.isUnary()) {
      s << args[0];
      node->printPart(s,1);
    } else {
      for (int i=0;i<args.size();++i) {
        s << args[i];
        node->printPart(s,1+i);
      }
    }
    
    return s.str();
  }

  void substituteInPlace(const std::vector<MX>& v, std::vector<MX>& vdef, bool reverse){
    // Empty vector
    vector<MX> ex;
    substituteInPlace(v,vdef,ex,reverse);
  }

  void substituteInPlace(const std::vector<MX>& v, std::vector<MX>& vdef, std::vector<MX>& ex, bool reverse){
    casadi_assert_message(v.size()==vdef.size(),"Mismatch in the number of expression to substitute.");
    for(int k=0; k<v.size(); ++k){
      casadi_assert_message(isSymbolic(v[k]),"Variable " << k << " is not symbolic");
      casadi_assert_message(v[k].sparsity() == vdef[k].sparsity(), "Inconsistent sparsity for variable " << k << ".");
    }
    casadi_assert_message(reverse==false,"Not implemented");

    // quick return if nothing to replace
    if(v.empty()) return;

    // Function inputs
    std::vector<MX> f_in = v;

    // Function outputs
    std::vector<MX> f_out = vdef;
    f_out.insert(f_out.end(),ex.begin(),ex.end());
  
    // Write the mapping function
    MXFunction f(f_in,f_out);
    f.init();

    // Get references to the internal data structures
    std::vector<MXAlgEl>& algorithm = f->algorithm_;
    vector<MX> work(f.getWorkSize());
    MXPtrV input_p, output_p;
    MXPtrVV dummy_p;

    for(vector<MXAlgEl>::iterator it=algorithm.begin(); it!=algorithm.end(); ++it){
      switch(it->op){
      case OP_INPUT:
        work.at(it->res.front()) = vdef.at(it->arg.front());
        break;
      case OP_PARAMETER:
      case OP_CONST:
        work.at(it->res.front()) = it->data;
        break;
      case OP_OUTPUT:
        if(it->res.front()<vdef.size()){
          vdef.at(it->res.front()) = work.at(it->arg.front());
        } else {
          ex.at(it->res.front()-vdef.size()) = work.at(it->arg.front());
        }
        break;
      default:
        {
          input_p.resize(it->arg.size());
          for(int i=0; i<input_p.size(); ++i){
            int el = it->arg[i];
            input_p[i] = el<0 ? 0 : &work.at(el);
          }
        
          output_p.resize(it->res.size());
          for(int i=0; i<output_p.size(); ++i){
            int el = it->res[i];
            output_p[i] = el<0 ? 0 : &work.at(el);
          }
        
          it->data->evaluateMX(input_p,output_p,dummy_p,dummy_p,dummy_p,dummy_p,false);
        }
      }
    }  
  }

  std::vector<MX> msym(const std::string& name, const CRSSparsity& sp, int p){
    std::vector<MX> ret(p);
    for(int k=0; k<p; ++k){
      stringstream ss;
      ss << name << "_" << k;
      ret[k] = msym(ss.str(),sp);
    }
    return ret;
  }

  std::vector<std::vector<MX> > msym(const std::string& name, const CRSSparsity& sp, int p, int r){
    std::vector<std::vector<MX> > ret(r);
    for(int k=0; k<r; ++k){
      stringstream ss;
      ss << name << "_" << k;
      ret[k] = msym(ss.str(),sp,p);
    }
    return ret;
  }

  std::vector<MX> msym(const std::string& name, int n, int m, int p){
    return msym(name,sp_dense(n,m),p);
  }

  std::vector<std::vector<MX> > msym(const std::string& name, int n, int m, int p, int r){
    return msym(name,sp_dense(n,m),p,r);
  }

  MX substitute(const MX &ex, const MX& v, const MX& vdef){
    return substitute(vector<MX>(1,ex),vector<MX>(1,v),vector<MX>(1,vdef)).front();
  }

  std::vector<MX> substitute(const std::vector<MX> &ex, const std::vector<MX> &v, const std::vector<MX> &vdef){
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
  
    // Otherwise, evaluate symbolically     
    MXFunction F(v,ex);
    F.init();
    return F.evalMX(vdef);
  }

  void extractShared(std::vector<MX>& ex, std::vector<MX>& v, std::vector<MX>& vdef, const std::string& v_prefix, const std::string& v_suffix){
  
    // Sort the expression
    MXFunction f(vector<MX>(),ex);
    f.init();

    // Get references to the internal data structures
    const vector<MXAlgEl>& algorithm = f.algorithm();
    vector<MX> work(f.getWorkSize());
  
    // Count how many times an expression has been used
    vector<int> usecount(work.size(),0);
  
    // Remember the origin of every calculation
    vector<pair<int,int> > origin(work.size(),make_pair(-1,-1));
    
    // Which evaluations to replace
    vector<pair<int,int> > replace;
  
    // Evaluate the algorithm to identify which evaluations to replace
    int k=0;
    for(vector<MXAlgEl>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it, ++k){
      // Increase usage counters
      switch(it->op){
      case OP_CONST:
      case OP_PARAMETER:
        break;
      default: // Unary operation, binary operation or output
        for(int c=0; c<it->arg.size(); ++c){
          if(usecount[it->arg[c]]==0){
            usecount[it->arg[c]]=1;
          } else if(usecount[it->arg[c]]==1){
            replace.push_back(origin[it->arg[c]]);
            usecount[it->arg[c]]=-1; // Extracted, do not extract again
          }
        }
      }
        
      // Perform the operation
      switch(it->op){
      case OP_OUTPUT: 
        break;
      case OP_CONST:
      case OP_PARAMETER:
        usecount[it->res.front()] = -1; // Never extract since it is a primitive type
        break;
      default:
        for(int c=0; c<it->res.size(); ++c){
          if(it->res[c]>=0){
            work[it->res[c]] = it->data.getOutput(c);
            usecount[it->res[c]] = 0; // Not (yet) extracted
            origin[it->res[c]] = make_pair(k,c);
          }
        }
        break;
      }
    }
  
    // New variables and definitions
    v.clear();
    v.reserve(replace.size());
    vdef.clear();
    vdef.reserve(replace.size());
  
    // Quick return
    if(replace.empty()) return;
  
    // Sort the elements to be replaced in the order of appearence in the algorithm
    sort(replace.begin(),replace.end());
    vector<pair<int,int> >::const_iterator replace_it=replace.begin();
  
    // Name of intermediate variables
    stringstream v_name;
  
    // Arguments for calling the atomic operations
    MXPtrV input_p, output_p;
    MXPtrVV dummy_p;
  
    // Evaluate the algorithm
    k=0;
    for(vector<MXAlgEl>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it, ++k){
      switch(it->op){
      case OP_OUTPUT:     ex[it->res.front()] = work[it->arg.front()];      break;
      case OP_CONST:
      case OP_PARAMETER:  work[it->res.front()] = it->data; break;
      default:
        {
          // Pointers to the arguments of the evaluation
          input_p.resize(it->arg.size());
          for(int i=0; i<input_p.size(); ++i){
            int el = it->arg[i]; // index of the argument
            input_p[i] = el<0 ? 0 : &work[el];
          }
        
          // Pointers to the result of the evaluation
          output_p.resize(it->res.size());
          for(int i=0; i<output_p.size(); ++i){
            int el = it->res[i]; // index of the output
            output_p[i] = el<0 ? 0 : &work[el];
          }
        
          // Evaluate atomic operation
          const_cast<MX&>(it->data)->evaluateMX(input_p,output_p,dummy_p,dummy_p,dummy_p,dummy_p,false);
        
          // Possibly replace results with new variables
          for(int c=0; c<it->res.size(); ++c){
            int ind = it->res[c];
            if(ind>=0 && replace_it->first==k && replace_it->second==c){
              // Store the result
              vdef.push_back(work[ind]);
            
              // Create a new variable
              v_name.str(string());
              v_name << v_prefix << v.size() << v_suffix;
              v.push_back(MX(v_name.str()));
            
              // Use in calculations
              work[ind] = v.back();
            
              // Go to the next element to be replaced
              replace_it++;
            }
          }
        }
      }
    }
  }

  void printCompact(const MX& ex, std::ostream &stream){
    // Extract shared subexpressions from ex
    vector<MX> v,vdef;
    vector<MX> ex_extracted(1,ex);
    extractShared(ex_extracted,v,vdef,"@","");
  
    // Print the expression without shared subexpressions
    ex_extracted.front().print(stream);
  
    // Print the shared subexpressions
    if(!v.empty()){
      stream << endl << "where:" << endl;
      for(int i=0; i<v.size(); ++i){
        stream << v[i] << " := " << vdef[i] << endl;
      }
    }
  }

  MX jacobian(const MX& ex, const MX &arg) {
    MXFunction temp(arg,ex); // make a runtime
    temp.init();
    return temp.jac();
  }

  MX gradient(const MX& ex, const MX &arg) {
    MXFunction temp(arg,ex); // make a runtime
    temp.init();
    return temp.grad();
  }

  MX blockcat(const std::vector< std::vector<MX > > &v) {
    std::vector< MX > ret;
    for(int i=0; i<v.size(); ++i)
      ret.push_back(horzcat(v[i]));
    return vertcat(ret);
  }
  
  MX blockcat(const MX &A,const MX &B,const MX &C,const MX &D) {
    return vertcat(horzcat(A,B),horzcat(C,D));
  }

  MX det(const MX& A){
    return A->getDeterminant();
  }

  MX inv(const MX& A){
    return A->getInverse();
  }
  
  bool isRegular(const MX& ex) {
    if (ex.isConstant()) {
      return isRegular(ex.getMatrixValue());
    } else {
      casadi_error("Cannot check regularity for symbolic MX");
    }
  }
  
  std::vector<MX> getSymbols(const MX& e) {
    MXFunction f(std::vector<MX>(),e);
    f.init();
    return f.getFree();
  }


} // namespace CasADi

