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

#include "sx_function_internal.hpp"
#include <cassert>
#include <limits>
#include <stack>
#include <deque>
#include <fstream>
#include <sstream>
#include "../stl_vector_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../sx/sx_node.hpp"
#include "../mx/evaluation.hpp"
#include "../casadi_types.hpp"
#include "../matrix/crs_sparsity_internal.hpp"

namespace CasADi{

using namespace std;


SXFunctionInternal::SXFunctionInternal(const vector<Matrix<SX> >& inputv, const vector<Matrix<SX> >& outputv) : inputv_(inputv),  outputv_(outputv) {
  setOption("name","unnamed_sx_function");
  addOption("live_variables",OT_BOOLEAN,false,"Reuse variables in the work vector");

  casadi_assert(!outputv_.empty());

  // Check that inputs are symbolic
  for(int i=0; i<inputv.size(); ++i) {
    if (!isSymbolicSparse(inputv[i])) {
      stringstream ss;
      ss << "SXFunctionInternal::SXFunctionInternal: SXfunction input arguments must be purely symbolic." << endl;
      ss << "Argument #" << i << " is not symbolic." << endl;
      throw CasadiException(ss.str());
    }
  }
  
  // Input dimensions
  input_.resize(inputv_.size());
  for(int i=0; i<inputv_.size(); ++i){
    // References
    const Matrix<SX>& ip = inputv_[i];
    input(i) = Matrix<double>(ip.sparsity());
  }

  // Output dimensions
  output_.resize(outputv_.size());
  for(int i=0; i<outputv_.size(); ++i){
    // References
    const Matrix<SX>& op = outputv_[i];
    output(i) = Matrix<double>(op.size1(),op.size2(),op.col(),op.rowind());
  }
}

SXFunctionInternal::~SXFunctionInternal(){
}


void SXFunctionInternal::evaluate(int nfdir, int nadir){
  if(!free_vars_.empty()){
    stringstream ss;
    ss << "Cannot evaluate \"";
    repr(ss);
    ss << "\" since variables " << free_vars_ << " are free";
    throw CasadiException(ss.str());
  }
  
  // Copy the function arguments to the work vector
  for(int ind=0; ind<getNumInputs(); ++ind){
    const Matrix<double> &arg = input(ind);
    for(int i=0; i<arg.size(); ++i){
      work_[input_ind_[ind][i]] = arg.data()[i];
    }
  }
  
  // Evaluate the algorithm
  if(nfdir==0 && nadir==0){
    // without taping
    for(vector<AlgEl>::iterator it=algorithm_.begin(); it<algorithm_.end(); ++it){
      casadi_math<double>::funNew(it->op,work_[it->ch[0]],work_[it->ch[1]],work_[it->ind]);
    }
  } else {
    // with taping
    vector<AlgEl>::iterator it = algorithm_.begin();
    vector<AlgElData>::iterator it1 = pder_.begin();
    for(; it<algorithm_.end(); ++it, ++it1){
      casadi_math<double>::derFNew(it->op,work_[it->ch[0]],work_[it->ch[1]],work_[it->ind],it1->d);
    }
  }

  // Get the results
  for(int ind=0; ind<getNumOutputs(); ++ind){
    Matrix<double> &res = output(ind);
    for(int i=0; i<res.size(); ++i){
      res.data()[i] = work_[output_ind_[ind][i]];
    }
  }

  // Loop over all forward directions
  for(int dir=0; dir<nfdir; ++dir){

    // Clear the seeds (not necessary if constants and parameters have zero value!)
    fill(dwork_.begin(),dwork_.end(),0);
    
    // Copy the function arguments to the work vector
    for(int ind=0; ind<input_.size(); ++ind){
      const Matrix<double> &fseed = fwdSeed(ind,dir);
      for(int i=0; i<fseed.size(); ++i){
        dwork_[input_ind_[ind][i]] = fseed.data()[i];
      }
    }
  
    // Evaluate the algorithm_ for the sensitivities
    vector<AlgEl>::const_iterator it = algorithm_.begin();
    vector<AlgElData>::const_iterator it2 = pder_.begin();
    for(; it!=algorithm_.end(); ++it, ++it2){
      dwork_[it->ind] = it2->d[0] * dwork_[it->ch[0]] + it2->d[1] * dwork_[it->ch[1]];
    }
  
    // Get the forward sensitivities
    for(int ind=0; ind<output_.size(); ++ind){
      Matrix<double> &fsens = fwdSens(ind,dir);
      for(int i=0; i<output_ind_[ind].size(); ++i){
        fsens.data()[i] = dwork_[output_ind_[ind][i]];
      }
    }
  }
  
  // Clear the seeds since it might have been set during forward sensitivity calculations
  if(nadir>0)
    fill(dwork_.begin(),dwork_.end(),0);
  
  for(int dir=0; dir<nadir; ++dir){

    // Pass the output seeds
    for(int ind=0; ind<output_.size(); ++ind){
      const Matrix<double> &aseed = adjSeed(ind,dir);
      
      // Add seed to each entry, note that multiple outputs may point to the same variable
      for(int i=0; i<output_ind_[ind].size(); ++i){
        dwork_[output_ind_[ind][i]] += aseed.data()[i];
      }
    }
    
    #if 0 // this didn't work on a mac (why?)
    vector<AlgEl>::const_reverse_iterator it = algorithm_.rbegin();
    vector<AlgElData>::const_reverse_iterator it2 = pder_.rbegin();

    for(; it!=algorithm_.rend(); ++it, ++it2){
      // copy the seed and clear the cache entry
      double seed = dwork_[it->ind];
      dwork_[it->ind] = 0;
      dwork_[it->ch[0]] += it2->d[0] * seed;
      dwork_[it->ch[1]] += it2->d[1] * seed;
    }
    #else
    for(int i=algorithm_.size()-1; i>=0; --i){
      const AlgEl& ae = algorithm_[i];
      const AlgElData& aed = pder_[i];
      
      // copy the seed and clear the cache entry
      double seed = dwork_[ae.ind];
      dwork_[ae.ind] = 0;

      dwork_[ae.ch[0]] += aed.d[0] * seed;
      dwork_[ae.ch[1]] += aed.d[1] * seed;
    }
    #endif
    // Collect the adjoint sensitivities
    for(int ind=0; ind<getNumInputs(); ++ind){
      Matrix<double> &asens = adjSens(ind,dir);
      for(int i=0; i<input_ind_[ind].size(); ++i){
        asens.data()[i] = dwork_[input_ind_[ind][i]];
        dwork_[input_ind_[ind][i]] = 0;
      }
    }
  }
}

vector<Matrix<SX> > SXFunctionInternal::jac(const vector<pair<int,int> >& jblocks){
  if(verbose()) cout << "SXFunctionInternal::jac begin" << endl;
  
  // Create return object
  vector<Matrix<SX> > ret(jblocks.size());
  
  // Which blocks are involved in the Jacobian
  vector<pair<int,int> > jblocks_no_f;
  vector<int> jblock_ind;
  
  // Add the information we already know
  for(int i=0; i<ret.size(); ++i){
    // Get input/output indices for the block
    int oind = jblocks[i].first;
    int iind = jblocks[i].second;

    // Check if nondifferentiated variable
    if(iind<0){
      // Save to output
      ret[i] = outputv_.at(oind);
    } else { // Jacobian block
      // Mark block for evaluation
      jblocks_no_f.push_back(jblocks[i]);
      jblock_ind.push_back(i);
      
      // Make sure that the function as well as variables are vectors
      // -> Why? ticket #50   - disabling for now
      //assert(input(iind).size2()==1);
      //assert(output(oind).size2()==1);
      
      // Save sparsity
      ret[i] = SXMatrix(jacSparsity(iind,oind));
      if(verbose()){
        cout << "SXFunctionInternal::jac Block " << i << " has " << ret[i].size() << " nonzeros out of " << ret[i].numel() << " elements" << endl;
      }
    }
  }

  if(verbose()) cout << "SXFunctionInternal::jac allocated return value" << endl;
  
  // Quick return if no jacobians to be calculated
  if(jblocks_no_f.empty()){
    if(verbose()) cout << "SXFunctionInternal::jac end 1" << endl;
    return ret;
  }
  
  // Get a bidirectional partition
  vector<CRSSparsity> D1(jblocks_no_f.size()), D2(jblocks_no_f.size());
  getPartition(jblocks_no_f,D1,D2);
  
  // Bugfix for dealing with sparse jacobians
  vector<int> nz;
  if(!D2[0].isNull()){
    casadi_assert(jblocks_no_f.size()==1);
    int oind = jblocks_no_f.front().first;
    
    // Output sparsity
    const CRSSparsity& sp = output(oind).sparsity();
    
    // The nonzero element for each element
    nz.resize(sp.numel(),-1);
    for(int i=0; i<sp.size1(); ++i){
      for(int el=sp.rowind(i); el<sp.rowind(i+1); ++el){
        int c = sp.col(el);
        nz[c + i*sp.size2()] = el;
      }
    }
  }
  
  // Calculate the partial derivatives
  vector<SX> der1, der2;
  der1.reserve(algorithm_.size());
  der2.reserve(algorithm_.size());
  SX tmp[2];
  vector<SX>::const_iterator f_it=binops_.begin();
  for(vector<AlgEl>::const_iterator it = algorithm_.begin(); it!=algorithm_.end(); ++it, ++f_it){
      const SX& f = *f_it;
      const SX& x = f->dep(0);
      const SX& y = f->dep(1);
      casadi_math<SX>::derNew(it->op,x,y,f,tmp);
      if(!x->isConstant())  der1.push_back(tmp[0]);
      else                  der1.push_back(0);

      if(!y->isConstant())  der2.push_back(tmp[1]);
      else                  der2.push_back(0);
  }

  // Gradient (this is also the working array)
  vector<SX> g(swork_.size(),casadi_limits<SX>::zero);

  // Get the number of forward and adjoint sweeps
  int nfwd = D1.front().isNull() ? 0 : D1.front().size1();
  int nadj = D2.front().isNull() ? 0 : D2.front().size1();
  
  // Get transposes and mappings for all jacobian sparsity patterns if we are using forward mode
  vector<vector<int> > mapping;
  vector<CRSSparsity> sp_trans;
  if(nfwd>0){
    mapping.resize(jblock_ind.size());
    sp_trans.resize(jblock_ind.size());
    for(int i=0; i<jblock_ind.size(); ++i){
      sp_trans[i] = ret[jblock_ind[i]].sparsity().transpose(mapping[i]);
    }
  }

  // Carry out the forward sweeps
  for(int sweep=0; sweep<nfwd; ++sweep){
    
    // Remove Seeds
    fill(g.begin(),g.end(),0);
    
    // For all the input variables
    for(int v=0; v<D1.size(); ++v){

      // Get the input index
      int iind = jblocks_no_f[v].second;

      // For all the directions
      for(int el = D1[v].rowind(sweep); el<D1[v].rowind(sweep+1); ++el){

        // Get column of the Jacobian (i.e. input non-zero)
        int c = D1[v].col(el);

        // Give a seed in the direction
        g[input_ind_[iind][c]] = casadi_limits<SX>::one;
      }
    }
    
    // forward sweep
    for(int k = 0; k<algorithm_.size(); ++k){
      const AlgEl& ae = algorithm_[k];
      
      // First argument
      if(!g[ae.ch[0]]->isZero() &&  !der1[k]->isZero()){
        g[ae.ind] += der1[k] * g[ae.ch[0]];
      }
      
      // Second argument
      if(!g[ae.ch[1]]->isZero() &&  !der2[k]->isZero()){
        g[ae.ind] += der2[k] * g[ae.ch[1]];
      }
    }
   
    // For all the input variables
    for(int v=0; v<D1.size(); ++v){

      // Get the output index
      int oind = jblocks_no_f[v].first;

      // For all the columns of the Jacobian treated in the sweep
      for(int el = D1[v].rowind(sweep); el<D1[v].rowind(sweep+1); ++el){

        // Get column of the Jacobian
        int c = D1[v].col(el);
        
        // Loop over the nonzero elements in column c
        for(int el_out = sp_trans[v].rowind(c); el_out<sp_trans[v].rowind(c+1); ++el_out){
          
          // Get the row
          int r_out = sp_trans[v].col(el_out);
          
          // The nonzero of the Jacobian now treated
          int elJ = mapping[v][el_out];
          
          // Get the output seed
          casadi_assert(oind<output_ind_.size());
          casadi_assert(r_out<output_ind_[oind].size());
          
          ret[jblock_ind[v]].data()[elJ] = g[output_ind_[oind][r_out]];
        }
      }
    }
  }

  // And now the adjoint sweeps
  for(int sweep=0; sweep<nadj; ++sweep){
    
    // Remove Seeds
    fill(g.begin(),g.end(),0);

    // For all the output blocks
    for(int v=0; v<D2.size(); ++v){

      // Get the output index
      int oind = jblocks_no_f[v].first;

      // For all the directions
      for(int el = D2[v].rowind(sweep); el<D2[v].rowind(sweep+1); ++el){

        // Get the direction
        int c = D2[v].col(el);

        // Give a seed in the direction
        if(nz[c]>=0) // FIXME: ugly trick
          g[output_ind_[oind][nz[c]]] = casadi_limits<SX>::one;
      }
    }
    
    // backward sweep
    for(int k = algorithm_.size()-1; k>=0; --k){
      const AlgEl& ae = algorithm_[k];
        
      // Skip if the seed is zero
      if(!g[ae.ind]->isZero()){
        
        if(!der1[k]->isZero())
          g[ae.ch[0]] += der1[k] * g[ae.ind];
          
        if(!der2[k]->isZero())
          g[ae.ch[1]] += der2[k] * g[ae.ind];
          
        // Remove the seed
        g[ae.ind] = casadi_limits<SX>::zero;
      }
    }
    
    // For all the output blocks
    for(int v=0; v<D2.size(); ++v){

      // Get the input index
      int iind = jblocks_no_f[v].second;

      // For all the rows of the Jacobian treated in the sweep
      for(int el = D2[v].rowind(sweep); el<D2[v].rowind(sweep+1); ++el){

        // Get row of the Jacobian
        int r = D2[v].col(el);
        if(nz[r]<0) continue; // FIXME: ugly trick

        // Loop over the nonzero elements in row r
        for(int elJ = ret[jblock_ind[v]].sparsity().rowind(r); elJ<ret[jblock_ind[v]].sparsity().rowind(r+1); ++elJ){
          
          // Get the input variable (i.e. column of the Jacobian)
          int c = ret[jblock_ind[v]].sparsity().col(elJ);
          
          // Get the input seed
          ret[jblock_ind[v]].data()[elJ] = g[input_ind_[iind][c]];
        }
      }
    }
  }
  
  // Return
  if(verbose()) cout << "SXFunctionInternal::jac end" << endl;
  return ret;
}

bool SXFunctionInternal::isSmooth() const{
  casadi_assert(isInit());
    // Go through all nodes and check if any node is non-smooth
    for(vector<AlgEl>::const_iterator it = algorithm_.begin(); it!=algorithm_.end(); ++it){
      if(it->op == STEP || it->op == FLOOR )
        return false;
    }
    return true;
}

void SXFunctionInternal::print(ostream &stream) const{
 for(vector<AlgEl>::const_iterator it = algorithm_.begin(); it!=algorithm_.end(); ++it){
    int op = it->op;
    stringstream s,s0,s1;
    s << "a" << it->ind;

    int i0 = it->ch[0], i1 = it->ch[1];
    s0 << "a" << i0;
    s1 << "a" << i1;
    stream << s.str() << " = ";
    casadi_math<double>::print[op](stream,s0.str(),s1.str());
    stream << ";" << endl;
  }
}

void SXFunctionInternal::printVector(std::ostream &cfile, const std::string& name, const std::vector<int>& v){
  cfile << "int " << name << "[] = {";
  for(int i=0; i<v.size(); ++i){
    if(i!=0) cfile << ",";
    cfile << v[i];
  }
  cfile << "};" << endl;
}

void SXFunctionInternal::printOperation(std::ostream &stream, int i) const{
  const AlgEl& ae = algorithm_[i];
  const SX& f = binops_[i];
  casadi_math<double>::printPre[ae.op](stream);
  for(int c=0; c<casadi_math<double>::ndeps[ae.op]; ++c){
    if(c==1) casadi_math<double>::printSep[ae.op](stream);

    if(f->dep(c)->isConstant()){
      double v = f->dep(c)->getValue();
      if(v>=0){
        stream << v;
      } else {
        stream << "(" << v << ")";
      }
    } else if(f->dep(c)->hasDep() && refcount_[f->dep(c)->temp-1]==1) {
      printOperation(stream,f->dep(c)->temp-1);
    } else {
      stream << "a" << ae.ch[c];
    }
  }
  casadi_math<double>::printPost[ae.op](stream);
}

void SXFunctionInternal::generateCode(const string& src_name){
  casadi_assert(isInit());
  
   // Output
  if(verbose()){
    cout << "Generating: " << src_name << " (" << algorithm_.size() << " elementary operations)" << endl;
  }
  // Create the c source file
  ofstream cfile;
  cfile.open (src_name.c_str());
  cfile.precision(numeric_limits<double>::digits10+2);
  cfile << scientific; // This is really only to force a decimal dot, would be better if it can be avoided
  
  // Print header
  cfile << "/* This function was automatically generated by CasADi */" << endl;
  cfile << "#include <math.h>" << endl << endl;
  
  // Space saving macro
  cfile << "#define d double" << endl << endl;

  // Number of inputs/outputs
  int n_i = input_.size();
  int n_o = output_.size();
  int n_io = n_i + n_o;

  // Dimensions
  cfile << "int n_in_ = " << n_i << ";" << endl;
  cfile << "int n_out_ = " << n_o << ";" << endl;

  // Number of rows and columns
  vector<int> nrow(n_io), ncol(n_io);
  for(int i=0; i<n_i; ++i){
    nrow[i] = input(i).size1();
    ncol[i] = input(i).size2();
  }
  for(int i=0; i<n_o; ++i){
    nrow[i+n_i] = output(i).size1();
    ncol[i+n_i] = output(i).size2();
  }
  
  // Print to file
  printVector(cfile,"nrow_",nrow);
  printVector(cfile,"ncol_",ncol);
  
  // Print row offsets
  for(int i=0; i<n_io; ++i){
    stringstream name;
    name << "rowind_" << i << "_";
    const vector<int>& rowind = i<n_i ? input(i).rowind() : output(i-n_i).rowind();
    printVector(cfile,name.str(),rowind);
  }
  
  // Array of pointers to the arrays above
  cfile << "int *rowind_[] = {";
  for(int i=0; i<n_io; ++i){
    if(i!=0) cfile << ",";
    cfile << "rowind_" << i << "_"; 
  }
  cfile << "};" << endl;
  
  // Print columns
  for(int i=0; i<n_io; ++i){
    stringstream name;
    name << "col_" << i << "_";
    const vector<int>& col = i<n_i ? input(i).col() : output(i-n_i).col();
    printVector(cfile,name.str(),col);
  }
  
  // Array of pointers to the arrays above
  cfile << "int *col_[] = {";
  for(int i=0; i<n_io; ++i){
    if(i!=0) cfile << ",";
    cfile << "col_" << i << "_"; 
  }
  cfile << "};" << endl << endl;
  
  // Function to get dimensions
  cfile << "int init(int *n_in, int *n_out){" << endl;
  cfile << "  *n_in = n_in_;" << endl;
  cfile << "  *n_out = n_out_;" << endl;
  cfile << "  return 0;" << endl;
  cfile << "}" << endl << endl;

  // Input sizes
  cfile << "int getSparsity(int i, int *nrow, int *ncol, int **rowind, int **col){" << endl;
  cfile << "  *nrow = nrow_[i];" << endl;
  cfile << "  *ncol = ncol_[i];" << endl;
  cfile << "  *rowind = rowind_[i];" << endl;
  cfile << "  *col = col_[i];" << endl;
  cfile << "  return 0;" << endl;
  cfile << "}" << endl << endl;

  // Evaluate function
  cfile << "int evaluate(const double** x, double** r){" << endl;

  // Which variables have been declared
  vector<bool> declared(dwork_.size(),false);
  
  // Copy the function arguments to the work vector
  for(int ind=0; ind<input_.size(); ++ind){
    for(int i=0; i<input_ind_[ind].size(); ++i){
      int el = input_ind_[ind][i];
      cfile << "d a" << el << "=x[" << ind << "][" << i << "];" << endl;
      declared[el] = true;
    }
  }

 // For each node, mark its place in the algorithm
 for(int i=0; i<binops_.size(); ++i){
   binops_[i]->temp = i+1;
 }

 // Count how many times a node is referenced in the algorithm_
 refcount_.resize(binops_.size());
 fill(refcount_.begin(),refcount_.end(),0);
 for(int i=0; i<algorithm_.size(); ++i){
   for(int c=0; c<2; ++c){
     // Get the place in the algorithm of the child
     int i_ch = binops_[i]->dep(c)->temp-1;
     if(i_ch>=0){
       // If the child is a binary operation
       refcount_[i_ch]++;
     }
   }
 }
  
 // Mark the inputs with negative numbers indicating where they are stored in the work vector
 for(int ind=0; ind<input_.size(); ++ind){
   for(int el=0; el<inputv_[ind].size(); ++el){
     inputv_[ind].at(el)->temp = -1-input_ind_[ind][el];
   }
 }

 // Outputs get extra reference in order to prevent them from being inlined
 for(int ind=0; ind<output_.size(); ++ind){
   for(int el=0; el<outputv_[ind].size(); ++el){
     int i = outputv_[ind].at(el)->temp-1;
     if(i>=0){
       refcount_[i] += 2;
     }
   }
 }
 
 // Run the algorithm_
 vector<SX>::const_iterator f_it=binops_.begin();
 int ii=0;
 for(vector<AlgEl>::const_iterator it = algorithm_.begin(); it!=algorithm_.end(); ++it, ++f_it, ++ii){
   // Skip printing variables that can be inlined
   if(refcount_[ii]<2) continue;
   
   // Get a pointer to the expression
   const SX& f = *f_it;
   
    if(!declared[it->ind]){
      cfile << "d ";
      declared[it->ind]=true;
    }
    cfile << "a" << it->ind << "=";
    printOperation(cfile,ii);
    cfile  << ";" << endl;
  }

 // Clear the reference counter
 refcount_.clear();

  // Get the results
  for(int ind=0; ind<output_.size(); ++ind){
    for(int el=0; el<outputv_[ind].size(); ++el){
      cfile << "r[" << ind << "][" << el << "]=";
      const SX& f = outputv_[ind].at(el);
      int i=f->temp;
      if(i==0){ // constant
        cfile << f->getValue();
      } else if(i>0){ // binary node
        cfile << "a" << algorithm_[i-1].ind;
      } else { // input feedthrough
        cfile << "a" << (-i-1);
      }
      cfile << ";" << endl;
    }
  }

 // Unmark the binary operations
 for(int i=0; i<binops_.size(); ++i){
   binops_[i]->temp = 0;
 }
 
  // Unmark the inputs 
 for(int ind=0; ind<input_.size(); ++ind){
   for(int el=0; el<inputv_[ind].size(); ++el){
     inputv_[ind].at(el)->temp = 0;
   }
 }

  cfile << "return 0;" << endl;
  cfile << "}" << endl << endl;
  
  // Close the results file
  cfile.close();
}

void SXFunctionInternal::init(){
  // Call the init function of the base class
  XFunctionInternal::init();

    // Stack
  stack<SXNode*> s;

  // Add the inputs to the stack
  for(vector<Matrix<SX> >::const_reverse_iterator it = inputv_.rbegin(); it != inputv_.rend(); ++it)
    for(vector<SX>::const_reverse_iterator itc = it->rbegin(); itc != it->rend(); ++itc)
      s.push(itc->get());

  // Add the outputs to the stack
  for(vector<Matrix<SX> >::const_reverse_iterator it = outputv_.rbegin(); it != outputv_.rend(); ++it)
    for(vector<SX>::const_reverse_iterator itc = it->rbegin(); itc != it->rend(); ++itc)
      s.push(itc->get());

  // Order the nodes in the order of dependencies using a depth-first topological sorting
  std::vector<SXNode*> nodes;
  sort_depth_first(s,nodes);

  // Sort the nodes by type
  vector<SXNode*> bnodes, cnodes, snodes;
  for(vector<SXNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
    SXNode* t = *it;
    if(t->isConstant())
      cnodes.push_back(t);
    else if(t->isSymbolic())
      snodes.push_back(t);
    else
      bnodes.push_back(t);
  }
  
  // Get the sortign algorithm
  bool breadth_first_search;
  if(getOption("topological_sorting")=="breadth-first"){
    breadth_first_search = true;
  } else if(getOption("topological_sorting")=="depth-first"){
    breadth_first_search = false;
  } else {
    stringstream ss;
    ss << "Unrecongnized topological_sorting: " << getOption("topological_sorting") << endl;
    throw CasadiException(ss.str());
  }
  
  // Resort the nodes in a more cache friendly order (Kahn 1962)
  if(breadth_first_search){
    resort_breadth_first(bnodes);
  }

  // Set the temporary variables to be the corresponding place in the sorted graph
  for(int i=0; i<nodes.size(); ++i)
    nodes[i]->temp = i;

  // Place in the work vector for each of the nodes in the tree
  vector<int> place(nodes.size(),-1);

  // Use live variables?
  bool live_variables = getOption("live_variables");
  if(live_variables){
    if(verbose()){
      casadi_warning("Live variables is currently not compatible with symbolic calculations");
    }
    
    // Count the number of times each node is used
    vector<int> refcount(nodes.size(),0);
    for(int i=0; i<bnodes.size(); ++i){
      for(int c=0; c<2; ++c){
        refcount[bnodes[i]->dep(c)->temp]++;
      }
    }
    
    // Add artificial count to the outputs to avoid them being overwritten
    for(int ind=0; ind<output_.size(); ++ind){
      for(int el=0; el<outputv_[ind].size(); ++el){
        refcount[outputv_[ind].at(el)->temp]++;
      }
    }

    // Stack with unused elements in the work vector
    stack<int> unused;

    // Place the work vector
    int p=0;

    #if 1 // sort by type, then by location in node vector
    // Place the constants in the work vector
    for(vector<SXNode*>::iterator it=cnodes.begin(); it!=cnodes.end(); ++it){
      place[(**it).temp] = p++;
      refcount[(**it).temp]++; // make sure constants are not overwritten
    }
    
    // Then all symbolic variables
    for(vector<SXNode*>::iterator it=snodes.begin(); it!=snodes.end(); ++it){
      place[(**it).temp] = p++;
    }
   
     
    for(vector<SXNode*>::iterator it=bnodes.begin(); it!=bnodes.end(); ++it){
      int i = (**it).temp;
      
      // decrease reference count of children
      for(int c=1; c>=0; --c){ // reverse order so that the first argument will end up at the top of the stack
        int ch_ind = (*it)->dep(c)->temp;
        int remaining = --refcount[ch_ind];
        if(remaining==0) unused.push(place[ch_ind]);
      }
      
      // Try to reuse a variable from the stack if possible
      if(!unused.empty()){
        place[i] = unused.top();
        unused.pop();
      } else {
        place[i] = p++;
      }
    }
    
    #else
    // Use the same sorting as the node vector, but with live variables
    // FIXME:2011-08-09:Joel:Doesn't work, c_code_generation.py example fails.
    
    // Loop over the nodes
    for(vector<SXNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
      if((*it)->isConstant()){
        place[(**it).temp] = p++;
      } else if((*it)->isSymbolic()) {
        place[(**it).temp] = p++;
      } else {
        int i = (**it).temp;
        
        // decrease reference count of children
        for(int c=0; c<2; ++c){
          int ch_ind = (*it)->dep(c)->temp;
          int remaining = --refcount[ch_ind];
          if(remaining==0) unused.push(place[ch_ind]);
        }
        
        // Try to reuse a variable from the stack if possible
        if(!unused.empty()){
          place[i] = unused.top();
          unused.pop();
        } else {
          place[i] = p++;
        }
      }
    }
    #endif
    
    if(verbose()){
      cout << "Using live variables. Worksize is " << p << " instead of " << nodes.size() << "(" << bnodes.size() << " elementary operations)" << endl;
    }
    worksize_ = p;
  } else {
    worksize_ = nodes.size();
    for(int i=0; i<place.size(); ++i)
      place[i] = i;
  }

//  casadi_assert(0);

  // Allocate a vector containing expression corresponding to each binary operation
  binops_.resize(bnodes.size());
  for(int i=0; i<bnodes.size(); ++i){
    binops_[i] = SX::create(bnodes[i]);
  }
  bnodes.clear();
  
  // Add the binary operations
  algorithm_.clear();
  algorithm_.resize(binops_.size());
  pder_.resize(algorithm_.size());
  for(int i=0; i<algorithm_.size(); ++i){

    // Save the node index
    algorithm_[i].ind = place[binops_[i]->temp];
    
    // Save the indices of the children
    algorithm_[i].ch[0] = place[binops_[i]->dep(0).get()->temp];
    algorithm_[i].ch[1] = place[binops_[i]->dep(1).get()->temp];

    // Operation
    algorithm_[i].op = binops_[i]->getOp();
  }

  // Indices corresponding to the inputs
  input_ind_.resize(inputv_.size());
  for(int i=0; i<inputv_.size(); ++i){
    // References
    const Matrix<SX>& ip = inputv_[i];

    // Allocate space for the indices
    vector<int>& ii = input_ind_[i];
    ii.resize(ip.size());

    // save the indices
    for(int j=0; j<ip.size(); ++j){
      ii[j] = place[ip.data()[j].get()->temp];
    }
  }

  // Indices corresponding to each non-zero outputs
  output_ind_.resize(output_.size());
  for(int i=0; i<outputv_.size(); ++i){
    // References
    const Matrix<SX>& op = outputv_[i];
    
    // Allocate space for the indices
    vector<int>& oi = output_ind_[i];  
    oi.resize(op.size());

    // save the indices
    for(int j=0; j<op.size(); ++j){
      oi[j] = place[op.data()[j].get()->temp];
    }
  }

  // Allocate work vectors (symbolic/numeric)
  work_.resize(worksize_,numeric_limits<double>::quiet_NaN());
  swork_.resize(worksize_);

  // Save the constants to the work vector
  for(vector<SXNode*>::iterator it=cnodes.begin(); it!=cnodes.end(); ++it){
    work_[place[(**it).temp]] = (**it).getValue();
    swork_[place[(**it).temp]] = SX::create(*it);
  }

  // Save the symbolic variables to the symbolic work vector
  for(vector<SXNode*>::iterator it=snodes.begin(); it!=snodes.end(); ++it){
    swork_[place[(**it).temp]] = SX::create(*it);
  }

  // Reset the temporary variables
  for(int i=0; i<nodes.size(); ++i){
    nodes[i]->temp = 0;
  }


  // Mark all the variables
  for(vector<SXMatrix>::iterator i=inputv_.begin(); i!=inputv_.end(); ++i){
    for(vector<SX>::iterator j=i->begin(); j!=i->end(); ++j){
      j->setTemp(1);
    }
  }
  
  // Collect free varables
  free_vars_.clear();
  for(int el=0; el<snodes.size(); ++el){
    if(!snodes[el]->temp){
      free_vars_.push_back(SX::create(snodes[el]));
    }
  }

  // Unmark variables
  for(vector<SXMatrix>::iterator i=inputv_.begin(); i!=inputv_.end(); ++i){
    for(vector<SX>::iterator j=i->begin(); j!=i->end(); ++j){
      j->setTemp(0);
    }
  }

  if(nfdir_>0 || nadir_>0){
    dwork_.resize(worksize_,numeric_limits<double>::quiet_NaN());
  }
  
  // Get the full Jacobian already now
  if(jac_for_sens_){
    getFullJacobian();
  }
  
  // Print
  if(verbose()){
    cout << "SXFunctionInternal::init Initialized " << getOption("name") << " (" << algorithm_.size() << " elementary operations)" << endl;
  }
}

FX SXFunctionInternal::hessian(int iind, int oind){
  if(output(oind).numel() != 1)
    throw CasadiException("SXFunctionInternal::hess: function must be scalar");
  
  // Reverse mode to calculate gradient
  if(verbose()){
    cout << "SXFunctionInternal::hessian: calculating gradient " << endl;
  }
  SXFunction this_;
  this_.assignNode(this);
  Matrix<SX> g = this_.grad(iind,oind);

  if(verbose()){
    cout << "SXFunctionInternal::hessian: calculating gradient done " << endl;
  }

  makeDense(g);
  if(verbose()){
    cout << "SXFunctionInternal::hessian: made gradient dense (workaround!) " << endl;
  }

  // Numeric or symbolic hessian
  bool numeric_hessian = getOption("numeric_hessian");
  
  // Hessian function
  FX hfcn;
  
  if(numeric_hessian){
    // Create function including all inputs
    SXFunction gfcn(inputv_,g);
    gfcn.setOption("verbose",getOption("verbose"));
    gfcn.setOption("numeric_jacobian",true);
    gfcn.setOption("number_of_fwd_dir",getOption("number_of_fwd_dir"));
    gfcn.setOption("number_of_adj_dir",getOption("number_of_adj_dir"));
    gfcn.init();
    
    if(verbose()) cout << "SXFunctionInternal::hessian: calculating numeric Jacobian " << endl;
    hfcn = static_cast<FX&>(gfcn).jacobian(iind,0); // TODO: SXFunction::Jacobian cannot return SXFunction!!!
    
  } else {
    // Create function including only input to be differentiated
    SXFunction gfcn(inputv_.at(iind),g);
    gfcn.setOption("verbose",getOption("verbose"));
    gfcn.setOption("numeric_jacobian",false);
    gfcn.init();
    
    if(verbose()) cout << "SXFunctionInternal::hessian: calculating symbolic Jacobian" << endl;
    SXMatrix h = gfcn.jac(0,0);
    hfcn = SXFunction(inputv_,h);
  }
  
  // Calculate jacobian of gradient
  if(verbose()) cout << "SXFunctionInternal::hessian: calculating Hessian done" << endl;
  
  // Return jacobian of the gradient
  return hfcn;
}

void SXFunctionInternal::evaluateSX(const vector<Matrix<SX> >& input_s, vector<Matrix<SX> >& output_s, bool eliminate_constants){
  casadi_assert_message(inputv_.size() == input_s.size(),"SXFunctionInternal::evaluateSX: wrong number of inputs");
  for(int i=0; i<inputv_.size(); ++i){
    casadi_assert_message(input_s[i].size()==inputv_[i].size(), "SXFunctionInternal::evaluateSX: argument nonzero number does not match");
    // casadi_assert_message(input_s[i].sparsity()==inputv_[i].sparsity(), "SXFunctionInternal::evaluateSX: argument sparsity does not match");
  }
  
  // Assert input dimension
  assert(input_.size() == input_s.size());

  // Copy the function arguments to the work vector
  for(int ind=0; ind<input_s.size(); ++ind){
    for(int i=0; i<input_ind_[ind].size(); ++i){
      swork_[input_ind_[ind][i]] = input_s[ind].data()[i];
    }
  }
  
  // Evaluate the algorithm
  for(vector<AlgEl>::const_iterator it=algorithm_.begin(); it<algorithm_.end(); ++it){

    // Get the arguments
    SX x = swork_[it->ch[0]];
    SX y = swork_[it->ch[1]];
    if(eliminate_constants && x.isConstant() && y.isConstant()){
      // Check if both arguments are constants
      double temp;
      casadi_math<double>::funNew(it->op,x.getValue(),y.getValue(),temp);
      swork_[it->ind] = temp;
    } else {
      casadi_math<SX>::funNew(it->op,x,y,swork_[it->ind]);
    }
  }

  // Get the results
  for(int ind=0; ind<output_.size(); ++ind){
    for(int i=0; i<output_ind_[ind].size(); ++i){
      output_s[ind].data()[i] = swork_[output_ind_[ind][i]];
    }
  }
}

SXFunctionInternal* SXFunctionInternal::clone() const{
  return new SXFunctionInternal(*this);
}


void SXFunctionInternal::clearSymbolic(){
  inputv_.clear();
  outputv_.clear();
  swork_.clear();
}

FX SXFunctionInternal::jacobian(const vector<pair<int,int> >& jblocks){
  // Jacobian blocks
  vector<SXMatrix> jac_out(jblocks.size());
  jac_out.reserve(jblocks.size());
  for(int el=0; el<jac_out.size(); ++el){
    int oind = jblocks[el].first;
    int iind = jblocks[el].second;
    
    if(iind==-1){
      // Undifferentiated function
      jac_out[el] = outputv_.at(oind);
    } else {
      // Jacobian (workaround)
      vector<pair<int,int> > jblocks_local(1,jblocks[el]);
      jac_out[el] = jac(jblocks_local).front();
    }
  }
  
  // Return function
  return SXFunction(inputv_,jac_out);
}


CRSSparsity SXFunctionInternal::getJacSparsity(int iind, int oind){
  if(verbose()) cout << "SXFunctionInternal::getJacSparsity begin (iind == " << iind <<", oind == " << oind << ")" << endl;

  // Number of input variables (columns of the Jacobian)
  int n_in = input(iind).numel();
  
  // Number of output variables (rows of the Jacobian)
  int n_out = output(oind).numel();
  
  // Number of nonzero inputs
  int nz_in = input(iind).size();
  
  // Number of nonzero outputs
  int nz_out = output(oind).size();
    
  // Make sure that dwork_, which we will now use, has been allocated
  if(dwork_.size() < worksize_) dwork_.resize(worksize_);
  
  // We need a work array containing unsigned long rather than doubles. Since the two datatypes have the same size (64 bits)
  // we can save overhead by reusing the double array
  bvec_t *iwork = get_bvec_t(dwork_);
  fill_n(iwork,dwork_.size(),0);

  // Number of forward sweeps we must make
  int nsweep_fwd = nz_in/bvec_size;
  if(nz_in%bvec_size>0) nsweep_fwd++;
  
  // Number of adjoint sweeps we must make
  int nsweep_adj = nz_out/bvec_size;
  if(nz_out%bvec_size>0) nsweep_adj++;

  // Sparsity of the output
  const CRSSparsity& oind_sp = output(oind).sparsity();
  int oind_d1 = oind_sp.size1();
  int oind_d2 = oind_sp.size2();
  const vector<int>& oind_rowind = oind_sp.rowind();
  const vector<int>& oind_col = oind_sp.col();
  
  // Nonzero offset
  int offset = 0;

  // Return sparsity
  CRSSparsity ret;
  
  // Progress
  int progress = -10;

  // Temporary vectors
  vector<int> detected_row, detected_col;
  vector<int> temp(bvec_size+1); // needed in a binsort of the rows
  
  // We choose forward or adjoint based on whichever requires less sweeps
  if(nsweep_fwd <= nsweep_adj){ // forward mode
    if(verbose()) cout << "SXFunctionInternal::getJacSparsity: using forward mode: " << nsweep_fwd << " sweeps needed for " << nz_in << " directions" << endl;
    
    // Return value (sparsity pattern)
    CRSSparsity ret_trans(0,n_out);
    ret_trans.reserve(std::max(nz_in,nz_out),n_in);
    vector<int>& ret_trans_rowind = ret_trans.rowindRef();
    vector<int>& ret_trans_col = ret_trans.colRef();
    
    // Loop over the variables, ndir variables at a time
    for(int s=0; s<nsweep_fwd; ++s){
      // Print progress
      if(verbose()){
        int progress_new = (s*100)/nsweep_fwd;
        // Print when entering a new decade
        if(progress_new / 10 > progress / 10){
          progress = progress_new;
          cout << progress << " %%"  << endl;
        }
      }
      
      // Integer seed for each direction
      bvec_t b = 1;
      
      // Give seeds to a set of directions
      for(int i=0; i<bvec_size && offset+i<nz_in; ++i){
        iwork[input_ind_[iind][offset+i]] = b;
        b <<= 1;
      }

      // Propagate the dependencies
      for(vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it){
        iwork[it->ind] = iwork[it->ch[0]] | iwork[it->ch[1]];
      }
      
      // Get the number of nonzeros before adding the new ones
      int nz_offset = ret_trans_col.size();
      
      // Reset the vectors holding the nonzero elements
      detected_row.clear();
      detected_col.clear();

      // Number of local seed directions
      int ndir_local = std::min(bvec_size,nz_in-offset);
      
      // Loop over the rows of the output
      for(int ii=0; ii<oind_d1; ++ii){
          
        // Loop over the nonzeros of the output
        for(int el=oind_rowind[ii]; el<oind_rowind[ii+1]; ++el){
          
          // If there is a dependency in any of the directions
          if(0 != iwork[output_ind_[oind][el]]){
          
            // Dependency to be checked
            b = 1;
      
            // Loop over seed directions
            for(int i=0; i<ndir_local; ++i){
              
              // If dependents on the variable
              if(b & iwork[output_ind_[oind][el]]){
                
                // Column
                int jj = oind_col[el];
                
                // Add to pattern
                detected_row.push_back(i);
                detected_col.push_back(jj + ii*oind_d2);
              }
              
              // Go to next dependency
              b <<= 1;
            }
          }
        }
      }
      
      // Now check how many elements there are in each row
      fill_n(temp.begin(),ndir_local+1,0);
      for(int i=0; i<detected_row.size(); ++i){
        temp[1+detected_row[i]]++;
      }
      
      // Now make a cumsum to get offset for each row
      for(int i=0; i<ndir_local; ++i){
        temp[i+1] += temp[i];
      }
      
      // Loop over the elements again, sorting the nodes
      ret_trans_col.resize(nz_offset+detected_col.size());
      for(int i=0; i<detected_row.size(); ++i){
        // Get the place to put the element
        int el = temp[detected_row[i]]++;
        
        // Get the column
        int j = detected_col[i];
        
        // Save to sparsity pattern
        ret_trans_col[nz_offset+el] = j;
      }
      
      // Now add the offset to the rows
      for(int i=0; i<ndir_local; ++i){
        temp[i] += nz_offset;
      }
      
      // And finally add to the rowind vector and update dimension
      ret_trans_rowind.insert(ret_trans_rowind.end(),temp.begin(),temp.begin()+ndir_local);
      ret_trans->nrow_ = ret_trans_rowind.size()-1;
      
      // Remove the seeds
      for(int i=0; i<bvec_size && offset+i<nz_in; ++i){
        iwork[input_ind_[iind][offset+i]] = 0;
      }

      // Update offset
      offset += bvec_size;
    }
    
    // Return sparsity pattern
    vector<int> mapping;
    ret = ret_trans.transpose(mapping);
  } else { // Adjoint mode
    if(verbose()) cout << "SXFunctionInternal::getJacSparsity: using adjoint mode: " << nsweep_adj << " sweeps needed for " << nz_out << " directions" << endl;
    
    // Return value (sparsity pattern)
    ret = CRSSparsity(0,nz_in);
    ret.reserve(std::max(nz_in,nz_out),n_out);
    vector<int>& ret_rowind = ret.rowindRef();
    vector<int>& ret_col = ret.colRef();
    
    // Loop over the variables, ndir variables at a time
    for(int s=0; s<nsweep_adj; ++s){
      
      // Print progress
      if(verbose()){
        int progress_new = (s*100)/nsweep_adj;
        // Print when entering a new decade
        if(progress_new / 10 > progress / 10){
          progress = progress_new;
          cout << progress << " %%"  << endl;
        }
      }
      
      // Integer seed for each direction
      bvec_t b = 1;
     
      // Remove all seeds
      fill_n(iwork,dwork_.size(),0);

      // Give seeds to a set of directions
      for(int i=0; i<bvec_size && offset+i<nz_out; ++i){
        iwork[output_ind_[oind][offset+i]] |= b; // note that we may have several nonzeros using the same entry in the work vector, therefore |=
        b <<= 1;
      }
      
      // Propagate the dependencies
      // for(vector<AlgEl>::reverse_iterator it=algorithm_.rbegin(); it!=algorithm_.rend(); ++it) // commented out due to bug(?) in Mac using old gcc
      for(int i=algorithm_.size()-1; i>=0; --i){
        AlgEl *it = &algorithm_[i];
        iwork[it->ch[0]] |= iwork[it->ind];
        iwork[it->ch[1]] |= iwork[it->ind];
      }
      
      // Get the number of nonzeros before adding the new ones
      int nz_offset = ret_col.size();
      
      // Reset the vectors holding the nonzero elements
      detected_row.clear();
      detected_col.clear();

      // Number of local seed directions
      int ndir_local = std::min(bvec_size,nz_out-offset);
      
      // Loop over the nonzeros of the input
      for(int el=0; el<nz_in; ++el){

        // If there is a dependency in any of the directions
        if(0 != iwork[input_ind_[iind][el]]){
          
          // Output dependency to be checked
          b = 1;
          
          // Loop over seed directions
          for(int i=0; i<ndir_local ; ++i){
              
            // If the output is influenced by the variable
            if(b & iwork[input_ind_[iind][el]]){
            
              // Add to pattern
              detected_row.push_back(i);
              detected_col.push_back(el);
            }

            // Go to next dependency
            b <<= 1;
          }
        }
      }
          
      // Now check how many elements there are in each row
      fill_n(temp.begin(),ndir_local+1,0);
      for(int i=0; i<detected_row.size(); ++i){
        temp[1+detected_row[i]]++;
      }
      
      // Now make a cumsum to get offset for each row
      for(int i=0; i<ndir_local; ++i){
        temp[i+1] += temp[i];
      }
      
      // Loop over the elements again, sorting the nodes
      ret_col.resize(nz_offset+detected_col.size());
      for(int i=0; i<detected_row.size(); ++i){
        // Get the place to put the element
        int el = temp[detected_row[i]]++;
        
        // Get the column
        int j = detected_col[i];
        
        // Save to sparsity pattern
        ret_col[nz_offset+el] = j;
      }
      
      // Now add the offset to the rows
      for(int i=0; i<ndir_local; ++i){
        temp[i] += nz_offset;
      }
      
      // And finally add to the rowind vector and update dimension
      ret_rowind.insert(ret_rowind.end(),temp.begin(),temp.begin()+ndir_local);
      ret->nrow_ = ret_rowind.size()-1;
      
      // Update offset
      offset += bvec_size;
    }
  }


  // Enlarge if sparse output
  if(n_out!=ret.size1()){
    casadi_assert(ret.size1()==nz_out);
    
    // New row for each old row
    vector<int> row_map = output(oind).sparsity().getElementMapping();

    // Insert rows
    ret.enlargeRows(n_out,row_map);
  }

  
  // Enlarge if sparse input
  if(n_in!=ret.size2()){
    casadi_assert(ret.size2()==nz_in);
    
    // New column for each old column
    vector<int> col_map = input(iind).sparsity().getElementMapping();
    
    // Insert columns
    ret.enlargeColumns(n_in,col_map);
  }
  
  // Return sparsity pattern
  if(verbose()) cout << "SXFunctionInternal::getJacSparsity end " << endl;
  return ret;
}

} // namespace CasADi

