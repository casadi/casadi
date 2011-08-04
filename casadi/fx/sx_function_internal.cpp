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

namespace CasADi{

using namespace std;


SXFunctionInternal::SXFunctionInternal(const vector<Matrix<SX> >& inputv_, const vector<Matrix<SX> >& outputv_) : inputv(inputv_),  outputv(outputv_) {
  addOption("symbolic_jacobian",OT_BOOLEAN,true); // generate jacobian symbolically by source code transformation
  setOption("name","unnamed_sx_function");
  
  casadi_assert(!outputv.empty());

  // Check that inputs are symbolic
  for(int i=0; i<inputv_.size(); ++i) {
    if (!isSymbolicSparse(inputv_[i])) {
      stringstream ss;
      ss << "SXFunctionInternal::SXFunctionInternal: SXfunction input arguments must be purely symbolic." << endl;
      ss << "Argument #" << i << " is not symbolic." << endl;
      throw CasadiException(ss.str());
    }
  }

  // Stack
  stack<SXNode*> s;

  // Add the inputs to the stack
  for(vector<Matrix<SX> >::const_iterator it = inputv.begin(); it != inputv.end(); ++it)
    for(vector<SX>::const_iterator itc = it->begin(); itc != it->end(); ++itc)
      s.push(itc->get());

  // Add the outputs to the stack
  for(vector<Matrix<SX> >::const_iterator it = outputv.begin(); it != outputv.end(); ++it)
    for(vector<SX>::const_iterator itc = it->begin(); itc != it->end(); ++itc)
      s.push(itc->get());

  // Order the nodes in the order of dependencies using a depth-first topological sorting
  nodes.clear();
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
  
  // Resort the nodes in a more cache friendly order (Kahn 1962)
  resort_bredth_first(bnodes);

  // Set the temporary variables to be the corresponding place in the sorted graph
  for(int i=0; i<nodes.size(); ++i)
    nodes[i]->temp = i;

  // Place in the work vector for each of the nodes in the tree
  vector<int> place(nodes.size(),-1);

  worksize = nodes.size();
  for(int i=0; i<place.size(); ++i)
    place[i] = i;


  // Add the binary operations
  algorithm.clear();
  algorithm.resize(bnodes.size());
  pder1.resize(algorithm.size());
  pder2.resize(algorithm.size());
  for(int i=0; i<algorithm.size(); ++i){

    // Locate the element
    SXNode* bnode = bnodes[i];

    // Save the node index
    algorithm[i].ind = place[bnode->temp];
    
    // Save the indices of the children
    algorithm[i].ch[0] = place[bnode->dep(0).get()->temp];
    algorithm[i].ch[1] = place[bnode->dep(1).get()->temp];

    // Operation
    algorithm[i].op = bnode->getOp();
  }

  // Indices corresponding to the inputs
  input_.resize(inputv.size());
  input_ind.resize(inputv.size());
  for(int i=0; i<inputv.size(); ++i){
    // References
    const Matrix<SX>& ip = inputv[i];
    input(i) = Matrix<double>(ip.sparsity());

    // Allocate space for the indices
    vector<int>& ii = input_ind[i];
    ii.resize(ip.size());

    // save the indices
    for(int j=0; j<ip.size(); ++j){
      ii[j] = place[ip.data()[j].get()->temp];
    }
  }

  // Indices corresponding to each non-zero outputs
  output_.resize(outputv.size());
  output_ind.resize(output_.size());
  for(int i=0; i<outputv.size(); ++i){
    // References
    const Matrix<SX>& op = outputv[i];
    output(i) = Matrix<double>(op.size1(),op.size2(),op.col(),op.rowind());
    
    // Allocate space for the indices
    vector<int>& oi = output_ind[i];  
    oi.resize(op.size());

    // save the indices
    for(int j=0; j<op.size(); ++j){
      oi[j] = place[op.data()[j].get()->temp];
    }
  }

  // Allocate a work vector
  work.resize(worksize,numeric_limits<double>::quiet_NaN());

  // Save the constants to the work vector
  for(vector<SXNode*>::iterator it=cnodes.begin(); it!=cnodes.end(); ++it)
    work[(**it).temp] = (**it).getValue();

  // Reset the temporary variables
  for(int i=0; i<nodes.size(); ++i){
    nodes[i]->temp = 0;
  }


  // Mark all the variables
  for(vector<SXMatrix>::iterator i=inputv.begin(); i!=inputv.end(); ++i){
    for(vector<SX>::iterator j=i->begin(); j!=i->end(); ++j){
      j->setTemp(1);
    }
  }
  
  // Collect free varables
  free_vars.clear();
  for(int el=0; el<snodes.size(); ++el){
    if(!snodes[el]->temp){
      free_vars.push_back(SX(snodes[el]));
    }
  }

  // Unmark variables
  for(vector<SXMatrix>::iterator i=inputv.begin(); i!=inputv.end(); ++i){
    for(vector<SX>::iterator j=i->begin(); j!=i->end(); ++j){
      j->setTemp(0);
    }
  }
}

SXFunctionInternal::~SXFunctionInternal(){
}


void SXFunctionInternal::evaluate(int nfdir, int nadir){
  if(!free_vars.empty()){
    stringstream ss;
    ss << "Cannot evaluate \"";
    repr(ss);
    ss << "\" since variables " << free_vars << " are free";
    throw CasadiException(ss.str());
  }
  
  // Copy the function arguments to the work vector
  for(int ind=0; ind<getNumInputs(); ++ind){
    const Matrix<double> &arg = input(ind);
    for(int i=0; i<arg.size(); ++i){
      work[input_ind[ind][i]] = arg.data()[i];
    }
  }
  
  // Evaluate the algorithm
  if(nfdir==0 && nadir==0){
    // without taping
    for(vector<AlgEl>::iterator it=algorithm.begin(); it<algorithm.end(); ++it){
      // Get the arguments
      double x = work[it->ch[0]];
      double y = work[it->ch[1]];
      casadi_math<double>::fun[it->op](x,y,work[it->ind]);
    }
  } else {
    // with taping
    vector<AlgEl>::iterator it = algorithm.begin();
    vector<AlgElData<1> >::iterator it1 = pder1.begin();
    for(; it<algorithm.end(); ++it, ++it1){
      // Get the arguments
      double x = work[it->ch[0]];
      double y = work[it->ch[1]];
      casadi_math<double>::fun[it->op](x,y,work[it->ind]);
      casadi_math<double>::der[it->op](x,y,work[it->ind],it1->d);
    }
  }
  
  // Get the results
  for(int ind=0; ind<getNumOutputs(); ++ind){
    Matrix<double> &res = output(ind);
    for(int i=0; i<res.size(); ++i){
      res.data()[i] = work[output_ind[ind][i]];
    }
  }

  // Loop over all forward directions
  for(int dir=0; dir<nfdir; ++dir){

    // Clear the seeds (not necessary if constants and parameters have zero value!)
    fill(dwork.begin(),dwork.end(),0);
    
    // Copy the function arguments to the work vector
    for(int ind=0; ind<input_.size(); ++ind){
      const Matrix<double> &fseed = fwdSeed(ind,dir);
      for(int i=0; i<fseed.size(); ++i){
        dwork[input_ind[ind][i]] = fseed.data()[i];
      }
    }
  
    // Evaluate the algorithm for the sensitivities
    vector<AlgEl>::const_iterator it = algorithm.begin();
    vector<AlgElData<1> >::const_iterator it2 = pder1.begin();
    for(; it!=algorithm.end(); ++it, ++it2){
      dwork[it->ind] = it2->d[0] * dwork[it->ch[0]] + it2->d[1] * dwork[it->ch[1]];
    }
  
    // Get the forward sensitivities
    for(int ind=0; ind<output_.size(); ++ind){
      Matrix<double> &fsens = fwdSens(ind,dir);
      for(int i=0; i<output_ind[ind].size(); ++i){
        fsens.data()[i] = dwork[output_ind[ind][i]];
      }
    }
  }
  
  for(int dir=0; dir<nadir; ++dir){

  // Clear the seeds (should not be necessary)
  fill(dwork.begin(),dwork.end(),0);

    // Pass the output seeds
    for(int ind=0; ind<output_.size(); ++ind){
      const Matrix<double> &aseed = adjSeed(ind,dir);
      for(int i=0; i<output_ind[ind].size(); ++i){
        dwork[output_ind[ind][i]] += aseed.data()[i];
      }
    }

  for(int i=algorithm.size()-1; i>=0; --i){
    const AlgEl& ae = algorithm[i];
    const AlgElData<1>& aed = pder1[i];
    
    // copy the seed and clear the cache entry
    double seed = dwork[ae.ind];
    dwork[ae.ind] = 0;

    dwork[ae.ch[0]] += aed.d[0] * seed;
    dwork[ae.ch[1]] += aed.d[1] * seed;
  }

  // Collect the adjoint sensitivities
  for(int ind=0; ind<getNumInputs(); ++ind){
    Matrix<double> &asens = adjSens(ind,dir);
    for(int i=0; i<input_ind[ind].size(); ++i){
      asens.data()[i] = dwork[input_ind[ind][i]];
    }
  }

  // Should clean up any zero added to constants and parameters!

  }
}

vector<Matrix<SX> > SXFunctionInternal::jac(const vector<pair<int,int> >& jblocks){
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
      ret[i] = outputv.at(oind);
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
    }
  }
  
  // Quick return if no jacobians to be calculated
  if(jblocks_no_f.empty())
    return ret;
  
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
  der1.reserve(algorithm.size());
  der2.reserve(algorithm.size());
  SX tmp[2];
  for(vector<AlgEl>::const_iterator it = algorithm.begin(); it!=algorithm.end(); ++it){
      SX f = SX(nodes[it->ind]);
      SX ch[2];
      ch[0] = SX(nodes[it->ch[0]]);
      ch[1] = SX(nodes[it->ch[1]]);
      casadi_math<SX>::der[it->op](ch[0],ch[1],f,tmp);
      if(!ch[0]->isConstant())  der1.push_back(tmp[0]);
      else                      der1.push_back(0);

      if(!ch[1]->isConstant())  der2.push_back(tmp[1]);
      else                      der2.push_back(0);
  }

  // Gradient (this is also the working array)
  vector<SX> g(nodes.size(),casadi_limits<SX>::zero);

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
        g[input_ind[iind][c]] = casadi_limits<SX>::one;
      }
    }
    
    // forward sweep
    for(int k = 0; k<algorithm.size(); ++k){
      const AlgEl& ae = algorithm[k];
      
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
          ret[jblock_ind[v]].data()[elJ] = g[output_ind[oind][r_out]];
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
          g[output_ind[oind][nz[c]]] = casadi_limits<SX>::one;
      }
    }
    
    // backward sweep
    for(int k = algorithm.size()-1; k>=0; --k){
      const AlgEl& ae = algorithm[k];
        
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
          ret[jblock_ind[v]].data()[elJ] = g[input_ind[iind][c]];
        }
      }
    }
  }
  
  // Return
  return ret;
}

bool SXFunctionInternal::isSmooth() const{
    // Go through all nodes and check if any node is non-smooth
    for(vector<AlgEl>::const_iterator it = algorithm.begin(); it!=algorithm.end(); ++it){
      if(it->op == STEP || it->op == FLOOR )
        return false;
    }
    return true;
}

void SXFunctionInternal::print(ostream &stream) const{
 for(vector<AlgEl>::const_iterator it = algorithm.begin(); it!=algorithm.end(); ++it){
    int op = it->op;
    stringstream s,s0,s1;
    s << "a" << it->ind;

    int i0 = it->ch[0], i1 = it->ch[1];
#if 0
    if(nodes[i0]->hasDep())  s0 << "a" << i0;
    else                      s0 << SX(nodes[i0]);
    if(nodes[i1]->hasDep())  s1 << "a" << i1;
    else                      s1 << SX(nodes[i1]);
#else
    s0 << "a" << i0;
    s1 << "a" << i1;
#endif

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

std::string SXFunctionInternal::printOperation(int i){
  stringstream s;
  if(nodes.at(i)->isConstant()){
    double v = nodes.at(i)->getValue();
    if(v>=0){
      s << v;
    } else {
      s << "(" << v << ")";
    }
  } else if(refcount[i]==1){
    // Call can be inlined
    int j = nodes[i]->temp;
    string ss[2];
    for(int c=0; c<2; ++c){
      ss[c] = printOperation(algorithm[j].ch[c]);
    }
    int op = algorithm[j].op;
    casadi_math<double>::print[op](s,ss[0],ss[1]);
  } else {
    s << "a" << i;
  }
  return s.str();
}

void SXFunctionInternal::generateCode(const string& src_name){
  casadi_assert(isInit());
  
   // Output
  cout << "Generating: " << src_name << " (" << algorithm.size() << " elementary operations)" << endl;
  // Create the c source file
  ofstream cfile;
  cfile.open (src_name.c_str());

  // Print header
  cfile << "/* Automatically generated function. */" << endl;
//   cfile << "#include <stdio.h>" << endl << endl;
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
  vector<bool> declared(dwork.size(),false);
  
  // Copy the function arguments to the work vector
  for(int ind=0; ind<input_.size(); ++ind){
    for(int i=0; i<input_ind[ind].size(); ++i){
      int el = input_ind[ind][i];
      cfile << "d a" << el << "=x[" << ind << "][" << i << "];" << endl;
      declared[el] = true;
    }
  }
  
 // Count how many times a node is referenced in the algorithm
 refcount.clear();
 refcount.resize(nodes.size(),0);
 for(vector<AlgEl>::const_iterator it = algorithm.begin(); it!=algorithm.end(); ++it){
   for(int c=0; c<2; ++c) refcount[it->ch[c]]++;
 }
  
 // Inputs get value two in order to prevent them from being inlined
 for(int ind=0; ind<input_.size(); ++ind){
   for(int i=0; i<input_ind[ind].size(); ++i){
     int el = input_ind[ind][i];
     refcount[el] = 2;
   }
 }

 // Outputs get value two in order to prevent them from being inlined
 for(int ind=0; ind<output_.size(); ++ind){
   for(int i=0; i<output_ind[ind].size(); ++i){
     int el = output_ind[ind][i];
     refcount[el] = 2;
   }
 }
 
 // Mark the node its place in the algorithm
 for(int i=0; i<algorithm.size(); ++i){
   nodes[algorithm[i].ind]->temp = i;
 }
 
 // Run the algorithm
 for(vector<AlgEl>::const_iterator it = algorithm.begin(); it!=algorithm.end(); ++it){
   // Skip printing variables that can be inlined
   if(refcount[it->ind]<2) continue;
   
    if(!declared[it->ind]){
      cfile << "d ";
      declared[it->ind]=true;
    }
    cfile << "a" << it->ind << "=";
    string s[2];
    for(int c=0; c<2; ++c){
      s[c] = printOperation(it->ch[c]);
    }
    int op = it->op;
    casadi_math<double>::print[op](cfile,s[0],s[1]);
    cfile  << ";" << endl;
  }

 // Unmark the nodes
 for(int i=0; i<algorithm.size(); ++i){
   nodes[algorithm[i].ind]->temp = 0;
 }

 // Clear the reference counter
 refcount.clear();

  // Get the results
  for(int ind=0; ind<output_.size(); ++ind){
    for(int i=0; i<output_ind[ind].size(); ++i){
      cfile << "r[" << ind << "][" << i << "]=";
      int el = output_ind[ind][i];
      if(nodes[el]->isConstant())
        cfile << nodes[el]->getValue() << ";" << endl;
      else
        cfile << "a" << el << ";" << endl;
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

  // allocate a vector with the values at the nodes, the first vector also contains the partial derivatives
  
  if(nfdir_>0 || nadir_>0)
    dwork.resize(worksize,numeric_limits<double>::quiet_NaN());
}

FX SXFunctionInternal::hessian(int iind, int oind){
  SXFunction f;
  f.assignNode(this);
  Matrix<SX> H = f.hess(iind,oind);
  return SXFunction(inputv,H);
}

void SXFunctionInternal::evaluateSX(const vector<Matrix<SX> >& input_s, vector<Matrix<SX> >& output_s, bool eliminate_constants){
  casadi_assert_message(inputv.size() == input_s.size(),"SXFunctionInternal::evaluateSX: wrong number of inputs");
  for(int i=0; i<inputv.size(); ++i){
    casadi_assert_message(input_s[i].size()==inputv[i].size(), "SXFunctionInternal::evaluateSX: argument nonzero number does not match");
    // casadi_assert_message(input_s[i].sparsity()==inputv[i].sparsity(), "SXFunctionInternal::evaluateSX: argument sparsity does not match");
  }
  
  // Assert input dimension
  assert(input_.size() == input_s.size());
  
  // Create a symbolic work vector if not existing
  if(swork.size() != nodes.size()){
    swork.resize(nodes.size());
    for(int i=0; i<nodes.size(); ++i){
      if(!nodes[i]->hasDep())
        swork[i] = SX(nodes[i]);
    }
  }
  
  // Copy the function arguments to the work vector
  for(int ind=0; ind<input_s.size(); ++ind){
    for(int i=0; i<input_ind[ind].size(); ++i){
      swork[input_ind[ind][i]] = input_s[ind].data()[i];
    }
  }
  
  // Evaluate the algorithm
  for(vector<AlgEl>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it){
      // Get the arguments
    SX x = swork[it->ch[0]];
    SX y = swork[it->ch[1]];
    if(eliminate_constants && x.isConstant() && y.isConstant()){
      // Check if both arguments are constants
      double temp;
      casadi_math<double>::fun[it->op](x.getValue(),y.getValue(),temp);
      swork[it->ind] = temp;
    } else {
      casadi_math<SX>::fun[it->op](x,y,swork[it->ind]);
    }
  }

  // Get the results
  for(int ind=0; ind<output_.size(); ++ind){
    for(int i=0; i<output_ind[ind].size(); ++i)
      output_s[ind].data()[i] = swork[output_ind[ind][i]];
  }
}

SXFunctionInternal* SXFunctionInternal::clone() const{
  return new SXFunctionInternal(*this);
}


void SXFunctionInternal::clearSymbolic(){
  inputv.clear();
  outputv.clear();
  swork.clear();
}

FX SXFunctionInternal::jacobian(const vector<pair<int,int> >& jblocks){
  // Jacobian blocks
  vector<SXMatrix> jac_out(jblocks.size());
  jac_out.reserve(jblocks.size());
  for(int el=0; el<jac_out.size(); ++el){
    if(jblocks[el].second==-1){
      // Undifferentiated function
      jac_out[el] = outputv[jblocks[el].first];
    } else {
      // Jacobian
      vector<pair<int,int> > jblocks_local(1,jblocks[el]);
      jac_out[el] = jac(jblocks_local).front();
    }
  }
  
  // Return function
  return SXFunction(inputv,jac_out);
}


CRSSparsity SXFunctionInternal::getJacSparsity(int iind, int oind){

  // Number of input variables (columns of the Jacobian)
  int n_in = input(iind).numel();
  
  // Number of output variables (rows of the Jacobian)
  int n_out = output(oind).numel();
  
  // Number of nonzero inputs
  int nz_in = input(iind).size();
  
  // Number of nonzero outputs
  int nz_out = output(oind).size();
    
  // Make sure that dwork, which we will now use, has been allocated
  if(dwork.size() < worksize) dwork.resize(worksize);
  
  // We need a work array containing unsigned long rather than doubles. Since the two datatypes have the same size (64 bits)
  // we can save overhead by reusing the double array
  bvec_t *iwork = get_bvec_t(dwork);
  fill_n(iwork,dwork.size(),0);

  // Number of forward sweeps we must make
  int nsweep_fwd = nz_in/bvec_size;
  if(nz_in%bvec_size>0) nsweep_fwd++;
  
  // Number of adjoint sweeps we must make
  int nsweep_adj = nz_out/bvec_size;
  if(nz_out%bvec_size>0) nsweep_adj++;

  // Sparsity of the output
  const CRSSparsity& oind_sp = output(oind).sparsity();

  // Nonzero offset
  int offset = 0;

  // Return sparsity
  CRSSparsity ret;
  
  // We choose forward or adjoint based on whichever requires less sweeps
  if(nsweep_fwd <= nsweep_adj){ // forward mode
    
    // Return value (sparsity pattern)
    CRSSparsity ret_trans(nz_in,n_out);
    
    // Loop over the variables, ndir variables at a time
    for(int s=0; s<nsweep_fwd; ++s){
      
      // Integer seed for each direction
      bvec_t b = 1;
      
      // Give seeds to a set of directions
      for(int i=0; i<bvec_size && offset+i<nz_in; ++i){
        iwork[input_ind[iind][offset+i]] = b;
        b <<= 1;
      }
      
      // Propagate the dependencies
      for(vector<AlgEl>::iterator it=algorithm.begin(); it!=algorithm.end(); ++it){
        iwork[it->ind] = iwork[it->ch[0]] | iwork[it->ch[1]];
      }

      // Dependency to be checked
      b = 1;
    
      // Loop over seed directions
      for(int i=0; i<bvec_size && offset+i<nz_in; ++i){
        
        // Loop over the rows of the output
        for(int ii=0; ii<oind_sp.size1(); ++ii){
          
          // Loop over the nonzeros of the output
          for(int el=oind_sp.rowind(ii); el<oind_sp.rowind(ii+1); ++el){
            
            // If dependents on the variable
            if(b & iwork[output_ind[oind][el]]){
              
              // Column
              int jj = oind_sp.col(el);
              
              // Add to pattern
              ret_trans.getNZ(offset+i,jj + ii*oind_sp.size2());
            }
          }
        }
        
        // Go to next dependency
        b <<= 1;
      }
      
      // Remove the seeds
      for(int i=0; i<bvec_size && offset+i<nz_in; ++i){
        iwork[input_ind[iind][offset+i]] = 0;
      }

      // Update offset
      offset += bvec_size;
    }
    
    // Return sparsity pattern
    vector<int> mapping;
    ret = ret_trans.transpose(mapping);
  } else { // Adjoint mode
    
    // Return value (sparsity pattern)
    ret = CRSSparsity(n_out,nz_in);
    
    // Loop over the variables, ndir variables at a time
    for(int s=0; s<nsweep_adj; ++s){
      
      // Integer seed for each direction
      bvec_t b = 1;
     
      // Remove all seeds
      fill_n(iwork,dwork.size(),0);

      // Give seeds to a set of directions
      for(int i=0; i<bvec_size && offset+i<nz_out; ++i){
        iwork[output_ind[oind][offset+i]] |= b; // note that we may have several nonzeros using the same entry in the work vector, therefore |=
        b <<= 1;
      }
      
      // Propagate the dependencies
      // for(vector<AlgEl>::reverse_iterator it=algorithm.rbegin(); it!=algorithm.rend(); ++it) // commented out due to bug(?) in Mac using old gcc
      for(int i=algorithm.size()-1; i>=0; --i){
        AlgEl *it = &algorithm[i];
        iwork[it->ch[0]] |= iwork[it->ind];
        iwork[it->ch[1]] |= iwork[it->ind];
      }
      
      // Output dependency to be checked
      b = 1;
    
      // Loop over seed directions
      for(int i=0; i<bvec_size && offset+i<nz_out; ++i){
        
        // Loop over the nonzeros of the input
        for(int el=0; el<nz_in; ++el){
          
          // If the output is influenced by the variable
          if(b & iwork[input_ind[iind][el]]){
            
            // Add to pattern
            ret.getNZ(offset+i,el);
            
          }
        }
        
        // Go to next dependency
        b <<= 1;
      }

      // Update offset
      offset += bvec_size;
    }
  }
  
  // Enlarge if sparse input
  if(n_in!=nz_in){
    
    // New column for each old column
    vector<int> col_map(input(iind).size());
    
    // Loop over rows of the input
    for(int i=0; i<input(iind).size1(); ++i){
      
      // Loop over the nonzeros
      for(int el=input(iind).rowind(i); el<input(iind).rowind(i+1); ++el){
        
        // Get column
        int j = input(iind).col(el);
        
        // Get the column of the jacobian
        col_map[el] = j+i*input(iind).size2();
      }
    }
    
    // Enlarge matrix
    ret.resize(ret.size1(),n_in);
    
    // Swap columns
    vector<int>& col = ret.colRef();
    for(vector<int>::iterator it=col.begin(); it!=col.end(); ++it){
      *it = col_map[*it];
    }
  }
  
  // Return sparsity pattern
  return ret;
}

} // namespace CasADi

