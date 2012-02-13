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


SXFunctionInternal::SXFunctionInternal(const vector<Matrix<SX> >& inputv, const vector<Matrix<SX> >& outputv) : 
  XFunctionInternalCommon<SXFunctionInternal,Matrix<SX>,SXNode>(inputv,outputv) {
  setOption("name","unnamed_sx_function");
  addOption("live_variables",OT_BOOLEAN,false,"Reuse variables in the work vector");
  addOption("inplace",OT_BOOLEAN,false,"Evaluate with inplace operations (experimental)");

  casadi_assert(!outputv_.empty());
}

SXFunctionInternal::~SXFunctionInternal(){
}


void SXFunctionInternal::evaluate(int nfdir, int nadir){
  // NOTE: The implementation of this function is very delicate. Small changes in the class structure
  // can cause large performance losses. For this reason, the preprocessor macros are used below
  
  if (!free_vars_.empty()) {
    std::stringstream ss;
    repr(ss);
    casadi_error("Cannot evaluate \"" << ss.str() << "\" since variables " << free_vars_ << " are free.");
  }
  // Copy the function arguments to the work vector
  for(int ind=0; ind<getNumInputs(); ++ind){
    const Matrix<double> &arg = input(ind);
    for(int i=0; i<arg.size(); ++i){
      work_[input_ind_[ind][i]] = arg.data()[i];
    }
  }

  // Shorthands
  #define x work_[it->ch[0]]
  #define y work_[it->ch[1]]
  #define f work_[it->ind]

  // Evaluate the algorithm
  if(nfdir==0 && nadir==0){ // without taping
    if(evaluate_inplace_){ // allow inplace operations
      
      for(vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it){
        // NOTE: This is equivalent to casadi_math<double>::inplacefun(it->op,x,y,f);
        // but forces the function to be inlined, which is important for speed here
        CASADI_MATH_INPLACEFUN(double,it->op,x,y,f)
      }
    } else { // no inplace operations
    
      for(vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it){
        // NOTE: This is equivalent to casadi_math<double>::fun(it->op,x,y,f);
        // but forces the function to be inlined, which is important for speed here
        CASADI_MATH_FUN(double,it->op,x,y,f)
      }
    }
  } else { // with taping
    
    vector<AlgEl>::iterator it = algorithm_.begin();
    vector<AlgElData>::iterator it1 = pder_.begin();
    for(; it!=algorithm_.end(); ++it, ++it1){
      // NOTE: This is equivalent to casadi_math<double>::derF(it->op,x,y,f,it1->d);
      // but forces the function to be inlined, which is important for speed here
      CASADI_MATH_DERF(double,it->op,x,y,f,it1->d)
    }
  }
  
  // Undefine shorthands
  #undef x
  #undef y
  #undef f

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

vector<Matrix<SX> > SXFunctionInternal::jac(const vector<pair<int,int> >& jblocks, bool compact, const vector<bool>& symmetric_block){
  return jacGen(jblocks,compact,symmetric_block);
}

bool SXFunctionInternal::isSmooth() const{
  assertInit();
    // Go through all nodes and check if any node is non-smooth
    for(vector<AlgEl>::const_iterator it = algorithm_.begin(); it!=algorithm_.end(); ++it){
      if(it->op == STEP || it->op == FLOOR )
        return false;
    }
    return true;
}

void SXFunctionInternal::print(ostream &stream) const{
 FXInternal::print(stream);
 
 for(vector<AlgEl>::const_iterator it = algorithm_.begin(); it!=algorithm_.end(); ++it){
    int op = it->op;
    int ip = op/NUM_BUILT_IN_OPS;
    op -= NUM_BUILT_IN_OPS*ip;
    
    stringstream s,s0,s1;
    s << "a" << it->ind;

    int i0 = it->ch[0], i1 = it->ch[1];
    s0 << "a" << i0;
    s1 << "a" << i1;
    
    stream << s.str();
    switch(ip){
      case 0:  stream << " = "; break;
      case 1:  stream << " += "; break;
      case 2:  stream << " -= "; break;
      case 3:  stream << " *= "; break;
      case 4:  stream << " /= "; break;
    }
    
    casadi_math<double>::print(op,stream,s0.str(),s1.str());
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
  casadi_math<double>::printPre(ae.op,stream);
  for(int c=0; c<casadi_math<double>::ndeps(ae.op); ++c){
    if(c==1) casadi_math<double>::printSep(ae.op,stream);

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
  casadi_math<double>::printPost(ae.op,stream);
}

void SXFunctionInternal::generateCode(const string& src_name){
  assertInit();
  
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

  // The sign function
  cfile << "double sign(double x){ return x<0 ? -1 : x>0 ? 1 : x;}" << endl << endl;
  
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

 // Count how many times a node is referenced in the algorithm
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

  // Stack used to sort the computational graph
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
    casadi_error("Unrecongnized topological_sorting: " << getOption("topological_sorting"));
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

  // Evaluate with inplace operations (experimental)
  evaluate_inplace_ = getOption("inplace");
  vector<int> refcount_copy;
  if(evaluate_inplace_){
    casadi_assert_message(live_variables,"\"evaluate_inplace\" requires \"live_variables\"");
  }

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

    // Make a copy of refcount for later
    if(evaluate_inplace_){
      refcount_copy = refcount;
    }

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

  // Allocate a vector containing expression corresponding to each binary operation
  binops_.resize(bnodes.size());
  for(int i=0; i<bnodes.size(); ++i){
    binops_[i] = SX::create(bnodes[i]);
  }
  bnodes.clear();
  
  // Number of inplace functions
  int num_inplace=0;
  
  // Index currently assigned to a node
  vector<int> curr_assign;
  if(evaluate_inplace_){
    curr_assign.resize(worksize_,-1);
  }
  
  // Add the binary operations
  algorithm_.clear();
  algorithm_.resize(binops_.size());
  pder_.resize(algorithm_.size());
  vector<AlgEl>::iterator it=algorithm_.begin();
  for(int i=0; i<algorithm_.size(); ++i, ++it){

    // Save the node index
    it->ind = place[binops_[i]->temp];
    
    // Save the indices of the children
    it->ch[0] = place[binops_[i]->dep(0).get()->temp];
    it->ch[1] = place[binops_[i]->dep(1).get()->temp];
    
    // Make sure that the first argument is equal to the index, if possible
//     if(it->ch[1]==it->ind && 
//       casadi_math<double>::isCommutative(it->op) && 
//       ndeps==2){
//       std::swap(it->ch[0],it->ch[1]);
//     }

    // Operation
    it->op = binops_[i]->getOp();

    // Replace inplace operations
    if(evaluate_inplace_){
      // Check if the operation can be performed inplace
      if(it->ch[0]==it->ind){
        int ip;
        switch(it->op){
          case ADD: ip=1; break;
          case SUB: ip=2; break;
          case MUL: ip=3; break;
          case DIV: ip=4; break;
          default:  ip=0; break;
        }
        
        // if the operation can be performed inplace
        if(ip!=0){
          
          // Reference to the child node
          SX& c = binops_[i]->dep(1);
          
          // If the dependent variable is a binary variable
          bool dep_is_binary = c.isBinary();
          if(dep_is_binary){

            // Get the indices of the children of the child
            int c_temp = c.get()->temp;
            int c_temp0 = c->dep(0).get()->temp;
            int c_temp1 = c->dep(1).get()->temp;
            
            // Index currently assigned to the storage location
            int c_curr = curr_assign[place[c_temp]];
            int c_curr0 = curr_assign[place[c_temp0]];
            int c_curr1 = curr_assign[place[c_temp1]];
            
            // Number of dependencies
            int c_ndeps = casadi_math<double>::ndeps(c.getOp());
            
            // Check if the node is used exactly one time
            bool ok_replace = refcount_copy.at(c_temp)==1;
            
            // Check if the values of the children are still available in the work vector or was assigned during the child operation
            ok_replace = ok_replace && (c_curr0<0 || c_curr0==c_temp0 || c_curr0==c_curr);
            ok_replace = ok_replace && (c_ndeps<2 ||  c_curr1<1 || c_curr1==c_temp1 || c_curr1==c_curr);
            
            if(ok_replace){
              
              // Mark the node negative so that we can remove it later
              refcount_copy.at(c.get()->temp) = -1;
              
              // Save the indices of the children
              it->ch[0] = place[c_temp0];
              it->ch[1] = place[c_temp1];

              // Replace operation with an inplace operation
              it->op = ip*NUM_BUILT_IN_OPS + c->getOp();
              
              // Increase the inplace counter
              num_inplace++;
            }
          }
        }
      }
      // Save index
      curr_assign[it->ind] = binops_[i]->temp;
    }
  }
  
  if(num_inplace>0){
    if(verbose()){
      cout << "SXFunctionInternal::init Located " << num_inplace << " inplace operators" << endl;
    }
    
    vector<AlgEl> algorithm_new;
    algorithm_new.reserve(binops_.size()-num_inplace);
    
    for(int i=0; i<binops_.size(); ++i, ++it){
      
      // Check if marked
      bool marked = refcount_copy[binops_[i].get()->temp]<0;
      
      // Add if not marked
      if(marked){
        num_inplace--;
      } else {
        algorithm_new.push_back(algorithm_[i]);
      }
    }
    
    // Consistency check
    casadi_assert(num_inplace==0);
    
    // Save modified algorithm
    algorithm_new.swap(algorithm_);
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

  // Allocate memory for directional derivatives
  SXFunctionInternal::updateNumSens(false);
  
  // Get the full Jacobian already now
  if(jac_for_sens_){
    getFullJacobian();
  }
  
  // Print
  if(verbose()){
    cout << "SXFunctionInternal::init Initialized " << getOption("name") << " (" << algorithm_.size() << " elementary operations)" << endl;
  }
}

void SXFunctionInternal::updateNumSens(bool recursive){
  // Call the base class if needed
  if(recursive) XFunctionInternal::updateNumSens(recursive);
  
  // Allocate a working array
  if(nfdir_>0 || nadir_>0){
    dwork_.resize(worksize_,numeric_limits<double>::quiet_NaN());
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
    SXMatrix h = gfcn.jac(0,0,false,true);
    
    if(false){ // Does not appear to help
      // Calculate the transpose of the sparsity pattern
      vector<int> mapping;
      CRSSparsity h_sp_trans = h.sparsity().transpose(mapping);
      
      // Make sure that the Hesssian is indeed symmetric
      casadi_assert(h_sp_trans == h.sparsity());
      
      // Access sparsity pattern
      int h_nrow=h.sparsity().size1();
      const vector<int> h_rowind =h.sparsity().rowind();
      const vector<int> h_col =h.sparsity().col();
      vector<SX>& h_data = h.data();
      
      // Loop over the rows of the hessian
      for(int i=0; i<h_nrow; ++i){
        // Loop over the strictly lower triangular part of the matrix
        for(int el=h_rowind[i]; el<h_rowind[i+1] && h_col[el]<i; ++el){
          // Mirror upper triangular part
          // h_data[el] = h_data[mapping[el]];
          h_data[mapping[el]] = h_data[el];
        }
      }
    }
    
    // Create the hessian function
    hfcn = SXFunction(inputv_,h);
  }
  
  // Calculate jacobian of gradient
  if(verbose()) cout << "SXFunctionInternal::hessian: calculating Hessian done" << endl;

  // Return jacobian of the gradient
  return hfcn;
}

void SXFunctionInternal::evalSX(const std::vector<SXMatrix>& input, std::vector<SXMatrix>& output, 
                                const std::vector<std::vector<SXMatrix> >& fwdSeed, std::vector<std::vector<SXMatrix> >& fwdSens, 
                                const std::vector<std::vector<SXMatrix> >& adjSeed, std::vector<std::vector<SXMatrix> >& adjSens,
                                bool output_given, bool eliminate_constants){
  
  if(verbose()) cout << "SXFunctionInternal::eval begin" << endl;
  
  casadi_assert_message(inputv_.size() == input.size(),"SXFunctionInternal::eval: wrong number of inputs." << std::endl << "Expecting " << inputv_.size() << " inputs, but got " << input.size() << " instead.");
  for(int i=0; i<input.size(); ++i){
    casadi_assert_message(input[i].sparsity()==inputv_[i].sparsity(), "SXFunctionInternal::eval: argument sparsity inconsistent");
  }
  
  casadi_assert_message(outputv_.size() == output.size(),"SXFunctionInternal::eval: wrong number of outputs." << std::endl << "Expecting " << outputv_.size() << " inputs, but got " << output.size() << " instead.");
  for(int i=0; i<output.size(); ++i){
    casadi_assert_message(output[i].sparsity()==outputv_[i].sparsity(), "SXFunctionInternal::evals: result sparsity inconsistent");
  }
  
  if(output_given){
    // Iterator to the binary operations
    vector<SX>::const_iterator f_it=binops_.begin();
    
    // Assign to the symbolic work vector
    for(vector<AlgEl>::const_iterator it = algorithm_.begin(); it!=algorithm_.end(); ++it, ++f_it){
      swork_[it->ind] = *f_it;
    }
    
  } else {
    // Copy the function arguments to the work vector
    for(int ind=0; ind<input.size(); ++ind){
      const vector<SX>& idata = input[ind].data();
      for(int i=0; i<input_ind_[ind].size(); ++i){
        swork_[input_ind_[ind][i]] = idata[i];
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
        casadi_math<double>::fun(it->op,x.getValue(),y.getValue(),temp);
        swork_[it->ind] = temp;
      } else {
        casadi_math<SX>::fun(it->op,x,y,swork_[it->ind]);
      }
    }

    // Get the results
    for(int ind=0; ind<output.size(); ++ind){
      vector<SX>& odata = output[ind].data();
      for(int i=0; i<output_ind_[ind].size(); ++i){
        odata[i] = swork_[output_ind_[ind][i]];
      }
    }
  }

  // Get the number of forward and adjoint sweeps
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();
  
/*  std::cout << "nfwd = " <<  nfwd << std::endl;
  std::cout << "nadj = " <<  nadj << std::endl;*/
  

  
  
  // Quick return if no sensitivities
  if(nfwd==0 && nadj==0) return;
  
  // Calculate the partial derivatives
  vector<SX> der1, der2;
  der1.reserve(algorithm_.size());
  der2.reserve(algorithm_.size());
  SX tmp[2];
  for(vector<AlgEl>::const_iterator it = algorithm_.begin(); it!=algorithm_.end(); ++it){
      const SX& f = swork_[it->ind];
      const SX& x = f->dep(0);
      const SX& y = f->dep(1);
      casadi_math<SX>::der(it->op,x,y,f,tmp);
      if(!x->isConstant())  der1.push_back(tmp[0]);
      else                  der1.push_back(0);

      if(!y->isConstant())  der2.push_back(tmp[1]);
      else                  der2.push_back(0);
  }

  // Gradient (this is also the working array)
  vector<SX> g(swork_.size(),casadi_limits<SX>::zero);
  if(verbose())   cout << "SXFunctionInternal::evalNew gradient working array size " << swork_.size() << endl;

  // Carry out the forward sweeps
  for(int dir=0; dir<nfwd; ++dir){
    
    // Remove Seeds
    fill(g.begin(),g.end(),0);
    
    // Pass the forward seeds
    for(int iind=0; iind<input.size(); ++iind){
      // Assert sparsity
      casadi_assert_message(input[iind].sparsity()==fwdSeed[dir][iind].sparsity(), "SXFunctionInternal::eval: sparsity inconsistent");

      // Copy seeds to work vector
      const vector<SX>& fdata = fwdSeed[dir][iind].data();
      for(int k=0; k<fdata.size(); ++k){
        g[input_ind_[iind][k]] = fdata[k];
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
    
    // Get the forward sensitivities
    for(int oind=0; oind<output.size(); ++oind){
      // Assert sparsity
      casadi_assert_message(output[oind].sparsity()==fwdSens[dir][oind].sparsity(), "SXFunctionInternal::eval: sparsity inconsistent");
      
      // Copy sens from work vector
      vector<SX>& fdata = fwdSens[dir][oind].data();
      for(int k=0; k<fdata.size(); ++k){
        fdata[k] = g[output_ind_[oind][k]];
      }
    }
  }
  
  // And now the adjoint sweeps
  for(int dir=0; dir<nadj; ++dir){
    
    // Remove Seeds
    fill(g.begin(),g.end(),0);
  
    // Pass the adjoint seeds
    for(int oind=0; oind<output.size(); ++oind){
      // Assert sparsity
      casadi_assert_message(output[oind].sparsity()==adjSeed[dir][oind].sparsity(), "SXFunctionInternal::eval: sparsity inconsistent");
      
      // Copy sens from work vector
      const vector<SX>& adata = adjSeed[dir][oind].data();
      for(int k=0; k<adata.size(); ++k){
        g[output_ind_[oind][k]] = adata[k];
      }
    }
    
    // backward sweep
    for(int k = algorithm_.size()-1; k>=0; --k){
      const AlgEl& ae = algorithm_[k];

      // Get the seed
      SX seed = g[ae.ind];
      
      // Clear the seed
      g[ae.ind] = 0;

      // Propagate if the seed is not zero
      if(!seed->isZero()){
        if(!der1[k]->isZero())
          g[ae.ch[0]] += der1[k] * seed;
        if(!der2[k]->isZero())
          g[ae.ch[1]] += der2[k] * seed;
      }
    }
    
    // Get the adjoint sensitivities
    for(int iind=0; iind<input.size(); ++iind){
      // Assert sparsity
      casadi_assert_message(input[iind].sparsity()==adjSens[dir][iind].sparsity(), "SXFunctionInternal::eval: sparsity inconsistent");
      
      // Copy sens from work vector
      vector<SX>& adata = adjSens[dir][iind].data();
      for(int k=0; k<adata.size(); ++k){
        adata[k] = g[input_ind_[iind][k]];
      }
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
  
  vector<SXMatrix> inputv_v(inputv_);
  
  // Return function
  return SXFunction(inputv_v,jac_out);
}

void SXFunctionInternal::spProp(bool fwd){
  if(fwd){
    for(std::vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it){
      iwork_[it->ind] = iwork_[it->ch[0]] | iwork_[it->ch[1]];
    }
  } else {
    // The following is commented out due to bug(?) in Mac using old gcc
    // for(vector<AlgEl>::reverse_iterator it=algorithm_.rbegin(); it!=algorithm_.rend(); ++it)
    for(int i=algorithm_.size()-1; i>=0; --i){ // workaround
      AlgEl *it = &algorithm_[i];
      
      // Get the seed
      bvec_t seed = iwork_[it->ind];
      
      // Clear the seed
      iwork_[it->ind] = 0;
      
      // Propagate seeds
      iwork_[it->ch[0]] |= seed;
      iwork_[it->ch[1]] |= seed;
    }
  }
}

void SXFunctionInternal::spReset(int iind, int oind){
  // Make sure that dwork_, which we will now use, has been allocated
  if(dwork_.size() < worksize_) dwork_.resize(worksize_);
  
  // We need a work array containing unsigned long rather than doubles. Since the two datatypes have the same size (64 bits)
  // we can save overhead by reusing the double array
  iwork_ = get_bvec_t(dwork_);
  fill_n(iwork_,dwork_.size(),0);
}

} // namespace CasADi

