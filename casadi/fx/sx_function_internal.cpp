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
#include "../casadi_types.hpp"
#include "../matrix/crs_sparsity_internal.hpp"

#ifdef WITH_LLVM
#include "llvm/DerivedTypes.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/ExecutionEngine/JIT.h"
#include "llvm/LLVMContext.h"
#include "llvm/Module.h"
#include "llvm/PassManager.h"
#include "llvm/Analysis/Verifier.h"
#include "llvm/Target/TargetData.h"
#include "llvm/Target/TargetSelect.h"
#include "llvm/Transforms/Scalar.h"
#include "llvm/Support/IRBuilder.h"

llvm::IRBuilder<> builder(llvm::getGlobalContext());
#endif // WITH_LLVM

namespace CasADi{

using namespace std;


SXFunctionInternal::SXFunctionInternal(const vector<Matrix<SX> >& inputv, const vector<Matrix<SX> >& outputv) : 
  XFunctionInternal<SXFunctionInternal,Matrix<SX>,SXNode>(inputv,outputv) {
  setOption("name","unnamed_sx_function");
  addOption("live_variables",OT_BOOLEAN,false,"Reuse variables in the work vector");
  addOption("inplace",OT_BOOLEAN,false,"Evaluate with inplace operations (experimental)");
  addOption("just_in_time",OT_BOOLEAN,false,"Just-in-time compilation for numeric evaluation (experimental)");

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
  
  #ifdef WITH_LLVM
  if(just_in_time_ && nfdir==0 && nadir==0){
    // Evaluate the jitted function
    jitfcn_(getPtr(input_ref_),getPtr(output_ref_));
    return;
  }
  #endif // WITH_LLVM
  
  
  // Copy the function arguments to the work vector
  for(int ind=0; ind<getNumInputs(); ++ind){
    const Matrix<double> &arg = input(ind);
    for(int i=0; i<arg.size(); ++i){
      work_[input_ind_[ind][i]] = arg.data()[i];
    }
  }

  // Evaluate the algorithm
  if(nfdir==0 && nadir==0){ // without taping
    for(vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it){
      // NOTE: This is equivalent to casadi_math<double>::inplacefun(it->op,x,y,f);
      // but forces the function to be inlined, which is important for speed here
      switch(it->op){
        // Start by adding all of the built operations
        CASADI_MATH_FUN_ALL_BUILTIN(work_[it->arg.i[0]],work_[it->arg.i[1]],work_[it->res])
        
        // Now all other operations
	
	// Constant
        case OP_CONST: work_[it->res] = it->arg.d; break;
	
	// Get input
// 	case OP_INPUT: work_[it->res] = input(it->arg.i[0]).at(it->arg.i[1]); break;
// 	case OP_OUTPUT: input(it->arg.i[0]).at(it->arg.i[1]) = work_[it->res]; break;
      }
    }
  } else { // with taping
    
    vector<AlgEl>::iterator it = algorithm_.begin();
    vector<AlgElData>::iterator it1 = pder_.begin();
    for(; it!=algorithm_.end(); ++it, ++it1){
      // NOTE: This is equivalent to casadi_math<double>::derF(it->op,x,y,f,it1->d);
      // but forces the function to be inlined, which is important for speed here
      CASADI_MATH_DERF(double,it->op,work_[it->arg.i[0]],work_[it->arg.i[1]],work_[it->res],it1->d)
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
      dwork_[it->res] = it2->d[0] * dwork_[it->arg.i[0]] + it2->d[1] * dwork_[it->arg.i[1]];
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
      double seed = dwork_[it->res];
      dwork_[it->res] = 0;
      dwork_[it->ch[0]] += it2->d[0] * seed;
      dwork_[it->ch[1]] += it2->d[1] * seed;
    }
    #else
    for(int i=algorithm_.size()-1; i>=0; --i){
      const AlgEl& ae = algorithm_[i];
      const AlgElData& aed = pder_[i];
      
      // copy the seed and clear the cache entry
      double seed = dwork_[ae.res];
      dwork_[ae.res] = 0;

      dwork_[ae.arg.i[0]] += aed.d[0] * seed;
      dwork_[ae.arg.i[1]] += aed.d[1] * seed;
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

SXMatrix SXFunctionInternal::jac(int iind, int oind, bool compact, bool symmetric){
  return jacGen(iind,oind,compact,symmetric);
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

 // If JIT, dump LLVM IR
  #ifdef WITH_LLVM
  if(just_in_time_){
    jit_module_->dump();
    return;
  }
  #endif // WITH_LLVM
  
  // Normal, interpreted output
  for(vector<AlgEl>::const_iterator it = algorithm_.begin(); it!=algorithm_.end(); ++it){
    int op = it->op;
    int ip = op/NUM_BUILT_IN_OPS;
    op -= NUM_BUILT_IN_OPS*ip;
    
    stringstream s,s0,s1;
    s << "a" << it->res;

    int i0 = it->arg.i[0], i1 = it->arg.i[1];
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
      stream << "a" << ae.arg.i[c];
    }
  }
  casadi_math<double>::printPost(ae.op,stream);
}

void SXFunctionInternal::generateCode(const string& src_name){
  assertInit();

  // Make sure that there are no free variables
  if (!free_vars_.empty()) {
    casadi_error("Code generation is not possible since variables " << free_vars_ << " are free.");
  }
  
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
   
    if(!declared[it->res]){
      cfile << "d ";
      declared[it->res]=true;
    }
    cfile << "a" << it->res << "=";
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
        cfile << "a" << algorithm_[i-1].res;
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
  XFunctionInternal<SXFunctionInternal,Matrix<SX>,SXNode>::init();

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
    it->res = place[binops_[i]->temp];
    
    // Save the indices of the children
    it->arg.i[0] = place[binops_[i]->dep(0).get()->temp];
    it->arg.i[1] = place[binops_[i]->dep(1).get()->temp];
    
    // Make sure that the first argument is equal to the index, if possible
//     if(it->ch[1]==it->res && 
//       casadi_math<double>::isCommutative(it->op) && 
//       ndeps==2){
//       std::swap(it->ch[0],it->ch[1]);
//     }

    // Operation
    it->op = binops_[i]->getOp();

    // Replace inplace operations
    if(evaluate_inplace_){
      // Check if the operation can be performed inplace
      if(it->arg.i[0]==it->res){
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
          
	  // Number of dependencies
	  int c_ndep = c->ndep();
	  
          // If the dependent variable is a unary variable
          if(c_ndep==1){

            int c_temp = c.get()->temp;
            int c_curr = curr_assign[place[c_temp]];

            // Check if the node is used exactly one time
            bool ok_replace = refcount_copy.at(c_temp)==1;
	    	    
	    int c_temp0 = c->dep(0).get()->temp;
            int c_curr0 = curr_assign[place[c_temp0]];
            	                            
            // Check if the values of the children are still available in the work vector or was assigned during the child operation
            ok_replace = ok_replace && (c_curr0<0 || c_curr0==c_temp0 || c_curr0==c_curr);
            
            if(ok_replace){
              
              // Mark the node negative so that we can remove it later
              refcount_copy.at(c.get()->temp) = -1;
              
              // Save the indices of the children
              it->arg.i[0] = place[c_temp0];

              // Replace operation with an inplace operation
              it->op = ip*NUM_BUILT_IN_OPS + c->getOp();
              
              // Increase the inplace counter
              num_inplace++;
            }
          }
          	  
          // If the dependent variable is a binary variable
          if(c_ndep==2){

            int c_temp = c.get()->temp;
            int c_curr = curr_assign[place[c_temp]];

            // Check if the node is used exactly one time
            bool ok_replace = refcount_copy.at(c_temp)==1;
	    	    
	    int c_temp0 = c->dep(0).get()->temp;
            int c_curr0 = curr_assign[place[c_temp0]];
            	    
	    int c_temp1 = c->dep(1).get()->temp;
            int c_curr1 = curr_assign[place[c_temp1]];
                        
            // Check if the values of the children are still available in the work vector or was assigned during the child operation
            ok_replace = ok_replace && (c_curr0<0 || c_curr0==c_temp0 || c_curr0==c_curr);
            ok_replace = ok_replace && (c_curr1<1 || c_curr1==c_temp1 || c_curr1==c_curr);
            
            if(ok_replace){
              
              // Mark the node negative so that we can remove it later
              refcount_copy.at(c.get()->temp) = -1;
              
              // Save the indices of the children
              it->arg.i[0] = place[c_temp0];
              it->arg.i[1] = place[c_temp1];

              // Replace operation with an inplace operation
              it->op = ip*NUM_BUILT_IN_OPS + c->getOp();
              
              // Increase the inplace counter
              num_inplace++;
            }
          }
        }
      }
      // Save index
      curr_assign[it->res] = binops_[i]->temp;
    }
  }
  
  if(num_inplace>0){
    if(verbose()){
      cout << "SXFunctionInternal::init Replacing " << num_inplace << " inplace operators" << endl;
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
  
  // Initialize just-in-time compilation
  just_in_time_ = getOption("just_in_time");
  if(just_in_time_){
    #ifdef WITH_LLVM
    llvm::InitializeNativeTarget();

    // Function name
    stringstream ss;
    ss << "SXFunction: " << this;
    
    // Make the module, which holds all the code.
    jit_module_ = new llvm::Module(ss.str(), llvm::getGlobalContext());

    // Create the JIT.  This takes ownership of the module.
    std::string ErrStr;
    llvm::ExecutionEngine *TheExecutionEngine = llvm::EngineBuilder(jit_module_).setErrorStr(&ErrStr).create();
    casadi_assert(TheExecutionEngine!=0);
    llvm::FunctionPassManager OurFPM(jit_module_);

    // Set up the optimizer pipeline.  Start with registering info about how the
    // target lays out data structures.
    OurFPM.add(new llvm::TargetData(*TheExecutionEngine->getTargetData()));
    
    // Do simple "peephole" optimizations and bit-twiddling optzns.
    OurFPM.add(llvm::createInstructionCombiningPass());
    
    // Reassociate expressions.
    OurFPM.add(llvm::createReassociatePass());
    
    // Eliminate Common SubExpressions.
    OurFPM.add(llvm::createGVNPass());
    
    // Simplify the control flow graph (deleting unreachable blocks, etc).
    OurFPM.add(llvm::createCFGSimplificationPass());
    OurFPM.doInitialization();

    // Single argument
    std::vector<const llvm::Type*> unaryArg(1,llvm::Type::getDoubleTy(llvm::getGlobalContext()));

    // Two arguments
    std::vector<const llvm::Type*> binaryArg(2,llvm::Type::getDoubleTy(llvm::getGlobalContext()));
    
    // Unary operation
    llvm::FunctionType *unaryFun = llvm::FunctionType::get(llvm::Type::getDoubleTy(llvm::getGlobalContext()),unaryArg, false);

    // Binary operation
    llvm::FunctionType *binaryFun = llvm::FunctionType::get(llvm::Type::getDoubleTy(llvm::getGlobalContext()),binaryArg, false);

    // Declare all the CasADi built-in functions
    vector<llvm::Function*> builtins(NUM_BUILT_IN_OPS,0);
    builtins[POW] = builtins[CONSTPOW] = llvm::Function::Create(binaryFun, llvm::Function::ExternalLinkage, "pow", jit_module_);
    builtins[SQRT] = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "sqrt", jit_module_);
    builtins[SIN] = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "sin", jit_module_);
    builtins[COS] = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "cos", jit_module_);
    builtins[TAN] = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "tan", jit_module_);
    builtins[ASIN] = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "asin", jit_module_);
    builtins[ACOS] = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "acos", jit_module_);
    builtins[ATAN] = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "atan", jit_module_);
    builtins[FLOOR] = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "floor", jit_module_);
    builtins[CEIL] = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "ceil", jit_module_);
    builtins[FMIN] = llvm::Function::Create(binaryFun, llvm::Function::ExternalLinkage, "fmin", jit_module_);
    builtins[FMAX] = llvm::Function::Create(binaryFun, llvm::Function::ExternalLinkage, "fmax", jit_module_);
    builtins[SINH] = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "sinh", jit_module_);
    builtins[COSH] = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "cosh", jit_module_);
    builtins[TANH] = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "tanh", jit_module_);

    // Void type
    const llvm::Type* void_t = llvm::Type::getVoidTy(llvm::getGlobalContext());

    // Double type
    const llvm::Type* double_t = llvm::Type::getDoubleTy(llvm::getGlobalContext());
    
    // Double pointer type
    const llvm::Type* double_ptr_t = llvm::Type::getDoublePtrTy(llvm::getGlobalContext());

    // Double pointer pointer type
    const llvm::Type* double_ptr_ptr_t = llvm::PointerType::getUnqual(double_ptr_t);

    // A normal 32-bit integer
    const llvm::IntegerType *int32Ty = llvm::IntegerType::get(llvm::getGlobalContext(), 32);
    
    // Two arguments in and two references
    std::vector<const llvm::Type*> genArg(2);
    genArg[0] = double_ptr_ptr_t;
    genArg[1] = double_ptr_ptr_t;
    
    // More generic operation, return by reference
    llvm::FunctionType *genFun = llvm::FunctionType::get(void_t,genArg, false);

    // Declare my function
    jit_function_ = llvm::Function::Create(genFun, llvm::Function::ExternalLinkage, ss.str(), jit_module_);

    // Create a new basic block to start insertion into.
    llvm::BasicBlock *BB = llvm::BasicBlock::Create(llvm::getGlobalContext(), "entry", jit_function_);
    builder.SetInsertPoint(BB);

    // Set names for all arguments.
    llvm::Function::arg_iterator AI = jit_function_->arg_begin();
    AI->setName("x");  llvm::Value *x_ptr = AI++;
    AI->setName("r");  llvm::Value *r_ptr = AI++;

    // Allocate work vector
    std::vector<llvm::Value*> jwork(work_.size());

    // Add all constants to the work vector
    for(int k=0; k<jwork.size(); ++k){
      jwork[k] = llvm::ConstantFP::get(llvm::getGlobalContext(), llvm::APFloat(work_[k]));
    }
    
    // Copy the arguments to the work vector
    for(int ind=0; ind<input_ind_.size(); ++ind){
      // Get the vector-valued argument first
      llvm::Value *ind_v = llvm::ConstantInt::get(int32Ty, ind);
      llvm::Value* x_arg_v = builder.CreateLoad(builder.CreateGEP(x_ptr,ind_v));

      // Copy the scalar-valued arguments
      for(int k=0; k<input_ind_[ind].size(); ++k){
	llvm::Value *k_v = llvm::ConstantInt::get(int32Ty, k);
	jwork[input_ind_[ind][k]] = builder.CreateLoad(builder.CreateGEP(x_arg_v,k_v));
      }
    }
    
    // Build up the LLVM expression graphs
    for(vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it){
      
      // Argument of the operation
      std::vector<llvm::Value*> oarg(casadi_math<double>::ndeps(it->op));
      for(int d=0; d<oarg.size(); ++d){
	oarg[d] = jwork[it->arg.i[d]];
      }
      
      // Result
      llvm::Value* res = 0;
      
      switch(it->op){
	case ADD: res = builder.CreateFAdd(oarg[0],oarg[1]); break;
	case SUB: res = builder.CreateFSub(oarg[0],oarg[1]); break;
	case MUL: res = builder.CreateFMul(oarg[0],oarg[1]); break;
	case DIV: res = builder.CreateFDiv(oarg[0],oarg[1]); break;
	case NEG: res = builder.CreateFNeg(oarg[0]);         break;

	default:
	  casadi_assert_message(builtins[it->op]!=0, "No way to treat: " << it->op);
	  res = builder.CreateCall(builtins[it->op], oarg.begin(), oarg.end());
      }
      
      // Save to work vector
      jwork[it->res] = res;
    }
    
    // Store the results from the work vector
    for(int ind=0; ind<output_ind_.size(); ++ind){
      // Get the vector-valued argument first
      llvm::Value *ind_v = llvm::ConstantInt::get(int32Ty, ind);
      llvm::Value* r_arg_v = builder.CreateLoad(builder.CreateGEP(r_ptr,ind_v));

      // Store all scalar-valued arguments
      for(int k=0; k<output_ind_[ind].size(); ++k){
	llvm::Value *k_v = llvm::ConstantInt::get(int32Ty, k);
	builder.CreateStore(jwork[output_ind_[ind][k]],builder.CreateGEP(r_arg_v,k_v));
      }
    }
    
    // Finish off the function.
    builder.CreateRetVoid();

    // Validate the generated code, checking for consistency.
    verifyFunction(*jit_function_);

    // Optimize the function.
    OurFPM.run(*jit_function_);

    // JIT the function
    jitfcn_ = evaluateFcn(intptr_t(TheExecutionEngine->getPointerToFunction(jit_function_)));

    // Allocate references to input nonzeros
    input_ref_.resize(getNumInputs());
    for(int ind=0; ind<input_ref_.size(); ++ind){
      input_ref_[ind] = getPtr(input(ind).data());
    }
        
    // Allocate references to output nonzeros
    output_ref_.resize(getNumOutputs());
    for(int ind=0; ind<output_ref_.size(); ++ind){
      output_ref_[ind] = getPtr(output(ind).data());
    }
    
    #else // WITH_LLVM
    casadi_error("Option \"just_in_time\" true requires CasADi to have been compiled with WITH_LLVM=ON");
    #endif //WITH_LLVM
  }
  
  // Print
  if(verbose()){
    cout << "SXFunctionInternal::init Initialized " << getOption("name") << " (" << algorithm_.size() << " elementary operations)" << endl;
  }
}

void SXFunctionInternal::updateNumSens(bool recursive){
  // Call the base class if needed
  if(recursive) XFunctionInternal<SXFunctionInternal,Matrix<SX>,SXNode>::updateNumSens(recursive);
  
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
    casadi_assert_message(input[i].sparsity()==inputv_[i].sparsity() || (input[i].empty() && inputv_[i].empty()), "SXFunctionInternal::eval: sparsity of argument " << i << " inconsistent");
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
      swork_[it->res] = *f_it;
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
      SX x = swork_[it->arg.i[0]];
      SX y = swork_[it->arg.i[1]];
      if(eliminate_constants && x.isConstant() && y.isConstant()){
        // Check if both arguments are constants
        double temp;
        casadi_math<double>::fun(it->op,x.getValue(),y.getValue(),temp);
        swork_[it->res] = temp;
      } else {
        casadi_math<SX>::fun(it->op,x,y,swork_[it->res]);
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
    
  // Quick return if no sensitivities
  if(nfwd==0 && nadj==0) return;
  
  // Calculate the partial derivatives
  vector<SX> der1, der2;
  der1.reserve(algorithm_.size());
  der2.reserve(algorithm_.size());
  SX tmp[2];
  for(vector<AlgEl>::const_iterator it = algorithm_.begin(); it!=algorithm_.end(); ++it){
      const SX& f = swork_[it->res];
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
      casadi_assert_message(input[iind].sparsity()==fwdSeed[dir][iind].sparsity(), "SXFunctionInternal::eval: sparsity inconsistent for input " << iind << " and forward direction " << dir << ".");

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
      if(!g[ae.arg.i[0]]->isZero() &&  !der1[k]->isZero()){
        g[ae.res] += der1[k] * g[ae.arg.i[0]];
      }
      
      // Second argument
      if(!g[ae.arg.i[1]]->isZero() &&  !der2[k]->isZero()){
        g[ae.res] += der2[k] * g[ae.arg.i[1]];
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
      SX seed = g[ae.res];
      
      // Clear the seed
      g[ae.res] = 0;

      // Propagate if the seed is not zero
      if(!seed->isZero()){
        if(!der1[k]->isZero())
          g[ae.arg.i[0]] += der1[k] * seed;
        if(!der2[k]->isZero())
          g[ae.arg.i[1]] += der2[k] * seed;
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
      jac_out[el] = jac(iind,oind);
    }
  }
  
  vector<SXMatrix> inputv_v(inputv_);
  
  // Return function
  return SXFunction(inputv_v,jac_out);
}

void SXFunctionInternal::spInit(bool fwd){
  // Make sure that dwork_, which we will now use, has been allocated
  if(dwork_.size() < worksize_) dwork_.resize(worksize_);
  
  // We need a work array containing unsigned long rather than doubles. Since the two datatypes have the same size (64 bits)
  // we can save overhead by reusing the double array
  bvec_t *iwork = get_bvec_t(dwork_);
  fill_n(iwork,dwork_.size(),0);
}

void SXFunctionInternal::spEvaluate(bool fwd){
  // Get work array
  bvec_t *iwork = get_bvec_t(dwork_);

  if(fwd){ // Forward propagation
    
    // Pass input seeds
    for(int ind=0; ind<input_ind_.size(); ++ind){
      bvec_t* swork = get_bvec_t(input(ind).data());
      for(int k=0; k<input_ind_[ind].size(); ++k){
	iwork[input_ind_[ind][k]] = swork[k];
      }
    }
    
    // Propagate sparsity forward
    for(std::vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it){
      iwork[it->res] = iwork[it->arg.i[0]] | iwork[it->arg.i[1]];
    }
    
    // Get the output seeds
    for(int ind=0; ind<output_ind_.size(); ++ind){
      bvec_t* swork = get_bvec_t(output(ind).data());
      for(int k=0; k<output_ind_[ind].size(); ++k){
	swork[k] = iwork[output_ind_[ind][k]];
      }
    }
    
  } else { // Backward propagation

    // Pass output seeds
    for(int ind=0; ind<output_ind_.size(); ++ind){
      bvec_t* swork = get_bvec_t(output(ind).data());
      for(int k=0; k<output_ind_[ind].size(); ++k){
	iwork[output_ind_[ind][k]] = swork[k];
      }
    }

    // Propagate sparsity backward
    for(vector<AlgEl>::reverse_iterator it=algorithm_.rbegin(); it!=algorithm_.rend(); ++it){
      
      // Get the seed
      bvec_t seed = iwork[it->res];
      
      // Clear the seed
      iwork[it->res] = 0;
      
      // Propagate seeds
      iwork[it->arg.i[0]] |= seed;
      iwork[it->arg.i[1]] |= seed;
    }
    
    // Get the input seeds and clear it from the work vector
    for(int ind=0; ind<input_ind_.size(); ++ind){
      bvec_t* swork = get_bvec_t(input(ind).data());
      for(int k=0; k<input_ind_[ind].size(); ++k){
	swork[k] |= iwork[input_ind_[ind][k]];
	iwork[input_ind_[ind][k]] = 0;
      }
    }
  }
}

FX SXFunctionInternal::getDerivative(int nfwd, int nadj){
  // Forward seeds
  vector<vector<SXMatrix> > fseed(nfwd,inputv_);
  for(int dir=0; dir<nfwd; ++dir){
    // Replace symbolic inputs
    for(vector<SXMatrix>::iterator j=fseed[dir].begin(); j!=fseed[dir].end(); ++j){
      for(vector<SX>::iterator i=j->begin(); i!=j->end(); ++i){
	// Name of the forward seed
	std::stringstream ss;
	ss << "f";
	if(nfwd>1) ss << dir;
	ss << "_" << *i;
	
	// Save to matrix
	*i = SX(ss.str());
      }
    }
  }
  
  // Adjoint seeds
  vector<vector<SXMatrix> > aseed(nadj,outputv_);
  for(int dir=0; dir<nadj; ++dir){
    // Replace symbolic inputs
    int oind=0;
    for(vector<SXMatrix>::iterator j=aseed[dir].begin(); j!=aseed[dir].end(); ++j, ++oind){
      int k=0;
      for(vector<SX>::iterator i=j->begin(); i!=j->end(); ++i, ++k){
	// Name of the adjoint seed
	std::stringstream ss;
	ss << "a";
	if(nadj>1) ss << dir << "_";
	ss << oind << "_" << k;

	// Save to matrix
	*i = SX(ss.str());
      }
    }
  }
  
  // Evaluate symbolically
  vector<vector<SXMatrix> > fsens(nfwd,outputv_), asens(nadj,inputv_);
  evalSX(inputv_,outputv_,fseed,fsens,aseed,asens,true,false);

  // All inputs of the return function
  vector<SXMatrix> ret_in;
  ret_in.reserve(inputv_.size()*(1+nfwd) + outputv_.size()*nadj);
  ret_in.insert(ret_in.end(),inputv_.begin(),inputv_.end());
  for(int dir=0; dir<nfwd; ++dir)
    ret_in.insert(ret_in.end(),fseed[dir].begin(),fseed[dir].end());
  for(int dir=0; dir<nadj; ++dir)
    ret_in.insert(ret_in.end(),aseed[dir].begin(),aseed[dir].end());
  
  // All outputs of the return function
  vector<SXMatrix> ret_out;
  ret_out.reserve(outputv_.size()*(1+nfwd) + inputv_.size()*nadj);
  ret_out.insert(ret_out.end(),outputv_.begin(),outputv_.end());
  for(int dir=0; dir<nfwd; ++dir)
    ret_out.insert(ret_out.end(),fsens[dir].begin(),fsens[dir].end());
  for(int dir=0; dir<nadj; ++dir)
    ret_out.insert(ret_out.end(),asens[dir].begin(),asens[dir].end());

  // Assemble function and return
  SXFunction ret(ret_in,ret_out);
  return ret;
}

} // namespace CasADi

