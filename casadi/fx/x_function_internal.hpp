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

#ifndef X_FUNCTION_INTERNAL_HPP
#define X_FUNCTION_INTERNAL_HPP

#include <map>
#include <stack>
#include "fx_internal.hpp"
#include "../matrix/sparsity_tools.hpp"

namespace CasADi{

/** \brief  Internal node class for the base class of SXFunctionInternal and MXFunctionInternal (lacks a public counterpart)
    The design of the class uses the curiously recurring template pattern (CRTP) idiom
    \author Joel Andersson 
    \date 2011
*/
template<typename DerivedType, typename MatType, typename NodeType>
class XFunctionInternal : public FXInternal{
  public:
    
    /** \brief  Constructor  */
    XFunctionInternal(const std::vector<MatType>& inputv, const std::vector<MatType>& outputv);
    
    /** \brief  Destructor */
    virtual ~XFunctionInternal(){}

    /** \brief  Topological sorting of the nodes based on Depth-First Search (DFS) */
    static void sort_depth_first(std::stack<NodeType*>& s, std::vector<NodeType*>& nodes);

    /** \brief  Topological (re)sorting of the nodes based on Breadth-First Search (BFS) (Kahn 1962) */
    static void resort_breadth_first(std::vector<NodeType*>& algnodes);

    /** \brief  Topological (re)sorting of the nodes with the purpose of postponing every calculation as much as possible, as long as it does not influence a dependent node */
    static void resort_postpone(std::vector<NodeType*>& algnodes, std::vector<int>& lind);
             
    /** \brief  Construct a complete Jacobian by compression */
    std::vector<MatType> jacGen(const std::vector<std::pair<int,int> >& jblocks, bool compact, 
                          const std::vector<bool>& symmetric_block);

    // Data members (all public)
    
    /** \brief  Inputs of the function (needed for symbolic calculations) */
    std::vector<MatType> inputv_;

    /** \brief  Outputs of the function (needed for symbolic calculations) */
    std::vector<MatType> outputv_;
};

// Template implementations

template<typename DerivedType, typename MatType, typename NodeType>
XFunctionInternal<DerivedType,MatType,NodeType>::XFunctionInternal(
    const std::vector<MatType>& inputv, const std::vector<MatType>& outputv) : inputv_(inputv),  outputv_(outputv){
      addOption("topological_sorting",OT_STRING,"breadth-first","Topological sorting algorithm","depth-first|breadth-first");
  
  // Make sure that inputs are symbolic
  for(int i=0; i<inputv.size(); ++i){
    if (inputv[i].isNull() || inputv[i].empty()) {
      sym(inputv_[i],"empty",0,0);
    } else if(!isSymbolicSparse(inputv[i])){
      casadi_error("XFunctionInternal::XFunctionInternal: Xfunction input arguments must be purely symbolic." << std::endl << "Argument #" << i << " is not symbolic.");
    }
  }
      
  // Allocate space for inputs
  setNumInputs(inputv_.size());
  for(int i=0; i<input_.size(); ++i)
    input(i) = DMatrix(inputv_[i].sparsity());
  
  // Null output arguments become empty
  for(int i=0; i<outputv_.size(); ++i) {
    if (outputv_[i].isNull()) {
      outputv_[i] = MatType(0,0);
    }
  }

  // Allocate space for outputs
  setNumOutputs(outputv_.size());
  for(int i=0; i<output_.size(); ++i)
    output(i) = DMatrix(outputv_[i].sparsity());
}


template<typename DerivedType, typename MatType, typename NodeType>
void XFunctionInternal<DerivedType,MatType,NodeType>::sort_depth_first(std::stack<NodeType*>& s, std::vector<NodeType*>& nodes){

    while(!s.empty()){
      
      // Get the topmost element
      NodeType* t = s.top();
      
      // If the last element on the stack has not yet been added
      if (t && !t->temp){
        
        // If a dependent node was added to the stack
        bool added_dep = false;
        
        // Initialize the node
        t->init();
    
        // Add dependent nodes if not already added
        for(int i=0; i<t->ndep(); ++i){
          if(t->dep(i).get() !=0 && static_cast<NodeType*>(t->dep(i).get())->temp == 0) {
            // if the first child has not yet been added
            s.push(static_cast<NodeType*>(t->dep(i).get()));
            added_dep=true;
            break;
          }
        }
        
        if(!added_dep){  // if no dependencies were added
          // Add to algorithm
          nodes.push_back(t);

          // Mark the node as found
          t->temp = 1;

          // Remove from stack
          s.pop();
          
        }
      } else {
        // If the last element on the stack has already been added
        s.pop();
      }
    }

  // Reset the node counters
  for(typename std::vector<NodeType*>::iterator it=nodes.begin(); it!=nodes.end(); ++it){
    (**it).temp = 0;
  }
}

template<typename DerivedType, typename MatType, typename NodeType>
void XFunctionInternal<DerivedType,MatType,NodeType>::resort_postpone(std::vector<NodeType*>& algnodes, std::vector<int>& lind){

  // Number of levels
  int nlevels = lind.size()-1;

  // Set the counter to be the corresponding place in the algorithm
  for(int i=0; i<algnodes.size(); ++i)
    algnodes[i]->temp = i;

  // Save the level of each element
  std::vector<int> level(algnodes.size());
  for(int i=0; i<nlevels; ++i)
    for(int j=lind[i]; j<lind[i+1]; ++j)
      level[j] = i;

  // Count the number of times each node is referenced inside the algorithm
  std::vector<int> numref(algnodes.size(),0);
  for(int i=0; i<algnodes.size(); ++i){
    for(int c=0; c<algnodes[i]->ndep(); ++c){ // for both children
      NodeType* child = static_cast<NodeType*>(algnodes[i]->dep(c).get());
      if(child && child->hasDep())
        numref[child->temp]++;
    }
  }

  // Stacks of additional nodes at the current and previous level
  std::stack<int> extra[2];

  // Loop over the levels in reverse order
  for(int i=nlevels-1; i>=0; --i){

    // The stack for the current level (we are removing elements from this stack)
    std::stack<int>& extra_this = extra[i%2]; // i odd -> use extra[1]

    // The stack for the previous level (we are adding elements to this stack)
    std::stack<int>& extra_prev = extra[1-i%2]; // i odd -> use extra[0]

    // Loop over the nodes of the level
    for(int j=lind[i]; j<lind[i+1]; ++j){
      // element to be treated
      int el = j;

      // elements in the stack have priority
      if(!extra_this.empty()){
        // Replace the element with one from the stack
        el = extra_this.top();
        extra_this.pop();
        --j; // redo the loop
      }

      // Skip the element if belongs to a higher level (i.e. was already treated)
      if(level[el] > i) continue;

      // for both children
      for(int c=0; c<algnodes[el]->ndep(); ++c){

        NodeType* child = static_cast<NodeType*>(algnodes[el]->dep(c).get());

        if(child && child->hasDep()){
          // Decrease the reference count of the children
          numref[child->temp]--;

          // If this was the last time the child was referenced ...
          // ... and it is not the previous level...
          if(numref[child->temp]==0 && level[child->temp] != i-1){

            // ... then assign a new level ...
            level[child->temp] = i-1;

            // ... and add to stack
            extra_prev.push(child->temp);

          } // if no more references
        } // if binary
      } // for c = ...
    } // for j
  } // for i

  // Count the number of elements on each level
  for(std::vector<int>::iterator it=lind.begin(); it!=lind.end(); ++it)
    *it = 0;
  for(std::vector<int>::const_iterator it=level.begin(); it!=level.end(); ++it)
    lind[*it + 1]++;

  // Cumsum to get the index corresponding to the first element of each level
  for(int i=0; i<nlevels; ++i)
    lind[i+1] += lind[i];

  // New index for each element
  std::vector<int> runind = lind; // running index for each level
  std::vector<int> newind(algnodes.size());
  for(int i=0; i<algnodes.size(); ++i)
    newind[i] = runind[level[algnodes[i]->temp]]++;

  // Resort the algorithm and reset the temporary
  std::vector<NodeType*> oldalgnodes = algnodes;
  for(int i=0; i<algnodes.size(); ++i){
    algnodes[newind[i]] = oldalgnodes[i];
    oldalgnodes[i]->temp = 0;
  }

}

template<typename DerivedType, typename MatType, typename NodeType>
void XFunctionInternal<DerivedType,MatType,NodeType>::resort_breadth_first(std::vector<NodeType*>& algnodes){

  // We shall assign a "level" to each element of the algorithm. A node which does not depend on other binary nodes are assigned level 0 and for nodes that depend on other nodes of the algorithm, the level will be the maximum level of any of the children plus 1. Note that all nodes of a level can be evaluated in parallel. The level will be saved in the temporary variable

  // Total number of levels
  int nlevels = 0;  

  // Get the earliest posible level
  for(typename std::vector<NodeType*>::iterator it=algnodes.begin(); it!=algnodes.end(); ++it){
    // maximum level of any of the children
    int maxlevel = -1;
    for(int c=0; c<(*it)->ndep(); ++c){    // Loop over the children
      NodeType* child = static_cast<NodeType*>((*it)->dep(c).get());
      if(child->hasDep() && child->temp > maxlevel)
        maxlevel = child->temp;
    }

    // Save the level of this element
    (*it)->temp = 1 + maxlevel;

    // Save if new maximum reached
    if(1 + maxlevel > nlevels)
      nlevels = 1 + maxlevel;
  }
  nlevels++;

  // Index of the first node on each level
  std::vector<int> lind;

  // Count the number of elements on each level
  lind.resize(nlevels+1,0); // all zeros to start with
  for(int i=0; i<algnodes.size(); ++i)
    lind[algnodes[i]->temp+1]++;

  // Cumsum to get the index of the first node on each level
  for(int i=0; i<nlevels; ++i)
    lind[i+1] += lind[i];

  // Get a new index for each element of the algorithm
  std::vector<int> runind = lind; // running index for each level
  std::vector<int> newind(algnodes.size());
  for(int i=0; i<algnodes.size(); ++i)
    newind[i] = runind[algnodes[i]->temp]++;

  // Resort the algorithm accordingly and reset the temporary
  std::vector<NodeType*> oldalgnodes = algnodes;
  for(int i=0; i<algnodes.size(); ++i){
    algnodes[newind[i]] = oldalgnodes[i];
    oldalgnodes[i]->temp = 0;
  }

#if 0

 int maxl=-1;
  for(int i=0; i<lind.size()-1; ++i){
    int l = (lind[i+1] - lind[i]);
//if(l>10)    std::cout << "#level " << i << ": " << l << std::endl;
  std::cout << l << ",";
    if(l>maxl) maxl= l;
  }
    std::cout << std::endl << "maxl = " << maxl << std::endl;

  for(int i=0; i<algnodes.size(); ++i){
    algnodes[i]->temp = i;
  }


  maxl=-1;
  for(int i=0; i<lind.size()-1; ++i){
    int l = (lind[i+1] - lind[i]);
    std::cout << std::endl << "#level " << i << ": " << l << std::endl;

int ii = 0;

    for(int j=lind[i]; j<lind[i+1]; ++j){

  std::vector<NodeType*>::const_iterator it = algnodes.begin() + j;

std::cout << "  "<< ii++ << ": ";

    int op = (*it)->op;
    stringstream s,s0,s1;
    s << "i_" << (*it)->temp;

    int i0 = (*it)->child[0].get()->temp;
    int i1 = (*it)->child[1].get()->temp;

    if((*it)->child[0]->hasDep())  s0 << "i_" << i0;
    else                             s0 << (*it)->child[0];
    if((*it)->child[1]->hasDep())  s1 << "i_" << i1;
    else                             s1 << (*it)->child[1];

    std::cout << s.str() << " = ";
    print_c[op](std::cout,s0.str(),s1.str());
    std::cout << ";" << std::endl;




    }

  std::cout << l << ",";
    if(l>maxl) maxl= l;
  }
    std::cout << std::endl << "maxl (before) = " << maxl << std::endl;


  for(int i=0; i<algnodes.size(); ++i){
    algnodes[i]->temp = 0;
  }


#endif

  // Resort in order to postpone all calculations as much as possible, thus saving cache
 resort_postpone(algnodes,lind);


#if 0

  for(int i=0; i<algnodes.size(); ++i){
    algnodes[i]->temp = i;
  }



  maxl=-1;
  for(int i=0; i<lind.size()-1; ++i){
    int l = (lind[i+1] - lind[i]);
    std::cout << std::endl << "#level " << i << ": " << l << std::endl;

int ii = 0;

    for(int j=lind[i]; j<lind[i+1]; ++j){

  std::vector<NodeType*>::const_iterator it = algnodes.begin() + j;

std::cout << "  "<< ii++ << ": ";

    int op = (*it)->op;
    stringstream s,s0,s1;
    s << "i_" << (*it)->temp;

    int i0 = (*it)->child[0].get()->temp;
    int i1 = (*it)->child[1].get()->temp;

    if((*it)->child[0]->hasDep())  s0 << "i_" << i0;
    else                             s0 << (*it)->child[0];
    if((*it)->child[1]->hasDep())  s1 << "i_" << i1;
    else                             s1 << (*it)->child[1];

    std::cout << s.str() << " = ";
    print_c[op](std::cout,s0.str(),s1.str());
    std::cout << ";" << std::endl;




    }

  std::cout << l << ",";
    if(l>maxl) maxl= l;
  }
    std::cout << std::endl << "maxl = " << maxl << std::endl;


//  return;




  for(int i=0; i<algnodes.size(); ++i){
    algnodes[i]->temp = 0;
  }



/*assert(0);*/
#endif

}

template<typename DerivedType, typename MatType, typename NodeType>
std::vector<MatType> XFunctionInternal<DerivedType,MatType,NodeType>::jacGen(
    const std::vector<std::pair<int,int> >& jblocks, bool compact, const std::vector<bool>& symmetric_block){
  using namespace std;
  if(verbose()) std::cout << "XFunctionInternal::jacGen begin" << std::endl;
  
  // Create return object
  std::vector<MatType> ret(jblocks.size());
  
  // Which blocks are involved in the Jacobian
  std::vector<std::pair<int,int> > jblocks_no_f;
  std::vector<bool> symmetric_block_no_f;
  std::vector<int> jblock_ind;
  
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
      symmetric_block_no_f.push_back(!symmetric_block.empty() && symmetric_block[i]);
/*      symmetric_block_no_f.push_back(false);*/
      
      // Save sparsity
      ret[i] = MatType(jacSparsity(iind,oind,compact));
      if(verbose()){
        std::cout << "XFunctionInternal::jac Block " << i << " has " << ret[i].size() << " nonzeros out of " << ret[i].numel() << " elements" << std::endl;
      }
    }
  }

  if(verbose()) std::cout << "XFunctionInternal::jac allocated return value" << std::endl;
  
  // Quick return if no jacobians to be calculated
  if(jblocks_no_f.empty()){
    if(verbose()) std::cout << "XFunctionInternal::jac end 1" << std::endl;
    return ret;
  }
  casadi_assert_message(jblocks_no_f.size()==1,"Multiple Jacobian blocks are not supported (and won't be)");
  
  // Get a bidirectional partition
  CRSSparsity D1, D2;
  getPartition(jblocks_no_f.front().second,jblocks_no_f.front().first,D1,D2,true,symmetric_block_no_f.front());
  if(verbose()) std::cout << "XFunctionInternal::jac graph coloring completed" << std::endl;

  // Get the number of forward and adjoint sweeps
  int nfwd = D1.isNull() ? 0 : D1.size1();
  int nadj = D2.isNull() ? 0 : D2.size1();

  if(false)  { // commented out: too verbose
    std::cout << "XFunctionInternal::jac partitioning" << std::endl;
    std::cout << "   jblocks_no_f " << jblocks_no_f.front() << std::endl;
    std::cout << "             D1 " << D1 << std::endl;
    if (!D1.isNull()) {
      std::cout << "                  col        " << D1.col() << std::endl;
      std::cout << "                  rowind     " << D1.rowind() << std::endl;
      std::cout << "                  spy        " << DMatrix(D1,1) << std::endl;
    }
    std::cout << "             D2 " << D2 << std::endl;
  }
  
  // Forward seeds
  std::vector<std::vector<MatType> > fseed(nfwd);
  for(int dir=0; dir<nfwd; ++dir){
    // initialize to zero
    fseed[dir].resize(getNumInputs());
    for(int iind=0; iind<fseed[dir].size(); ++iind){
      fseed[dir][iind] = MatType(input(iind).sparsity(),0);
    }
    
    // Get the input index
    int iind = jblocks_no_f.front().second;
    
    // For all the directions
    for(int el = D1.rowind(dir); el<D1.rowind(dir+1); ++el){
      
      // Get the direction
      int c = D1.col(el);

      // Give a seed in the direction
      fseed[dir][iind].at(c) = 1;
    }
  }
  
  // Adjoint seeds
  std::vector<std::vector<MatType> > aseed(nadj);
  for(int dir=0; dir<nadj; ++dir){
    //initialize to zero
    aseed[dir].resize(getNumOutputs());
    for(int oind=0; oind<aseed[dir].size(); ++oind){
      aseed[dir][oind] = MatType(output(oind).sparsity(),0);
    }
    
    // Pass seeds
    for(int v=0; v<D2.size(); ++v){
      
      // Get the output index
      int oind = jblocks_no_f.front().first;
      
      // For all the directions
      for(int el = D2.rowind(dir); el<D2.rowind(dir+1); ++el){
        
        // Get the direction
        int c = D2.col(el);

        // Give a seed in the direction
        aseed[dir][oind].at(c) = 1; // NOTE: should be +=, right?
      }
    }
  }

  // Forward sensitivities
  std::vector<std::vector<MatType> > fsens(nfwd);
  for(int dir=0; dir<nfwd; ++dir){
    // initialize to zero
    fsens[dir].resize(getNumOutputs());
    for(int oind=0; oind<fsens[dir].size(); ++oind){
      fsens[dir][oind] = MatType(output(oind).sparsity(),0);
    }
  }

  // Adjoint sensitivities
  std::vector<std::vector<MatType> > asens(nadj);
  for(int dir=0; dir<nadj; ++dir){
    // initialize to zero
    asens[dir].resize(getNumInputs());
    for(int iind=0; iind<asens[dir].size(); ++iind){
      asens[dir][iind] = MatType(input(iind).sparsity(),0);
    }
  }
  
  // Evaluate symbolically
  eval(inputv_,outputv_,fseed,fsens,aseed,asens,true,false);

  // Get transposes and mappings for all jacobian sparsity patterns if we are using forward mode
  if(verbose())   std::cout << "XFunctionInternal::jac transposes and mapping" << std::endl;
  std::vector<std::vector<int> > mapping;
  std::vector<CRSSparsity> sp_trans;
  if(nfwd>0){
    mapping.resize(jblock_ind.size());
    sp_trans.resize(jblock_ind.size());
    for(int i=0; i<jblock_ind.size(); ++i){
      int oind = jblocks_no_f[i].first;
      int iind = jblocks_no_f[i].second;
      sp_trans[i] = jacSparsity(iind,oind,true).transpose(mapping[i]);
      if(false)  { // commented out, not interesting
        std::cout << "   mapping[" << i << "] " << mapping[i] << std::endl;
        std::cout << "   sp_trans[" << i << "] " << DMatrix(sp_trans[i],1) << std::endl;
        std::cout << "   sp_trans[" << i << "].col() " << sp_trans[i].col() << std::endl;
        std::cout << "   sp_trans[" << i << "].rowind() " << sp_trans[i].rowind() << std::endl;
      }
    }
  }

  // The nonzeros of the sensitivity matrix
  std::vector<int> nzmap, nzmap2;
  
  // A vector used to resolve collitions between directions
  std::vector<int> hits;
  
  // Carry out the forward sweeps
  bool symmetric2=false;
  for(int dir=0; dir<nfwd; ++dir){
    
    // Get the output index
    int oind = jblocks_no_f.front().first;
    int iind = jblocks_no_f.front().second;

    // Is the block symmetric
    bool symmetric = symmetric_block_no_f.front();
    symmetric2 = symmetric;
    
    // If symmetric, see how many times each output appears
    if(symmetric){
      // Initialize to zero
      hits.resize(output(oind).sparsity().size());
      fill(hits.begin(),hits.end(),0);
      
      // Get the sparsity of the Jacobian block
      const CRSSparsity& jsp = jacSparsity(iind,oind,true);
      const std::vector<int>& jsp_rowind = jsp.rowind();
      const std::vector<int>& jsp_col = jsp.col();

      // "Multiply" Jacobian sparsity by seed vector
      for(int el = D1.rowind(dir); el<D1.rowind(dir+1); ++el){
	
	// Get the input nonzero
	int c = D1.col(el);
	
	// Propagate dependencies
	for(int el_jsp=jsp_rowind[c]; el_jsp<jsp_rowind[c+1]; ++el_jsp){
	  hits[jsp_col[el_jsp]]++;
	}
      }
    }

    // Locate the nonzeros of the forward sensitivity matrix
    output(oind).sparsity().getElements(nzmap,false);
    fsens[dir][oind].sparsity().getNZInplace(nzmap);

    if(symmetric){
      input(iind).sparsity().getElements(nzmap2,false);
      fsens[dir][oind].sparsity().getNZInplace(nzmap2);
    }
    
    
    // For all the input nonzeros treated in the sweep
    for(int el = D1.rowind(dir); el<D1.rowind(dir+1); ++el){

      // Get the input nonzero
      int c = D1.col(el);
      int f2_out;
      if(symmetric){
	f2_out = nzmap2[c];
      }
      
      // Loop over the output nonzeros corresponding to this input nonzero
      for(int el_out = sp_trans.front().rowind(c); el_out<sp_trans.front().rowind(c+1); ++el_out){
	
	// Get the output nonzero
	int r_out = sp_trans.front().col(el_out);
	
	// Get the forward sensitivity nonzero
	int f_out = nzmap[r_out];
	if(f_out<0) continue; // Skip if structurally zero
	
	// The nonzero of the Jacobian now treated
	int elJ = mapping.front()[el_out];
	
	if(symmetric){
	  if(hits[r_out]==1){
	    ret[jblock_ind.front()].at(el_out) = fsens[dir][oind].at(f_out);
	    ret[jblock_ind.front()].at(elJ) = fsens[dir][oind].at(f_out);
	  }
	} else {
	  // Get the output seed
	  ret[jblock_ind.front()].at(elJ) = fsens[dir][oind].at(f_out);
	}
      }
    }
  }
      
  // Add elements to the Jacobian matrix
  for(int dir=0; dir<nadj; ++dir){
    
    // For all the Jacobian blocks
    for(int v=0; v<D2.size(); ++v){

      // Get the input and output index of the block
      int oind = jblocks_no_f.front().first;
      int iind = jblocks_no_f.front().second;
      
      // Locate the nonzeros of the adjoint sensitivity matrix
      input(iind).sparsity().getElements(nzmap,false);
      asens[dir][iind].sparsity().getNZInplace(nzmap);
      
      // Get the (compact) Jacobian sparsity pattern
      const CRSSparsity& sp = jacSparsity(iind,oind,true);

      // For all the output nonzeros treated in the sweep
      for(int el = D2.rowind(dir); el<D2.rowind(dir+1); ++el){

        // Get the output nonzero
        int r = D2.col(el);

        // Loop over the input nonzeros that influences this output nonzero
        for(int elJ = sp.rowind(r); elJ<sp.rowind(r+1); ++elJ){
          
          // Get the input nonzero
          int inz = sp.col(elJ);
          
          // Get the corresponding adjoint sensitivity nonzero
          int anz = nzmap[inz];
          if(anz<0) continue;
          
          // Get the input seed
          ret[jblock_ind.front()].at(elJ) = asens[dir][iind].at(anz);
        }
      }
    }
  }
  
  // Return
  if(verbose()) std::cout << "XFunctionInternal::jac end" << std::endl;
  return ret;
}

} // namespace CasADi

#endif // X_FUNCTION_INTERNAL_HPP
