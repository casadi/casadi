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

#include "mx_function_internal.hpp"

#include "../stl_vector_tools.hpp"

#include <stack>
#include <typeinfo>

using namespace std;

namespace CasADi{

MXFunctionInternal::MXFunctionInternal(const std::vector<MX>& inputv_, const std::vector<MX>& outputv_) : inputv(inputv_), outputv(outputv_){
  setOption("name", "unnamed_mx_function");

  // Allocate space for inputs
  setNumInputs(inputv.size());
  for(int i=0; i<input_.size(); ++i)
    input(i) = DMatrix(inputv[i].sparsity());

  // Allocate space for outputs
  setNumOutputs(outputv.size());
  for(int i=0; i<output_.size(); ++i)
    output(i) = DMatrix(outputv[i].sparsity());
}


MXFunctionInternal::~MXFunctionInternal(){
}

void MXFunctionInternal::makeAlgorithm(MXNode* root, vector<MXNode*> &nodes, map<const MXNode*,int>  &nodemap){
    // Quick return if root is null
    if(root==0)
      return;
  
    // Stack
    stack<MXNode*> s;
    s.push(root);

    while(!s.empty()){

      // If the last element on the stack has not yet been added
      if (nodemap.find(s.top()) == nodemap.end()){

        // Loop over the children of the topmost element
        bool all_dependencies_added = true;
        for(vector<MX>::iterator it = s.top()->dep_.begin(); it!=s.top()->dep_.end(); ++it){
          if(!it->isNull() && nodemap.find((MXNode*)it->get()) == nodemap.end()){ // a dependency has not been added
            // Push the dependency to the top of the stack
            s.push((MXNode*)it->get());
            all_dependencies_added = false;
            break;
          }
        }
  
        // Continue to the next loop if a dependency has been added
        if(!all_dependencies_added) continue;

        // Add to the vector of nodes
        nodes.push_back(s.top());

        // Save the index
        nodemap[s.top()] = nodes.size()-1;

        // Remove from stack
        s.pop();
    }
  }
}

void MXFunctionInternal::init(){
  log("MXFunctionInternal::init begin");
  
  // Call the init function of the base class
  FXInternal::init();

  // Clear the algorithm
  alg.clear();
  nodemap.clear();
  std::vector<MXNode*> nodes;
  outputv_ind.clear();
  inputv_ind.clear();

  // Begin by adding the inputs to the algorithm
  for(vector<MX>::iterator it = inputv.begin(); it!=inputv.end(); ++it){
    makeAlgorithm((MXNode*)it->get(), nodes, nodemap);
    inputv_ind.push_back(nodemap[(MXNode*)it->get()]);
  }

  // Order all nodes of a matrix syntax tree in the order of calculation
  for(vector<MX>::iterator it = outputv.begin(); it!=outputv.end(); ++it){
    makeAlgorithm((MXNode*)it->get(), nodes, nodemap);
    outputv_ind.push_back(nodemap[(MXNode*)it->get()]);
  }
  
  // Make sure that the output nodes are placed directly after the corresponding multiple output node
  
  
  // Create runtime elements for each node
  alg.resize(nodes.size());
  for(int ii=0; ii<alg.size(); ++ii){
    MX &m = alg[ii].mx;
    AlgEl* it = &alg[ii];
    
    // Save the node
    it->mx.assignNode(nodes[ii]);
    it->val.data = Matrix<double>(m->sparsity());
    it->val.dataF.resize(nfdir_);
    it->val.dataA.resize(nadir_);
    it->val.init();

    // Save the indices of the children nodes
    it->ch.resize(m->ndep());
    for(int i=0; i<m->ndep(); ++i){
      if(m->dep(i).isNull())
        it->ch[i] = -1;
      else
        it->ch[i] = nodemap[(MXNode*)m->dep(i).get()];
    }
    
    it->input.resize(m->ndep(),0);
    it->fwdSeed.resize(m->ndep());
    it->adjSens.resize(m->ndep());
    for(int i=0; i<m->ndep(); ++i){
      it->fwdSeed[i].resize(nfdir_,0);
      it->adjSens[i].resize(nadir_,0);
      if(it->ch[i]>=0){
        it->input[i] = &alg[it->ch[i]].val.data[0];
        for(int d=0; d<nfdir_; ++d)
          it->fwdSeed[i][d] = &alg[it->ch[i]].val.dataF[d][0];
        for(int d=0; d<nadir_; ++d)
          it->adjSens[i][d] = &alg[it->ch[i]].val.dataA[d][0];
      }
    }

    it->output = &it->val.data[0];
    it->fwdSens.resize(nfdir_);
    for(int d=0; d<nfdir_; ++d)
      it->fwdSens[d] = &it->val.dataF[d][0];
    it->adjSeed.resize(nadir_);
    for(int d=0; d<nadir_; ++d)
      it->adjSeed[d] = &it->val.dataA[d][0];
    
  }
  
  log("MXFunctionInternal::init end");
}

bool MXFunctionInternal::hasEl(const MX& mx) const{
  // Locate the node:
  map<const MXNode*,int>::const_iterator it = nodemap.find((MXNode*)mx.get());

  // Make sure that the node was indeed found
  return it != nodemap.end();
}


int MXFunctionInternal::findEl(const MX& mx) const{
  // Locate the node:
  map<const MXNode*,int>::const_iterator it = nodemap.find((MXNode*)mx.get());

  // Make sure that the node was indeed found
  if(it == nodemap.end()){
    stringstream ss;
    ss << "MXFunctionInternal::findEl: Node \"" << mx << "\" not in the algorithm.";
    throw CasadiException(ss.str());
  }
  // Return the index
  return it->second;
}

void MXFunctionInternal::evaluate(int fsens_order, int asens_order){
  log("MXFunctionInternal::evaluate begin");
  int nfdir = fsens_order ? nfdir_ : 0;
  int nadir = asens_order ? nadir_ : 0;
  
  // Pass the inputs
  for(int ind=0; ind<input_.size(); ++ind)
    alg[inputv_ind[ind]].val.data.set(input(ind));

  // Pass the forward seeds
  for(int dir=0; dir<nfdir; ++dir)
    for(int ind=0; ind<input_.size(); ++ind)
      alg[inputv_ind[ind]].val.dataF.at(dir).set(fwdSeed(ind,dir));
  
  // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
  for(vector<AlgEl>::iterator it=alg.begin(); it!=alg.end(); it++)
    it->mx->evaluate(it->input, it->output, it->fwdSeed, it->fwdSens, it->adjSeed, it->adjSens, nfdir, 0);
  
  log("MXFunctionInternal::evaluate evaluated forward");

  // Get the outputs
  for(int ind=0; ind<outputv.size(); ++ind)
    alg[outputv_ind[ind]].val.data.get(output(ind));

  // Get the forward sensitivities
  for(int dir=0; dir<nfdir; ++dir)
    for(int ind=0; ind<outputv.size(); ++ind)
      alg[outputv_ind[ind]].val.dataF.at(dir).get(fwdSens(ind,dir));
          
  if(nadir>0){

    // Clear the adjoint seeds
    for(vector<AlgEl>::iterator it=alg.begin(); it!=alg.end(); it++)
      for(int dir=0; dir<nadir; ++dir)
        fill(it->val.dataA.at(dir).begin(),it->val.dataA.at(dir).end(),0.0);

    // Pass the adjoint seeds
    for(int ind=0; ind<outputv.size(); ++ind)
      for(int dir=0; dir<nadir; ++dir)
        alg[outputv_ind[ind]].val.dataA.at(dir).set(adjSeed(ind,dir));

    // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
    for(vector<AlgEl>::reverse_iterator it=alg.rbegin(); it!=alg.rend(); it++){
      it->mx->evaluate(it->input, it->output, it->fwdSeed, it->fwdSens, it->adjSeed, it->adjSens, 0, nadir);
    }

    // Get the adjoint sensitivities
    for(int ind=0; ind<input_.size(); ++ind)
      for(int dir=0; dir<nadir; ++dir)
        alg[inputv_ind[ind]].val.dataA.at(dir).get(adjSens(ind,dir));
    
    log("MXFunctionInternal::evaluate evaluated adjoint");
  }
  log("MXFunctionInternal::evaluate end");
}

void MXFunctionInternal::print(ostream &stream) const{
  for(int i=0; i<alg.size(); ++i){
    stream << "i_" << i<< " =  ";
    vector<string> chname(alg[i].ch.size());
    for(int j=0; j<chname.size(); ++j){
      if(alg[i].ch[j]>=0){
        stringstream ss;
        ss << "i_" << alg[i].ch[j];
        chname[j] = ss.str();
      } else {
        chname[j] = "[]";
      }
    }
    alg[i].mx->print(stream, chname);
    stream << endl;
  }
}

MXFunctionInternal* MXFunctionInternal::clone() const{
  return new MXFunctionInternal(*this);
}


} // namespace CasADi

