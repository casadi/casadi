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

#include "../mx/mx_node.hpp"
#include "../stl_vector_tools.hpp"

#include <stack>
#include <typeinfo>

using namespace std;

namespace CasADi{

MXFunctionInternal::MXFunctionInternal(const std::vector<MX>& inputv_, const std::vector<MX>& outputv_) : inputv(inputv_), outputv(outputv_){
  setOption("ad_order",1); // one by default

  // Allocate space for inputs
  input_.resize(inputv.size());
  for(int i=0; i<input_.size(); ++i)
    input(i).resize(inputv[i].size1(),inputv[i].size2());

  // Allocate space for outputs
  output_.resize(outputv.size());
  for(int i=0; i<output_.size(); ++i)
    output(i).resize(outputv[i].size1(),outputv[i].size2());

}

MXFunctionInternal::~MXFunctionInternal(){
}

void MXFunctionInternal::makeAlgorithm(MXNode* root, vector<MXNode*> &nodes, map<const MXNode*,int>  &nodemap){

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

  // Order all nodes of a matrix syntax tree in the order of calculation
  for(vector<MX>::iterator it = outputv.begin(); it!=outputv.end(); ++it)
    makeAlgorithm((MXNode*)it->get(), nodes, nodemap);

  // Create runtime elements for each node
  alg.resize(nodes.size());
  for(int i=0; i<alg.size(); ++i){
    // Save the node
    alg[i].assignNode(nodes[i]);

    // Set derivative order
    nodes[i]->maxord_ = 2;
    
    // Number of derivative directions
    nodes[i]->nfdir_ = nfdir_;
    nodes[i]->nadir_ = nadir_;

    nodes[i]->init();
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
  
  // Pass the inputs
  for(int ind=0; ind<input_.size(); ++ind)
    inputv[ind]->output().set(input(ind));

  // Pass the inputs seeds
  if(fsens_order>0)
    for(int dir=0; dir<nfdir_; ++dir)
      for(int ind=0; ind<input_.size(); ++ind)
        inputv[ind]->fwdSens(dir).set(fwdSeed(ind,dir));
  
  // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
  for(vector<MX>::iterator it=alg.begin(); it!=alg.end(); it++)
    (*it)->evaluate(fsens_order,0);
  
  log("MXFunctionInternal::evaluate evaluated forward");

  // Get the outputs
  for(int ind=0; ind<outputv.size(); ++ind)
    outputv[ind]->output().get(output(ind));

  // Get the output seeds
  if(fsens_order>0)
    for(int dir=0; dir<nfdir_; ++dir)
      for(int ind=0; ind<outputv.size(); ++ind)
        outputv[ind]->fwdSens(dir).get(fwdSens(ind,dir));
          
  if(asens_order>0){

    // Clear the adjoint seeds
    for(vector<MX>::iterator it=alg.begin(); it!=alg.end(); it++){
      vector<double>& asens = (*it)->adjSeed();
      fill(asens.begin(),asens.end(),0.0);
    }

    // Pass the adjoint seeds
    for(int ind=0; ind<outputv.size(); ++ind)
      outputv[ind]->adjSeed().set(adjSeed(ind));

    // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
    for(vector<MX>::reverse_iterator it=alg.rbegin(); it!=alg.rend(); it++){
      (*it)->evaluate(0,asens_order);
    }

    // Get the adjoint sensitivities
    for(int ind=0; ind<input_.size(); ++ind)
      inputv[ind]->adjSeed().get(adjSens(ind));
    
    log("MXFunctionInternal::evaluate evaluated adjoint");
  }
  log("MXFunctionInternal::evaluate end");
}

void MXFunctionInternal::print(ostream &stream) const{
  int algcount = 0;
  for(vector<MX>::const_iterator it=alg.begin(); it!=alg.end(); ++it){
    stream << "[" << algcount << "]: " << *it << endl;
    algcount++;
  }
}

void MXFunctionInternal::repr(std::ostream &stream) const{
  stream << "MXFunction \"" << getOption("name") << "\"";
}

MXFunctionInternal* MXFunctionInternal::clone() const{
  return new MXFunctionInternal(*this);
}


} // namespace CasADi

