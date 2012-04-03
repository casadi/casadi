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
#include "../mx/jacobian_reference.hpp"
#include "../mx/evaluation.hpp"
#include "../mx/mapping.hpp"
#include "../mx/mx_tools.hpp"
#include "../sx/sx_tools.hpp"

#include "../stl_vector_tools.hpp"
#include "../casadi_types.hpp"

#include <stack>
#include <typeinfo>

using namespace std;

namespace CasADi{

MXFunctionInternal::MXFunctionInternal(const std::vector<MX>& inputv, const std::vector<MX>& outputv) :
  XFunctionInternal<MXFunctionInternal,MX,MXNode>(inputv,outputv) {
  
  setOption("name", "unnamed_mx_function");
  setOption("topological_sorting","depth-first"); // breadth-first not working
  setOption("numeric_jacobian", true);

  liftfun_ = 0;
  liftfun_ud_ = 0;
}


MXFunctionInternal::~MXFunctionInternal(){
}


double* getptr(Matrix<double>& x){
	if(x.size()>0)
		return &x.front();
	else
		return 0;
}

void MXFunctionInternal::init(){  
  log("MXFunctionInternal::init begin");
  
  // Call the init function of the base class
  XFunctionInternal::init();

  // Stack for nodes to be added to the list of nodes
  stack<MXNode*> s;

  // Add the inputs to the stack
  for(vector<MX>::iterator it = inputv_.begin(); it!=inputv_.end(); ++it)
    if(it->get()) s.push(static_cast<MXNode*>(it->get()));

  // Add the outputs to the stack
  for(vector<MX>::iterator it = outputv_.begin(); it!=outputv_.end(); ++it)
    if(it->get()) s.push(static_cast<MXNode*>(it->get()));

  // All evaluation nodes in the order of evaluation
  vector<MXNode*> nodes;

  // Order the nodes in the order of dependencies using a depth-first topological sorting
  sort_depth_first(s,nodes);

  // TODO: Make sure that the output nodes are placed directly after the corresponding multiple output node

  // Get the sorting algorithm
  bool breadth_first_search;
  if(getOption("topological_sorting")=="breadth-first"){
    breadth_first_search = true;
  } else if(getOption("topological_sorting")=="depth-first"){
    breadth_first_search = false;
  } else {
    casadi_error("Unrecongnized topological_sorting: " << getOption("topological_sorting"));
  }

  // Resort the nodes in a more cache friendly order (Kahn 1962)
  casadi_assert_message(breadth_first_search==false,"Breadth-first search currently not working for MX");
  // resort_breadth_first(nodes); // FIXME
  
  // Set the temporary variables to be the corresponding place in the sorted graph
  for(int i=0; i<nodes.size(); ++i)
    nodes[i]->temp = i;

  // An output/input pair corresponding to a jacobian block
  typedef pair<int,int> oipair;
  
  // The set of jacobian blocks for a single evaluation node
  typedef set<oipair> jblocks;
  
  // The sets of jacobian blocks for each evaluation node
  typedef map<void*,jblocks> jbmap;
  
  // The set of evaluation nodes for each function
  typedef map<void*,jbmap> evmap;
  
  // Which functions need to be differentiated with respect to which arguments
  evmap m;

  // Place in the work vector for each node
  vector<int> place_in_work;
  place_in_work.reserve(nodes.size());
  
  // Place in the algorithm for each node
  vector<int> place_in_alg;
  place_in_alg.reserve(nodes.size());
  
  // Work vector for the virtual machine
  work.resize(0);
  work.reserve(nodes.size());
  
  // The sequence of instructions for the virtual machine
  alg.resize(0);
  alg.reserve(nodes.size());
  
  // Process nodes
  for(vector<MXNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it){
    // Get pointer to node
    MXNode* n = *it;
    
    // Add an element to the algorithm
    if(n->isOutputNode()){
      
      // Get the output index
      int oind = n->getFunctionOutput();

      // Get the index of the parent node
      int pind = place_in_alg[n->dep(0)->temp];
      
      // Save output to the parent node
      alg[pind].i_res.at(oind) = work.size();
      
      // Save place in work vector
      place_in_work.push_back(work.size());

      // Not in algorithm
      place_in_alg.push_back(-1);

    } else if(n->isJacobian()){
      
      // Get a reference to the function
      FX& fcn = n->dep(0)->dep(0)->getFunction();

      // Get the input and output indices
      int iind = n->getFunctionInput();
      int oind = n->getFunctionOutput();

      // Find the set of outputs in the map
      jblocks& jb = m[fcn.get()][n->dep(0)->dep(0).get()];
      
      // If empty, we need to add the undifferentiated outputs
      if(jb.empty()){
        int n_out = fcn.getNumOutputs();
        for(int k=0; k<n_out; ++k){
          jb.insert(oipair(k,-1));
        }
      }

      // Save jacobian to map
      pair<jblocks::iterator,bool> jb_loc = jb.insert(oipair(oind,iind));
      casadi_assert_message(jb_loc.second==true,"Multiple identical Jacobian references currently not supported, easy to fix with a test case available"); 

      // Get the place where the jacobian was inserted
      int place = fcn.getNumOutputs();
      for(jblocks::iterator it=jb.begin(); it!=jb.end(); ++it){
        if(it->second>=0){
          if(iind==it->second && oind==it->first) break;
          place++;
        }
      }
      
      // Get the index of the parent node
      int pind = place_in_alg[n->dep(0)->dep(0)->temp];
      
      // Save output to the parent node (note that we add to the end, which is not the ultimate position, therefore a sorting takes place later in the code
      alg[pind].i_res.insert(alg[pind].i_res.begin()+place,work.size());
      
      // Save place in work vector
      place_in_work.push_back(work.size());

      // Not in algorithm
      place_in_alg.push_back(-1);
      
/*    } else if(n->isSymbolic()){
      
      // Save place in work vector
      place_in_work.push_back(work.size());

      // Not in algorithm
      place_in_alg.push_back(-1);*/
      
    } else {
      
      // Add an element to the algorithm
      place_in_alg.push_back(alg.size());
      alg.resize(alg.size()+1);
      alg.back().mx.assignNode(n);

      // Save the indices of the arguments
      alg.back().i_arg.resize(n->ndep());
      for(int i=0; i<n->ndep(); ++i){
        alg.back().i_arg[i] = n->dep(i).isNull() ? -1 : place_in_work[n->dep(i)->temp];
      }
      
      if(n->isMultipleOutput()){

        // Not in work vector
        place_in_work.push_back(-1);
       
        // Allocate space for outputs indices
        alg.back().i_res.resize(n->getNumOutputs(),-1);
        
        // No element in the work vector needed, thus continue
        continue;
      } else {
        
        // Save place in work vector
        place_in_work.push_back(work.size());

        // Allocate space for the outputs index
        alg.back().i_res.resize(1,work.size());
      }
    }
    
    // Allocate memory for an element in the work vector
    work.resize(work.size()+1);
    work.back().data = Matrix<double>(n->sparsity(),0);
  }
  
  // Allocate memory for directional derivatives
  MXFunctionInternal::updateNumSens(false);

  // Indices corresponding to the inputs
  input_ind_.resize(inputv_.size());
  for(int i=0; i<input_ind_.size(); ++i){
    input_ind_[i] = place_in_work[inputv_[i]->temp];
  }

  // Indices corresponding to the outputs
  output_ind_.resize(outputv_.size());
  for(int i=0; i<output_ind_.size(); ++i){
    output_ind_[i] = place_in_work[outputv_[i]->temp];
  }

  // Reset the temporary variables
  for(int i=0; i<nodes.size(); ++i){
    nodes[i]->temp = 0;
  }
  
  // Eliminate Jacobian if the computational graph contains jacobians
  if(m.empty()) return;

  // Jacobian functions for each function and set of jacobians
  map<void*,map<jblocks,FX> > jac;
    
  // Loop over undifferentiated functions
  for(evmap::iterator it=m.begin(); it!=m.end(); ++it){
    
    // Jacobian functions for each set of jacobians
    map<jblocks,FX>& jb = jac[it->first];
    
    // Loop over function evaluations for the function
    for(jbmap::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2){
      
      // Get a reference to the jacobian function
      FX& j = jb[it2->second];
      
      // If j is null, the Jacobian must be generated
      if(j.isNull()){
        
        // Create input arguments to the jacobian function
        vector<oipair> jblocks;
        jblocks.insert(jblocks.end(),it2->second.begin(),it2->second.end());
        
        // Generate jacobian functions
        FXInternal* f = reinterpret_cast<FXInternal*>(it->first);
        j = f->jacobian_switch(jblocks);
        j.init();
      }
    }
  }

  // Loop over the nodes again, replacing nodes
  for(vector<AlgEl>::iterator it=alg.begin(); it!=alg.end(); it++){

    // Get the old references to the function and evaluation node (since they will change)
    if(it->mx->isEvaluation()){
      void *fcn_ref = it->mx->getFunction().get();
      void *eval_ref = it->mx.get();
    
      // Check if the function should be replaced
      evmap::const_iterator ii=m.find(fcn_ref);
      if(ii!=m.end()){
        
        // Check if the function evaluation is to be replaced
        jbmap::const_iterator jj=ii->second.find(eval_ref);
        if(jj!=ii->second.end()){
          
          // Make the pointer unique so that other expressions won't be affected
          it->mx.makeUnique(false);
          
          // Replace the function if its an evaluation node
          it->mx->getFunction() = jac[fcn_ref][jj->second];
          
          // Make a copy of the result indices
          vector<int> i_res = it->i_res;
          
          // Update the index of the undifferentiated outputs
          int oind_old=0, oind_new=0;
          for(jblocks::const_iterator kk=jj->second.begin(); kk!=jj->second.end(); ++kk){
            if(kk->second<0){
              i_res[oind_new] = it->i_res[oind_old++];
            }
            oind_new++;
          }

          // Update the index of the jacobian blocks outputs
          oind_new=0;
          for(jblocks::const_iterator kk=jj->second.begin(); kk!=jj->second.end(); ++kk){
            if(kk->second>=0){
              i_res[oind_new] = it->i_res[oind_old++];
            }
            oind_new++;
          }
          
          // Update i_res
          i_res.swap(it->i_res);
        }
      }
    }
  }
  
  // Get the full Jacobian already now
  if(jac_for_sens_){
    getFullJacobian();
  }

log("MXFunctionInternal::init end");
}

void MXFunctionInternal::updateNumSens(bool recursive){
  // Call the base class if needed
  if(recursive) XFunctionInternal::updateNumSens(recursive);
  
  // Allocate work for directional derivatives
  for(vector<FunctionIO>::iterator it=work.begin(); it!=work.end(); it++){
    it->dataF.resize(nfdir_,it->data);
    it->dataA.resize(nadir_,it->data);
  }
}

void MXFunctionInternal::setLiftingFunction(LiftingFunction liftfun, void* user_data){
  liftfun_ = liftfun;
  liftfun_ud_ = user_data;
}

void MXFunctionInternal::updatePointers(const AlgEl& el, int nfwd, int nadj){
  mx_input_.resize(el.i_arg.size());
  mx_output_.resize(el.i_res.size());

  mx_fwdSeed_.resize(nfwd);
  mx_fwdSens_.resize(nfwd);
  for(int d=0; d<nfwd; ++d){
    mx_fwdSeed_[d].resize(mx_input_.size());
    mx_fwdSens_[d].resize(mx_output_.size());
  }

  mx_adjSens_.resize(nadj);
  mx_adjSeed_.resize(nadj);
  for(int d=0; d<nadj; ++d){
    mx_adjSens_[d].resize(mx_input_.size());
    mx_adjSeed_[d].resize(mx_output_.size());
  }
  
  for(int i=0; i<mx_input_.size(); ++i){
    if(el.i_arg[i]>=0){
      mx_input_[i] = &work[el.i_arg[i]].data;
      for(int d=0; d<nfwd; ++d) mx_fwdSeed_[d][i] = &work[el.i_arg[i]].dataF[d];
      for(int d=0; d<nadj; ++d) mx_adjSens_[d][i] = &work[el.i_arg[i]].dataA[d];
    } else {
      mx_input_[i] = 0;
      for(int d=0; d<nfwd; ++d) mx_fwdSeed_[d][i] = 0;
      for(int d=0; d<nadj; ++d) mx_adjSens_[d][i] = 0;
    }
  }
  
  for(int i=0; i<mx_output_.size(); ++i){
    if(el.i_res[i]>=0){
      mx_output_[i] = &work[el.i_res[i]].data;
      for(int d=0; d<nfwd; ++d) mx_fwdSens_[d][i] = &work[el.i_res[i]].dataF[d];
      for(int d=0; d<nadj; ++d) mx_adjSeed_[d][i] = &work[el.i_res[i]].dataA[d];
    } else {
      mx_output_[i] = 0;
      for(int d=0; d<nfwd; ++d) mx_fwdSens_[d][i] = 0;
      for(int d=0; d<nadj; ++d) mx_adjSeed_[d][i] = 0;
    }
  }
}

void MXFunctionInternal::evaluate(int nfdir, int nadir){
  log("MXFunctionInternal::evaluate begin");

  // Pass the inputs
  for(int ind=0; ind<input_.size(); ++ind){
    work[input_ind_[ind]].data.set(input(ind));
  }

  // Pass the forward seeds
  for(int dir=0; dir<nfdir; ++dir)
    for(int ind=0; ind<input_.size(); ++ind){
      work[input_ind_[ind]].dataF.at(dir).set(fwdSeed(ind,dir));
    }
  
  // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
  for(vector<AlgEl>::iterator it=alg.begin(); it!=alg.end(); it++){
    
    // Point pointers to the data corresponding to the element
    updatePointers(*it,nfdir,0);

    // Evaluate
    it->mx->evaluateD(mx_input_, mx_output_, mx_fwdSeed_, mx_fwdSens_, mx_adjSeed_, mx_adjSens_);
    // Lifting
    if(liftfun_ && it->mx->isNonLinear()){
      for(int i=0; i<it->i_res.size(); ++i){
        liftfun_(&mx_output_[i]->front(),it->mx.size(),liftfun_ud_);
      }
    }
  }
  
  log("MXFunctionInternal::evaluate evaluated forward");
  
  // Get the outputs
  for(int ind=0; ind<outputv_.size(); ++ind){
    work[output_ind_[ind]].data.get(output(ind));
  }

  // Get the forward sensitivities
  for(int dir=0; dir<nfdir; ++dir)
    for(int ind=0; ind<outputv_.size(); ++ind){
      work[output_ind_[ind]].dataF.at(dir).get(fwdSens(ind,dir));
    }
          
  if(nadir>0){
    log("MXFunctionInternal::evaluate evaluating adjoint");

    // Clear the adjoint seeds
    for(vector<FunctionIO>::iterator it=work.begin(); it!=work.end(); it++)
      for(int dir=0; dir<nadir; ++dir){
        fill(it->dataA.at(dir).begin(),it->dataA.at(dir).end(),0.0);
      }

    // Pass the adjoint seeds
    for(int ind=0; ind<outputv_.size(); ++ind)
      for(int dir=0; dir<nadir; ++dir){
        work[output_ind_[ind]].dataA.at(dir).set(adjSeed(ind,dir));
      }

    // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
    for(vector<AlgEl>::reverse_iterator it=alg.rbegin(); it!=alg.rend(); it++){
      
      // Point pointers to the data corresponding to the element
      updatePointers(*it,0,nadir);
      
      // Evaluate
      it->mx->evaluateD(mx_input_, mx_output_, mx_fwdSeed_, mx_fwdSens_, mx_adjSeed_, mx_adjSens_);
    }

    // Get the adjoint sensitivities
    for(int ind=0; ind<input_.size(); ++ind)
      for(int dir=0; dir<nadir; ++dir){
        work[input_ind_[ind]].dataA.at(dir).get(adjSens(ind,dir));
      }
    
    log("MXFunctionInternal::evaluate evaluated adjoint");
  }
  log("MXFunctionInternal::evaluate end");
}

void MXFunctionInternal::print(ostream &stream) const{
  FXInternal::print(stream);
  for(int i=0; i<alg.size(); ++i){
    stream << "i_" << i<< " =  ";
    alg[i].mx->printPart(stream,0);
    for(int j=0; j<alg[i].mx->ndep(); ++j){
      if(alg[i].i_arg[j]>=0){
        stream << "i_" << alg[i].i_arg[j];
      } else {
        stream << "[]";
      }
      alg[i].mx->printPart(stream,j+1);
    }
    stream << endl;
  }
}

MXFunctionInternal* MXFunctionInternal::clone() const{
  return new MXFunctionInternal(*this);
}

void MXFunctionInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
  XFunctionInternal::deepCopyMembers(already_copied);
  for(vector<AlgEl>::iterator it=alg.begin(); it!=alg.end(); ++it){
    if(it->mx->isEvaluation()){
      it->mx.makeUnique(already_copied,false);
      it->mx->getFunction() = deepcopy(it->mx->getFunction(),already_copied);
    }
  }
}

void MXFunctionInternal::spProp(bool fwd){
  if(fwd){
    for(vector<AlgEl>::iterator it=alg.begin(); it!=alg.end(); it++){
      if(it->mx->isSymbolic()) continue;
      
      // Point pointers to the data corresponding to the element
      updatePointers(*it,0,0);

      // Propagate sparsity forwards
      it->mx->propagateSparsity(mx_input_, mx_output_,true);
    }
  } else {
    for(vector<AlgEl>::reverse_iterator it=alg.rbegin(); it!=alg.rend(); it++){
      if(it->mx->isSymbolic()) continue;
      
      // Point pointers to the data corresponding to the element
      updatePointers(*it,0,0);
      
      // Propagate sparsity backwards
      it->mx->propagateSparsity(mx_input_, mx_output_,false);
      
      // Clear the seeds for the next sweep
      for(DMatrixPtrV::iterator it=mx_output_.begin(); it!=mx_output_.end(); ++it){
        DMatrix* seed = *it;
        if(seed){
          bvec_t *iseed = get_bvec_t(seed->data());
          fill_n(iseed,seed->size(),0);
        }
      }
    }
  }
}

void MXFunctionInternal::spReset(int iind, int oind){
  
  // Start by setting all elements of the work vector to zero
  for(vector<FunctionIO>::iterator it=work.begin(); it!=work.end(); ++it){
    //Get a pointer to the int array
    bvec_t *iwork = get_bvec_t(it->data.data());
    fill_n(iwork,it->data.size(),0);
  }

  // Pointer to the data vector for the input
  int el_in = input_ind_[iind];
  iwork_in_ = get_bvec_t(work[el_in].data.data());
  
  // Pointer to the data vector for the output
  int el_out = output_ind_[oind];
  iwork_out_ = get_bvec_t(work[el_out].data.data());
}

FX MXFunctionInternal::jacobian(const std::vector<std::pair<int,int> >& jblocks){
  // Make sure initialized
  assertInit();
  
  // Outputs of the Jacobian function
  vector<MX> j_out;
  j_out.reserve(jblocks.size());

  // The jacobians for each input
  vector<vector<MX> > jacs(getNumInputs());
  
  for(vector<pair<int,int> >::const_iterator it=jblocks.begin(); it!=jblocks.end(); ++it){
    // Get the block
    int iind = it->second;
    int oind = it->first;
    
    // If variable index is -1, we want nondifferentiated function output
    if(iind==-1){
      // Nondifferentiated function
      j_out.push_back(outputv_.at(oind));
      
    } else {
      
      // Create a Jacobian symbolically, if not already created
      if(jacs.at(iind).empty()){
        jacs.at(iind) = jac(iind);
      }
      
      // Add to outputs
      j_out.push_back(jacs.at(iind).at(oind));
    }
  }
  
  // Create function
  return MXFunction(inputv_,j_out);
}

vector<MX> MXFunctionInternal::jac(const vector<pair<int,int> >& jblocks, bool compact, const std::vector<bool>& symmetric_block){
  return jacGen(jblocks,compact,symmetric_block);
}

std::vector<MX> MXFunctionInternal::jac(int ider){
  assertInit();
  
  for (int i=0;i<inputv_.size();i++) {
    if (!inputv_[i].dense()) {
      casadi_error("MXFunctionInternal::jac is currently unsupported for sparse inputs.");
    }
  }
  
  casadi_assert_message(ider<input_.size(),"Index out of bounds");

  // Variable with respect to which we differentiate
  const MX& sder = inputv_[ider];

  // Number of forward directions
  int nfwd = sder.size();
  
  // Get the seed matrices
  vector<vector<MX> > fseed(nfwd);
  
  // Give zero seeds in all directions
  for(int d=0; d<nfwd; ++d){
    fseed[d].resize(input_.size());
    for(int iind=0; iind<input_.size(); ++iind){
      
      // Access the input expression
      const MX& s = inputv_[iind];
      if(d==0){
        fseed[d][iind] = MX::sparse(s.size1(),s.size2());
      } else {
        fseed[d][iind] = fseed[0][iind];
      }
    }
  }
  
  // Loop over the rows of the input
  for(int i=0; i<sder.size1(); ++i){
    // loop over the nonzeros of the input (i.e. derivative directions)
    for(int d=sder.sparsity().rowind(i); d<sder.sparsity().rowind(i+1); ++d){
      // get the column in the input matrix
      int j=sder.sparsity().col(d); 
  
      // Give a seed in the (i,j) direction
      fseed[d][ider](i,j) = 1;
    }
  }
  
  // Forward mode automatic differentiation, symbolically
  vector<vector<MX> > fsens = adFwd(fseed);

  // Collect the directions
  vector<MX> ret(outputv_.size());
  vector<MX> tmp(nfwd);
  for(int oind=0; oind<outputv_.size(); ++oind){
    for(int d=0; d<nfwd; ++d){
      tmp[d] = flatten(fsens[d][oind]);
    }
    ret[oind] = horzcat(tmp);
  }
  
  return ret;
}

std::vector<MX> MXFunctionInternal::grad(int igrad){
  assertInit();
  casadi_assert_message(igrad<output_.size(),"Index out of bounds");

  // Variable with respect to which we differentiate
  const MX& sgrad = outputv_[igrad];

  // Number of adjoint directions
  int nadj = sgrad.size();
  
  // Get the seed matrices
  vector<vector<MX> > aseed(nadj);
  
  // Give zero seeds in all directions
  for(int d=0; d<nadj; ++d){
    aseed[d].resize(output_.size());
    for(int oind=0; oind<output_.size(); ++oind){
      
      // Access the input expression
      const MX& s = outputv_[oind];
      if(d==0){
        aseed[d][oind] = MX::sparse(s.size1(),s.size2());
      } else {
        aseed[d][oind] = aseed[0][oind];
      }
    }
  }
  
  // Loop over the rows of the output
  for(int i=0; i<sgrad.size1(); ++i){
    // loop over the nonzeros of the output (i.e. derivative directions)
    for(int d=sgrad.sparsity().rowind(i); d<sgrad.sparsity().rowind(i+1); ++d){
      // get the column in the input matrix
      int j=sgrad.sparsity().col(d); 
  
      // Give a seed in the (i,j) direction
      aseed[d][igrad](i,j) = 1;
    }
  }
  
  // Forward mode automatic differentiation, symbolically
  vector<vector<MX> > asens = adAdj(aseed);

  // Collect the sensitivities
  vector<MX> ret(inputv_.size());
  vector<MX> tmp(nadj);
  for(int iind=0; iind<inputv_.size(); ++iind){
    for(int d=0; d<nadj; ++d){
      tmp[d] = trans(flatten(asens[d][iind]));
    }
    ret[iind] = vertcat(tmp);
  }
  
  return ret;
}

void MXFunctionInternal::evalMX(const std::vector<MX>& arg, std::vector<MX>& res, 
                                const std::vector<std::vector<MX> >& fseed, std::vector<std::vector<MX> >& fsens, 
                                const std::vector<std::vector<MX> >& aseed, std::vector<std::vector<MX> >& asens,
                                bool output_given, bool eliminate_constants){
  assertInit();
  
  // Symbolic work, non-differentiated
  std::vector<MX> swork(work.size());

  // Get the number of directions
  const int nfwd = fseed.size();
  const int nadj = aseed.size();

  // Arguments for the function evaluation
  MXPtrV input_p, output_p;
  MXPtrVV fseed_p(nfwd), fsens_p(nfwd);
  MXPtrVV aseed_p(nadj), asens_p(nadj);
  MXPtrVV dummy_p;

  // Work vector, forward derivatives
  std::vector<std::vector<MX> > dwork(work.size());
  fill(dwork.begin(),dwork.end(),std::vector<MX>(nfwd));

  // Pass the inputs
  if(!output_given){
    for(int iind=0; iind<arg.size(); ++iind){
      swork[input_ind_[iind]] = arg[iind];
    }
  }
  
  // Pass the forward seeds
  for(int iind=0; iind<input_.size(); ++iind){
    for(int d=0; d<nfwd; ++d){
      dwork[input_ind_[iind]][d] = fseed[d][iind];
    }
  }
    
  // Loop over computational nodes in forward order
  for(vector<AlgEl>::iterator it=alg.begin(); it!=alg.end(); ++it){
    
    // Pointers to the arguments of the evaluation
    input_p.resize(it->i_arg.size());
    for(int i=0; i<input_p.size(); ++i){
      int el = it->i_arg[i]; // index of the argument
      input_p[i] = el<0 ? 0 : &swork[el];
    }
        
    // Pointers to the result of the evaluation
    output_p.resize(it->i_res.size());
    for(int i=0; i<output_p.size(); ++i){
      int el = it->i_res[i]; // index of the output
      output_p[i] = el<0 ? 0 : &swork[el];
    }

    // Copy answer of the evaluation, if known
    if(output_given){
      if(it->mx->isMultipleOutput()){
        for(int oind=0; oind<output_p.size(); ++oind){
          if(output_p[oind]!=0){
            *output_p[oind] = MX::create(new OutputNode(it->mx,oind));
          }
        }
      } else {
        if(output_p[0]!=0){
          *output_p[0] = it->mx;
        }
      }
      
      // Quick continue if no evaluation needed
      if(nfwd==0) continue;
    }
    
    // Skip node if it is given
    if(it->mx->isSymbolic()) continue; 

    // Forward seeds and sensitivities
    for(int d=0; d<nfwd; ++d){
      fseed_p[d].resize(it->i_arg.size());
      for(int iind=0; iind<it->i_arg.size(); ++iind){
        int el = it->i_arg[iind];
        fseed_p[d][iind] = el<0 ? 0 : &dwork[el][d];
        
        // Give zero seed if null
        if(el>=0 && dwork[el][d].isNull()){
          if(d==0){
            dwork[el][d] = MX::sparse(input_p[iind]->size1(),input_p[iind]->size2());
          } else {
            dwork[el][d] = dwork[el][0];
          }
        }
      }

      fsens_p[d].resize(it->i_res.size());
      for(int oind=0; oind<it->i_res.size(); ++oind){
        int el = it->i_res[oind];
        fsens_p[d][oind] = el<0 ? 0 : &dwork[el][d];
      }
    }

    // Call the evaluation function
    it->mx->evaluateMX(input_p,output_p,fseed_p,fsens_p,dummy_p,dummy_p,output_given);
  }
  
  // Collect the results
  if(!output_given){
    res.resize(outputv_.size());
    for(int oind=0; oind<res.size(); ++oind){
      int el = output_ind_[oind];
      res[oind] = swork[el];
    }
  }

  // Collect the symbolic forward sensitivities
  fsens.resize(nfwd);
  for(int d=0; d<nfwd; ++d){
    fsens[d].resize(outputv_.size());
    for(int oind=0; oind<outputv_.size(); ++oind){
      int el = output_ind_[oind];
      fsens[d][oind] = dwork[el][d];
    }
  }
  
  // Work vector, adjoint derivatives
  fill(dwork.begin(),dwork.end(),std::vector<MX>(nadj));
  
  // Pass the adjoint seeds
  for(int oind=0; oind<output_.size(); ++oind){
    int el = output_ind_[oind];
    for(int d=0; d<nadj; ++d){
      if(dwork[el][d].isNull())
        dwork[el][d] = aseed[d][oind];
      else
        dwork[el][d] += aseed[d][oind];
    }
  }
  
  // Loop over computational nodes in reverse order
  if(nadj>0){
    for(vector<AlgEl>::reverse_iterator it=alg.rbegin(); it!=alg.rend(); ++it){
      
      // Get the arguments of the evaluation
      input_p.resize(it->i_arg.size());
      for(int i=0; i<input_p.size(); ++i){
        int el = it->i_arg[i]; // index of the argument
        input_p[i] = el<0 ? 0 : &swork[el];
      }
          
      // Result of the evaluation
      output_p.resize(it->i_res.size());
      for(int i=0; i<output_p.size(); ++i){
        int el = it->i_res[i]; // index of the output
        output_p[i] = el<0 ? 0 : &swork[el];
      }

      // Sensitivity arguments
      for(int d=0; d<nadj; ++d){
        aseed_p[d].resize(it->i_res.size());
        for(int oind=0; oind<it->i_res.size(); ++oind){
          int el = it->i_res[oind];
          aseed_p[d][oind] = el<0 ? 0 : &dwork[el][d];
          
          // Provide a zero seed if no seed exists
          if(el>=0 && dwork[el][d].isNull()){
            dwork[el][d] = MX::sparse(swork[el].size1(),swork[el].size2());
          }
        }

        asens_p[d].resize(it->i_arg.size());
        for(int iind=0; iind<it->i_arg.size(); ++iind){
          int el = it->i_arg[iind];
          asens_p[d][iind] = el<0 ? 0 : &dwork[el][d];
          
          // Set sensitivities to zero if not yet used
          if(el>=0 && dwork[el][d].isNull()){
            dwork[el][d] = MX::sparse(swork[el].size1(),swork[el].size2());
          }
        }
      }

      // Call the evaluation function
      it->mx->evaluateMX(input_p,output_p,dummy_p,dummy_p,aseed_p,asens_p,true);
    }
  }
  
  // Collect the symbolic adjoint sensitivities
  asens.resize(nadj);
  for(int d=0; d<nadj; ++d){
    asens[d].resize(inputv_.size());
    for(int iind=0; iind<inputv_.size(); ++iind){
      int el = input_ind_[iind];
      asens[d][iind] = dwork[el][d];
    }
  }
}

std::vector<std::vector<MX> > MXFunctionInternal::adFwd(const std::vector<std::vector<MX> > & fseed){
  // NOTE: Function to be removed
  std::vector<std::vector<MX> > ret;
  std::vector<std::vector<MX> > dummy;
  evalMX(inputv_,outputv_,fseed,ret,dummy,dummy,true,false);
  return ret;
}

std::vector<std::vector<MX> > MXFunctionInternal::adAdj(const std::vector<std::vector<MX> > & aseed){
  // NOTE: Function to be removed
  std::vector<std::vector<MX> > ret;
  std::vector<std::vector<MX> > dummy;
  evalMX(inputv_,outputv_,dummy,dummy,aseed,ret,true,false);
  return ret;
}

FX MXFunctionInternal::hessian(int iind, int oind) {
  // Assert initialized
  casadi_assert_message(isInit(),"Function not initialized.");
  
  // Make sure that the function is scalar valued
  //casadi_assert_message(output(oind).numel()==1,"Can only create hessians for scalar valued functions");
  
  // Get the Jacobian of the function
  MX J = jac(iind).at(oind);
  // MX J = grad(oind).at(iind);
  
  // Construct the gradient function
  MXFunction gfcn(inputv_,trans(J));
  gfcn.init();
  
  // And return the jacobian of that gradient function
  return gfcn.jacobian(iind,0);
}

void MXFunctionInternal::evalSX(const std::vector<SXMatrix>& input_s, std::vector<SXMatrix>& output_s, 
                                const std::vector<std::vector<SXMatrix> >& fwdSeed, std::vector<std::vector<SXMatrix> >& fwdSens, 
                                const std::vector<std::vector<SXMatrix> >& adjSeed, std::vector<std::vector<SXMatrix> >& adjSens,
                                bool output_given, bool eliminate_constants){
  casadi_assert_message(fwdSens.empty(),"Not implemented");
  casadi_assert_message(adjSeed.empty(),"Not implemented");
  
  // Create a work array
  vector<SXMatrix> swork(work.size());
  for(vector<AlgEl>::iterator it=alg.begin(); it!=alg.end(); it++){
    for(int i=0; i<it->i_res.size(); ++i){
      if (it->i_res[i]>=0)
        swork[it->i_res[i]] = SXMatrix(it->mx->sparsity(i));
    }
  }
  
  // Pass the inputs
  for(int ind=0; ind<input_.size(); ++ind){
    swork[input_ind_[ind]].set(input_s[ind]);
  }

  // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
  vector<SXMatrix*> sxarg;
  vector<SXMatrix*> sxres;
  for(vector<AlgEl>::iterator it=alg.begin(); it!=alg.end(); it++){
    sxarg.resize(it->i_arg.size());
    for(int c=0; c<sxarg.size(); ++c){
      int ind = it->i_arg[c];
      sxarg[c] = ind<0 ? 0 : &swork[ind];
    }
    sxres.resize(it->i_res.size());
    for(int c=0; c<sxres.size(); ++c){
      int ind = it->i_res[c];
      sxres[c] = ind<0 ? 0 : &swork[ind];
    }
    it->mx->evaluateSX(sxarg,sxres);
  }
  
  // Get the outputs
  for(int ind=0; ind<outputv_.size(); ++ind)
    swork[output_ind_[ind]].get(output_s[ind]);
}

SXFunction MXFunctionInternal::expand(const std::vector<SXMatrix>& inputvsx ){
  assertInit();
  
  // Create inputs with the same name and sparsity as the matrix valued symbolic inputs
  vector<SXMatrix> arg(inputv_.size());
  if(inputvsx.empty()){ // No symbolic input provided
    for(int i=0; i<arg.size(); ++i){
      arg[i] = ssym(inputv_[i]->getName(),inputv_[i].sparsity());
    }
  } else { // Use provided symbolic input
    // Make sure number of inputs matches
    casadi_assert(inputvsx.size()==inputv_.size());
    
    // Make sure that sparsity matches
    for(int i=0; i<inputvsx.size(); ++i){
      casadi_assert(inputvsx[i].sparsity() == inputv_[i].sparsity());
    }

    // Copy to argument vector
    copy(inputvsx.begin(),inputvsx.end(),arg.begin());
  }

  // Create output vector with correct sparsity
  vector<SXMatrix> res(outputv_.size());
  for(int i=0; i<res.size(); ++i){
    res[i] = SXMatrix(outputv_[i].sparsity());
  }
  
  // No sensitivities
  vector<vector<SXMatrix> > dummy;
  
  // Evaluate symbolically
  eval(arg,res,dummy,dummy,dummy,dummy,false,false);
  
  // Create function
  SXFunction f(arg,res);
  return f;
}


} // namespace CasADi

