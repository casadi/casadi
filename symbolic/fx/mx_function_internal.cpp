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
#include "../mx/evaluation_mx.hpp"
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
  setOption("numeric_jacobian", true);
  
  // Check if any inputs is a mapping
  bool has_mapping_inputs = false;
  for(vector<MX>::const_iterator it = inputv_.begin(); it!=inputv_.end(); ++it){
    has_mapping_inputs = has_mapping_inputs || it->isMapping();
  }
  
  // If one or more inputs is a mapping, elimination needed before creating function
  if(has_mapping_inputs){
    
    // Name of replaced variable
    stringstream rep_name;
    int ind = 0;

    // Temp vectors
    std::vector<int> inz, onz;
    
    // Find new variables and expressions for all inputs that needs to be eliminated
    std::vector<MX> v, vdef;
    for(vector<MX>::iterator it = inputv_.begin(); it!=inputv_.end(); ++it, ++ind){
      if(it->isMapping()){
        
        // Get the mapping node
        Mapping* n = dynamic_cast<Mapping*>(it->get());
        
        // Initialize the mapping, i.e. sort by input and output index
        n->init();
        
        // Create a new variable
        rep_name.clear();
        rep_name << "r_" << ind;
        MX new_var = msym(rep_name.str(),n->sparsity());
        
        // For all dependencies to be replaced
        for(int iind=0; iind<n->ndep(); ++iind){
          
          // Variable to be replaced
          MX v_dep = n->dep(iind);
          casadi_assert_message(v_dep.isSymbolic(),"Mapping inputs may only map to symbolic variables.");
          
          // Express this variable in terms of the new variable
          const vector<pair<int,int> >& assigns = n->index_output_sorted_[0][iind];
          inz.clear();
          onz.clear();
          for(vector<pair<int,int> >::const_iterator it_ass=assigns.begin(); it_ass!=assigns.end(); ++it_ass){
            inz.push_back(it_ass->second);
            onz.push_back(it_ass->first);
          }
          MX vdef_dep = MX::create(new Mapping(v_dep.sparsity()));
          vdef_dep->assign(new_var,inz,onz);
          
          // Save variable to list of variables to be replaced
          v.push_back(v_dep);
          vdef.push_back(vdef_dep);
        }
        
        // Replace variable
        *it = new_var;
      }
    }
    
    // Replace expressions
    outputv_ = substitute(outputv_,v,vdef);
  }
  
  // Check for duplicate entries among the input expressions
  bool has_duplicates = false;
  for(vector<MX>::iterator it = inputv_.begin(); it != inputv_.end(); ++it){
    has_duplicates = has_duplicates || it->getTemp()!=0;
    it->setTemp(1);
  }
  
  // Reset temporaries
  for(vector<MX>::iterator it = inputv_.begin(); it != inputv_.end(); ++it){
    it->setTemp(0);
  }
  casadi_assert_message(!has_duplicates, "The input expressions are not independent.");
  
  liftfun_ = 0;
  liftfun_ud_ = 0;
}


MXFunctionInternal::~MXFunctionInternal(){
}


void MXFunctionInternal::init(){  
  log("MXFunctionInternal::init begin");
      
  // Call the init function of the base class
  XFunctionInternal<MXFunctionInternal,MX,MXNode>::init();    

  // Stack for nodes to be added to the list of nodes
  stack<MXNode*> s;

  // Add the inputs to the stack
  for(vector<MX>::reverse_iterator it = inputv_.rbegin(); it!=inputv_.rend(); ++it)
    if(it->get()) s.push(static_cast<MXNode*>(it->get()));

  // Add the outputs to the stack
  for(vector<MX>::reverse_iterator it = outputv_.rbegin(); it!=outputv_.rend(); ++it)
    if(it->get()) s.push(static_cast<MXNode*>(it->get()));

  // All evaluation nodes in the order of evaluation
  vector<MXNode*> nodes;

  // Order the nodes in the order of dependencies using a depth-first topological sorting
  sort_depth_first(s,nodes);

  // Reset the node markers
  for(std::vector<MXNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it){
    (*it)->temp = 0;
  }
  
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

  // Place in the work vector for each node
  vector<int> place_in_work;
  place_in_work.reserve(nodes.size());
  
  // Place in the algorithm for each node
  vector<int> place_in_alg;
  place_in_alg.reserve(nodes.size());
  
  // Work vector for the virtual machine
  work_.resize(0);
  work_.reserve(nodes.size());
  
  // The sequence of instructions for the virtual machine
  algorithm_.resize(0);
  algorithm_.reserve(nodes.size());
  
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
      algorithm_[pind].res.at(oind) = work_.size();
      
      // Save place in work vector
      place_in_work.push_back(work_.size());

      // Not in algorithm
      place_in_alg.push_back(-1);

    } else {
      
      // Add an element to the algorithm
      place_in_alg.push_back(algorithm_.size());
      algorithm_.resize(algorithm_.size()+1);
      algorithm_.back().op.assignNode(n);

      // Save the indices of the arguments
      algorithm_.back().arg.resize(n->ndep());
      for(int i=0; i<n->ndep(); ++i){
        algorithm_.back().arg[i] = n->dep(i).isNull() ? -1 : place_in_work[n->dep(i)->temp];
      }
      
      if(n->isMultipleOutput()){

        // Not in work vector
        place_in_work.push_back(-1);
       
        // Allocate space for outputs indices
        algorithm_.back().res.resize(n->getNumOutputs(),-1);
        
        // No element in the work vector needed, thus continue
        continue;
      } else {
        
        // Save place in work vector
        place_in_work.push_back(work_.size());

        // Allocate space for the outputs index
        algorithm_.back().res.resize(1,work_.size());
      }
    }
    
    // Allocate memory for an element in the work vector
    work_.resize(work_.size()+1);
    work_.back().data = Matrix<double>(n->sparsity(),0);
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
  
  // Get the free variables
  collectFree();

log("MXFunctionInternal::init end");
}

void MXFunctionInternal::updateNumSens(bool recursive){
  // Call the base class if needed
  if(recursive) XFunctionInternal<MXFunctionInternal,MX,MXNode>::updateNumSens(recursive);
  
  // Allocate work for directional derivatives
  for(vector<FunctionIO>::iterator it=work_.begin(); it!=work_.end(); it++){
    it->dataF.resize(nfdir_,it->data);
    it->dataA.resize(nadir_,it->data);
  }
}

void MXFunctionInternal::setLiftingFunction(LiftingFunction liftfun, void* user_data){
  liftfun_ = liftfun;
  liftfun_ud_ = user_data;
}

void MXFunctionInternal::updatePointers(const AlgEl& el, int nfwd, int nadj){
  mx_input_.resize(el.arg.size());
  mx_output_.resize(el.res.size());

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
    if(el.arg[i]>=0){
      mx_input_[i] = &work_[el.arg[i]].data;
      for(int d=0; d<nfwd; ++d) mx_fwdSeed_[d][i] = &work_[el.arg[i]].dataF[d];
      for(int d=0; d<nadj; ++d) mx_adjSens_[d][i] = &work_[el.arg[i]].dataA[d];
    } else {
      mx_input_[i] = 0;
      for(int d=0; d<nfwd; ++d) mx_fwdSeed_[d][i] = 0;
      for(int d=0; d<nadj; ++d) mx_adjSens_[d][i] = 0;
    }
  }
  
  for(int i=0; i<mx_output_.size(); ++i){
    if(el.res[i]>=0){
      mx_output_[i] = &work_[el.res[i]].data;
      for(int d=0; d<nfwd; ++d) mx_fwdSens_[d][i] = &work_[el.res[i]].dataF[d];
      for(int d=0; d<nadj; ++d) mx_adjSeed_[d][i] = &work_[el.res[i]].dataA[d];
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
    work_[input_ind_[ind]].data.set(input(ind));
  }

  // Pass the forward seeds
  for(int dir=0; dir<nfdir; ++dir)
    for(int ind=0; ind<input_.size(); ++ind){
      work_[input_ind_[ind]].dataF.at(dir).set(fwdSeed(ind,dir));
    }
  
  // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
  for(vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); it++){
    
    // Point pointers to the data corresponding to the element
    updatePointers(*it,nfdir,0);

    // Evaluate
    it->op->evaluateD(mx_input_, mx_output_, mx_fwdSeed_, mx_fwdSens_, mx_adjSeed_, mx_adjSens_);
  
    // Lifting
    if(liftfun_ && it->op->isNonLinear()){
      for(int i=0; i<it->res.size(); ++i){
        liftfun_(&mx_output_[i]->front(),it->op.size(),liftfun_ud_);
      }
    }
  }
  
  log("MXFunctionInternal::evaluate evaluated forward");
  
  // Get the outputs
  for(int ind=0; ind<outputv_.size(); ++ind){
    work_[output_ind_[ind]].data.get(output(ind));
  }

  // Get the forward sensitivities
  for(int dir=0; dir<nfdir; ++dir)
    for(int ind=0; ind<outputv_.size(); ++ind){
      work_[output_ind_[ind]].dataF.at(dir).get(fwdSens(ind,dir));
    }
          
  if(nadir>0){
    log("MXFunctionInternal::evaluate evaluating adjoint");

    // Clear the adjoint seeds
    for(vector<FunctionIO>::iterator it=work_.begin(); it!=work_.end(); it++)
      for(int dir=0; dir<nadir; ++dir){
        fill(it->dataA.at(dir).begin(),it->dataA.at(dir).end(),0.0);
      }

    // Pass the adjoint seeds
    for(int ind=0; ind<outputv_.size(); ++ind)
      for(int dir=0; dir<nadir; ++dir){
        work_[output_ind_[ind]].dataA.at(dir).set(adjSeed(ind,dir));
      }

    // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
    for(vector<AlgEl>::reverse_iterator it=algorithm_.rbegin(); it!=algorithm_.rend(); it++){
      
      // Point pointers to the data corresponding to the element
      updatePointers(*it,0,nadir);
      
      // Evaluate
      it->op->evaluateD(mx_input_, mx_output_, mx_fwdSeed_, mx_fwdSens_, mx_adjSeed_, mx_adjSens_);
    }

    // Get the adjoint sensitivities
    for(int ind=0; ind<input_.size(); ++ind)
      for(int dir=0; dir<nadir; ++dir){
        work_[input_ind_[ind]].dataA.at(dir).get(adjSens(ind,dir));
      }
    
    log("MXFunctionInternal::evaluate evaluated adjoint");
  }
  log("MXFunctionInternal::evaluate end");
}

void MXFunctionInternal::print(ostream &stream) const{
  FXInternal::print(stream);
  for(vector<AlgEl>::const_iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it){
    stream << "{";
    for(int i=0; i<it->res.size(); ++i){
      if(i!=0) stream << ",";
      if(it->res[i]>=0){
        stream << "i_" << it->res[i];
      } else {
        stream << "NULL";
      }
    }
    stream << "} = ";
    it->op->printPart(stream,0);
    for(int i=0; i<it->arg.size(); ++i){
      if(it->arg[i]>=0){
        stream << "i_" << it->arg[i];
      } else {
        stream << "NULL";
      }
      it->op->printPart(stream,i+1);
    }
    stream << endl;
  }
}

MXFunctionInternal* MXFunctionInternal::clone() const{
  return new MXFunctionInternal(*this);
}

void MXFunctionInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
  XFunctionInternal<MXFunctionInternal,MX,MXNode>::deepCopyMembers(already_copied);
  for(vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it){
    if(it->op->isEvaluation()){
      it->op.makeUnique(already_copied,false);
      it->op->getFunction() = deepcopy(it->op->getFunction(),already_copied);
    }
  }
}

void MXFunctionInternal::spInit(bool fwd){
  // Start by setting all elements of the work vector to zero
  for(vector<FunctionIO>::iterator it=work_.begin(); it!=work_.end(); ++it){
    //Get a pointer to the int array
    bvec_t *iwork = get_bvec_t(it->data.data());
    fill_n(iwork,it->data.size(),0);
  }
}

void MXFunctionInternal::spEvaluate(bool fwd){
  if(fwd){ // Forward propagation
    
    // Pass input seeds
    for(int ind=0; ind<input_ind_.size(); ++ind){
      vector<double> &w = work_[input_ind_[ind]].data.data();
      bvec_t* iwork = get_bvec_t(w);
      bvec_t* swork = get_bvec_t(input(ind).data());
      for(int k=0; k<w.size(); ++k){
	iwork[k] = swork[k];
      }
    }
    
    // Propagate sparsity forward
    for(vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); it++){
      if(it->op->isSymbolic()) continue;
      
      // Point pointers to the data corresponding to the element
      updatePointers(*it,0,0);

      // Propagate sparsity forwards
      it->op->propagateSparsity(mx_input_, mx_output_,true);
    }
    
    // Get the output seeds
    for(int ind=0; ind<output_ind_.size(); ++ind){
      vector<double> &w = work_[output_ind_[ind]].data.data();
      bvec_t* iwork = get_bvec_t(w);
      bvec_t* swork = get_bvec_t(output(ind).data());
      for(int k=0; k<w.size(); ++k){
	swork[k] = iwork[k];
      }
    }
    
  } else { // Backward propagation

    // Pass output seeds
    for(int ind=0; ind<output_ind_.size(); ++ind){
      vector<double> &w = work_[output_ind_[ind]].data.data();
      bvec_t* iwork = get_bvec_t(w);
      bvec_t* swork = get_bvec_t(output(ind).data());
      for(int k=0; k<w.size(); ++k){
	iwork[k] = swork[k];
      }
    }

    // Propagate sparsity backwards
    for(vector<AlgEl>::reverse_iterator it=algorithm_.rbegin(); it!=algorithm_.rend(); it++){
      if(it->op->isSymbolic()) continue;
      
      // Point pointers to the data corresponding to the element
      updatePointers(*it,0,0);
      
      // Propagate sparsity backwards
      it->op->propagateSparsity(mx_input_, mx_output_,false);
      
      // Clear the seeds for the next sweep
      for(DMatrixPtrV::iterator it=mx_output_.begin(); it!=mx_output_.end(); ++it){
        DMatrix* seed = *it;
        if(seed){
          bvec_t *iseed = get_bvec_t(seed->data());
          fill_n(iseed,seed->size(),0);
        }
      }
    }
    
    // Get the input seeds and clear it from the work vector
    for(int ind=0; ind<input_ind_.size(); ++ind){
      vector<double> &w = work_[input_ind_[ind]].data.data();
      bvec_t* iwork = get_bvec_t(w);
      bvec_t* swork = get_bvec_t(input(ind).data());
      for(int k=0; k<w.size(); ++k){
	 swork[k] |= iwork[k];
	 iwork[k] = 0;
      }
    }
  }
}

FX MXFunctionInternal::getJacobian(int iind, int oind, bool compact, bool symmetric){
  // Create expressions for the Jacobian
  vector<MX> ret_out;
  ret_out.reserve(1+outputv_.size());
  ret_out.push_back(jac(iind,oind,compact,symmetric));
  ret_out.insert(ret_out.end(),outputv_.begin(),outputv_.end());
  
  // Return function
  return MXFunction(inputv_,ret_out);
}

FX MXFunctionInternal::getNumericJacobian(int iind, int oind, bool compact, bool symmetric){
  // Create expressions for the Jacobian
  vector<MX> ret_out;
  ret_out.reserve(1+outputv_.size());
  ret_out.push_back(jac(iind,oind,compact,symmetric,false,true));
  // ret_out.insert(ret_out.end(),outputv_.begin(),outputv_.end()); // NOTE: Commented out for similarity with the Jacobian class
  
  // Return function
  return MXFunction(inputv_,ret_out);  
}

MX MXFunctionInternal::jac(int iind, int oind, bool compact, bool symmetric, bool always_inline, bool never_inline){
  log("MXFunctionInternal::jac begin");
  return jacGen(iind,oind,compact,symmetric,always_inline,never_inline);
}

void MXFunctionInternal::evalMX(const std::vector<MX>& arg, std::vector<MX>& res, 
                                const std::vector<std::vector<MX> >& fseed, std::vector<std::vector<MX> >& fsens, 
                                const std::vector<std::vector<MX> >& aseed, std::vector<std::vector<MX> >& asens,
                                bool output_given){
  assertInit();
  casadi_assert_message(arg.size()==getNumInputs(),"Wrong number of input arguments");
  if(output_given){
    casadi_assert_message(res.size()==getNumOutputs(),"Wrong number of output arguments");
  }
  
  // Symbolic work, non-differentiated
  std::vector<MX> swork(work_.size());
  
  // Pass free variables
  for(int k=0; k<free_vars_.size(); ++k){
    swork[free_vars_ind_[k]] = free_vars_[k];
  }

  // Get the number of directions
  const int nfwd = fseed.size();
  const int nadj = aseed.size();

  MXPtrV input_p, output_p;
  MXPtrVV fseed_p(nfwd), fsens_p(nfwd);
  MXPtrVV aseed_p(nadj), asens_p(nadj);
  MXPtrVV dummy_p;

  // Work vector, forward derivatives
  std::vector<std::vector<MX> > dwork(work_.size());
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
  for(vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it){
    
    // Pointers to the arguments of the evaluation
    input_p.resize(it->arg.size());
    for(int i=0; i<input_p.size(); ++i){
      int el = it->arg[i]; // index of the argument
      input_p[i] = el<0 ? 0 : &swork[el];
    }
        
    // Pointers to the result of the evaluation
    output_p.resize(it->res.size());
    for(int i=0; i<output_p.size(); ++i){
      int el = it->res[i]; // index of the output
      output_p[i] = el<0 ? 0 : &swork[el];
    }

    // Copy answer of the evaluation, if known
    if(output_given){
      if(it->op->isMultipleOutput()){
        for(int oind=0; oind<output_p.size(); ++oind){
          if(output_p[oind]){
            *output_p[oind] = MX::create(new OutputNode(it->op,oind));
          }
        }
      } else {
        if(output_p[0]){
          *output_p[0] = it->op;
        }
      }
      
      // Quick continue if no evaluation needed
      if(nfwd==0) continue;
    }
    
    // Skip node if it is given
    if(it->op->isSymbolic()){
      continue;
    }

    // Forward seeds and sensitivities
    for(int d=0; d<nfwd; ++d){
      fseed_p[d].resize(it->arg.size());
      for(int iind=0; iind<it->arg.size(); ++iind){
        int el = it->arg[iind];
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

      fsens_p[d].resize(it->res.size());
      for(int oind=0; oind<it->res.size(); ++oind){
        int el = it->res[oind];
        fsens_p[d][oind] = el<0 ? 0 : &dwork[el][d];
      }
    }

    // Call the evaluation function
    it->op->evaluateMX(input_p,output_p,fseed_p,fsens_p,dummy_p,dummy_p,output_given);
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
    for(vector<AlgEl>::reverse_iterator it=algorithm_.rbegin(); it!=algorithm_.rend(); ++it){
      
      // Get the arguments of the evaluation
      input_p.resize(it->arg.size());
      for(int i=0; i<input_p.size(); ++i){
        int el = it->arg[i]; // index of the argument
        input_p[i] = el<0 ? 0 : &swork[el];
      }
          
      // Result of the evaluation
      output_p.resize(it->res.size());
      for(int i=0; i<output_p.size(); ++i){
        int el = it->res[i]; // index of the output
        output_p[i] = el<0 ? 0 : &swork[el];
      }

      // Sensitivity arguments
      for(int d=0; d<nadj; ++d){
        aseed_p[d].resize(it->res.size());
        for(int oind=0; oind<it->res.size(); ++oind){
          int el = it->res[oind];
          aseed_p[d][oind] = el<0 ? 0 : &dwork[el][d];
          
          // Provide a zero seed if no seed exists
          if(el>=0 && dwork[el][d].isNull()){
            dwork[el][d] = MX::sparse(swork[el].size1(),swork[el].size2());
          }
        }

        asens_p[d].resize(it->arg.size());
        for(int iind=0; iind<it->arg.size(); ++iind){
          int el = it->arg[iind];
          asens_p[d][iind] = el<0 ? 0 : &dwork[el][d];
          
          // Set sensitivities to zero if not yet used
          if(el>=0 && dwork[el][d].isNull()){
            dwork[el][d] = MX::sparse(swork[el].size1(),swork[el].size2());
          }
        }
      }

      // Call the evaluation function
      it->op->evaluateMX(input_p,output_p,dummy_p,dummy_p,aseed_p,asens_p,true);
      
// For debugging #491
#if 0
        cout << "{";
        for(int i=0; i<it->res.size(); ++i){
          if(i!=0) cout << ",";
          if(it->res[i]>=0){
            cout << "i_" << it->res[i];
          } else {
            cout << "NULL";
          }
        }
        cout << "} = ";
        it->op->printPart(cout,0);
        for(int i=0; i<it->arg.size(); ++i){
          if(it->arg[i]>=0){
            cout << "i_" << it->arg[i];
          } else {
            cout << "NULL";
          }
          it->op->printPart(cout,i+1);
        }
        cout << endl;
      
        cout << "aseed" << endl;
      for(int d=0; d<aseed_p.size(); ++d){
        for(int i=0; i<aseed_p[d].size(); ++i){
          if(aseed_p[d][i]==0){
            cout << "NULL";
          } else {
            cout << *aseed_p[d][i];
          }
          cout << "," << endl;
        }
      }
        cout << endl;
        
        cout << "asens" << endl;
      for(int d=0; d<asens_p.size(); ++d){
        for(int i=0; i<asens_p[d].size(); ++i){
          if(asens_p[d][i]==0){
            cout << "NULL";
          } else {
            cout << *asens_p[d][i];
          }
          cout << ",";
        }
        cout << endl;
      }
#endif
        
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

FX MXFunctionInternal::hessian(int iind, int oind) {
  // Assert initialized
  casadi_assert_message(isInit(),"Function not initialized.");
  
  // Make sure that the function is scalar valued
  //casadi_assert_message(output(oind).numel()==1,"Can only create hessians for scalar valued functions");
  
  // Get the Jacobian of the function
  MX J = jac(iind,oind);
  // MX J = grad(oind).at(iind);
  
  // Construct the gradient function
  MXFunction gfcn(inputv_,trans(J));
  gfcn.init();
  
  // And return the jacobian of that gradient function
  return gfcn.jacobian(iind,0,false,true);
}

void MXFunctionInternal::evalSX(const std::vector<SXMatrix>& input_s, std::vector<SXMatrix>& output_s, 
                                const std::vector<std::vector<SXMatrix> >& fwdSeed, std::vector<std::vector<SXMatrix> >& fwdSens, 
                                const std::vector<std::vector<SXMatrix> >& adjSeed, std::vector<std::vector<SXMatrix> >& adjSens,
                                bool output_given, int offset_begin, int offset_end){
  casadi_assert_message(fwdSens.empty(),"Not implemented");
  casadi_assert_message(adjSeed.empty(),"Not implemented");
  casadi_assert_message(offset_begin==0,"Not implemented");
  casadi_assert_message(offset_end==0,"Not implemented");
  
  // Create a work array
  vector<SXMatrix> swork(work_.size());
  for(vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); it++){
    for(int i=0; i<it->res.size(); ++i){
      if (it->res[i]>=0)
        swork[it->res[i]] = SXMatrix(it->op->sparsity(i));
    }
  }
  
  // Pass the inputs
  for(int ind=0; ind<input_.size(); ++ind){
    swork[input_ind_[ind]].set(input_s[ind]);
  }

  // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
  vector<SXMatrix*> sxarg;
  vector<SXMatrix*> sxres;
  for(vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); it++){
    sxarg.resize(it->arg.size());
    for(int c=0; c<sxarg.size(); ++c){
      int ind = it->arg[c];
      sxarg[c] = ind<0 ? 0 : &swork[ind];
    }
    sxres.resize(it->res.size());
    for(int c=0; c<sxres.size(); ++c){
      int ind = it->res[c];
      sxres[c] = ind<0 ? 0 : &swork[ind];
    }
    it->op->evaluateSX(sxarg,sxres);
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
  evalSX(arg,res,dummy,dummy,dummy,dummy,false);
  
  // Create function
  SXFunction f(arg,res);
  return f;
}

void MXFunctionInternal::collectFree(){
  
  // Mark all the inputs
  for(vector<MX>::iterator i=inputv_.begin(); i!=inputv_.end(); ++i){
    i->setTemp(1);
  }
  
  // Collect free varables
  free_vars_.clear();
  free_vars_ind_.clear();
  for(vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); it++){
    if(it->op->isSymbolic() && it->op.getTemp()!=1){
      free_vars_.push_back(it->op);
      free_vars_ind_.push_back(it->res.front());
    }
  }

  // Unmark inputs
  for(vector<MX>::iterator i=inputv_.begin(); i!=inputv_.end(); ++i){
    i->setTemp(0);
  }
}

FX MXFunctionInternal::getDerivative(int nfwd, int nadj){
  // Forward seeds
  vector<vector<MX> > fseed(nfwd,inputv_);
  for(int dir=0; dir<nfwd; ++dir){
    // Replace symbolic inputs
    for(vector<MX>::iterator i=fseed[dir].begin(); i!=fseed[dir].end(); ++i){
      // Name of the forward seed
      std::stringstream ss;
      ss << "f";
      if(nfwd>1) ss << dir;
      ss << "_";
      i->print(ss);

      // Save to matrix
      *i = MX(ss.str(),i->sparsity());
      
    }
  }
  
  // Adjoint seeds
  vector<vector<MX> > aseed(nadj,outputv_);
  for(int dir=0; dir<nadj; ++dir){
    // Replace symbolic inputs
    int oind=0;
    for(vector<MX>::iterator i=aseed[dir].begin(); i!=aseed[dir].end(); ++i, ++oind){
      // Name of the adjoint seed
      std::stringstream ss;
      ss << "a";
      if(nadj>1) ss << dir << "_";
      ss << oind;

      // Save to matrix
      *i = MX(ss.str(),i->sparsity());
      
    }
  }
  
  // Evaluate symbolically
  vector<vector<MX> > fsens(nfwd,outputv_), asens(nadj,inputv_);
  evalMX(inputv_,outputv_,fseed,fsens,aseed,asens,true);

  // All inputs of the return function
  vector<MX> ret_in;
  ret_in.reserve(inputv_.size()*(1+nfwd) + outputv_.size()*nadj);
  ret_in.insert(ret_in.end(),inputv_.begin(),inputv_.end());
  for(int dir=0; dir<nfwd; ++dir)
    ret_in.insert(ret_in.end(),fseed[dir].begin(),fseed[dir].end());
  for(int dir=0; dir<nadj; ++dir)
    ret_in.insert(ret_in.end(),aseed[dir].begin(),aseed[dir].end());

  // All outputs of the return function
  vector<MX> ret_out;
  ret_out.reserve(outputv_.size()*(1+nfwd) + inputv_.size()*nadj);
  ret_out.insert(ret_out.end(),outputv_.begin(),outputv_.end());
  for(int dir=0; dir<nfwd; ++dir)
    ret_out.insert(ret_out.end(),fsens[dir].begin(),fsens[dir].end());
  for(int dir=0; dir<nadj; ++dir)
    ret_out.insert(ret_out.end(),asens[dir].begin(),asens[dir].end());

  // Assemble function and return
  MXFunction ret(ret_in,ret_out);
  ret.init();
  return ret;
}

} // namespace CasADi

