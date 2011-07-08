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

#include <stack>
#include <typeinfo>

using namespace std;

namespace CasADi{

MXFunctionInternal::MXFunctionInternal(const std::vector<MX>& inputv_, const std::vector<MX>& outputv_) : inputv(inputv_), outputv(outputv_){
  setOption("name", "unnamed_mx_function");

  liftfun_ = 0;
  liftfun_ud_ = 0;

  for(int i=0; i<inputv_.size(); ++i) {
    if (inputv_[i].isNull()) {
      stringstream ss;
      ss << "MXFunctionInternal::MXFunctionInternal: MXfunction input arguments cannot be null." << endl;
      ss << "Argument #" << i << " is null." << endl;
      throw CasadiException(ss.str());
    } else if (!inputv_[i]->isSymbolic()) {
      stringstream ss;
      ss << "MXFunctionInternal::MXFunctionInternal: MXfunction input arguments must be purely symbolic." << endl;
      ss << "Argument #" << i << " is not symbolic." << endl;
      throw CasadiException(ss.str());
    }
  }
  
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
  for(vector<MX>::iterator it = inputv.begin(); it!=inputv.end(); ++it)
    if(it->get()) s.push(static_cast<MXNode*>(it->get()));

  // Add the outputs to the stack
  for(vector<MX>::iterator it = outputv.begin(); it!=outputv.end(); ++it)
    if(it->get()) s.push(static_cast<MXNode*>(it->get()));

  // All evaluation nodes in the order of evaluation
  vector<MXNode*> nodes;

  // Order the nodes in the order of dependencies using a depth-first topological sorting
  sort_depth_first(s,nodes);
  
  // TODO: Make sure that the output nodes are placed directly after the corresponding multiple output node

  // Resort the nodes in a more cache friendly order (Kahn 1962)
  //  resort_bredth_first(nodes);
  
  // Set the temporary variables to be the corresponding place in the sorted graph
  for(int i=0; i<nodes.size(); ++i)
    nodes[i]->temp = i;
  
  // Indices corresponding to the inputs
  inputv_ind.resize(inputv.size());
  for(int i=0; i<inputv_ind.size(); ++i)
    inputv_ind[i] = inputv[i]->temp;

  // Indices corresponding to the outputs
  outputv_ind.resize(outputv.size());
  for(int i=0; i<outputv_ind.size(); ++i)
    outputv_ind[i] = outputv[i]->temp;
  
  // Create runtime elements for each node
  alg.resize(nodes.size());
  for(int ii=0; ii<alg.size(); ++ii){
    MX &m = alg[ii].mx;
    AlgEl* it = &alg[ii];
    
    // Save the node
    it->mx.assignNode(nodes[ii]);
    it->val.data = Matrix<double>(m->sparsity(),0);
    it->val.dataF.resize(nfdir_);
    it->val.dataA.resize(nadir_);
    it->val.init();

    // Save the indices of the children nodes
    it->ch.resize(m->ndep());
    for(int i=0; i<m->ndep(); ++i){
      if(m->dep(i).isNull())
        it->ch[i] = -1;
      else
        it->ch[i] = m->dep(i)->temp;
    }
    
    it->input.resize(m->ndep(),0);
    it->fwdSeed.resize(m->ndep());
    it->adjSens.resize(m->ndep());
    for(int i=0; i<m->ndep(); ++i){
      it->fwdSeed[i].resize(nfdir_,0);
      it->adjSens[i].resize(nadir_,0);
      if(it->ch[i]>=0){
		  it->input[i] = &alg[it->ch[i]].val.data;
        for(int d=0; d<nfdir_; ++d)
          it->fwdSeed[i][d] = &alg[it->ch[i]].val.dataF[d];
        for(int d=0; d<nadir_; ++d)
          it->adjSens[i][d] = &alg[it->ch[i]].val.dataA[d];
      }
    }

    it->fwdSens.resize(nfdir_);
    for(int d=0; d<nfdir_; ++d)
      it->fwdSens[d] = &it->val.dataF[d];
    it->adjSeed.resize(nadir_);
    for(int d=0; d<nadir_; ++d)
      it->adjSeed[d] = &it->val.dataA[d];
    
  }
  
  // Reset the temporary variables
  for(int i=0; i<nodes.size(); ++i){
    nodes[i]->temp = 0;
  }

  // Eliminate the jacobian nodes from the algorithm
  eliminateJacobian();

log("MXFunctionInternal::init end");
}

void MXFunctionInternal::eliminateJacobian(){
  // An output/input pair corresponding to a jacobian block
  typedef pair<int,int> oipair;
  
  // The set of jacobian blocks for a single evaluation node
  typedef set<oipair> jblocks;
  
  // The sets of jacobian blocks for each evaluation node
  typedef map<void*,jblocks> jbmap;
  
  // The set of evaluation nodes for each function
  typedef map<void*,jbmap> evmap;
  
  // Find out which functions need to be differentiated with respect to which arguments
  evmap m;
  for(vector<AlgEl>::iterator it=alg.begin(); it!=alg.end(); it++){
    // Check if the node is a Jacobian reference
    
    //JacobianReference* jref = dynamic_cast<JacobianReference*>(it->mx.get());
    if(it->mx->isJacobian()){
      // Get a reference to the function
      FX& fcn = it->mx->getFunction();

      // Get the input and output indices
      int iind = it->mx->getFunctionInput();
      int oind = it->mx->getFunctionOutput();

      // Save to map
      m[fcn.get()][it->mx->dep(0)->dep(0).get()].insert(oipair(oind,iind));
    }
  }

  // Quick end if no Jacobian references
  if(m.empty()) return;

  // Add function output to each of the function outputs
  for(evmap::iterator i=m.begin(); i!=m.end(); ++i){
    // Get number of outputs
    FXInternal* f = reinterpret_cast<FXInternal*>(i->first);
    int n_out = f->getNumOutputs();
  
    // For each evaluation node...
    for(jbmap::iterator j=i->second.begin(); j!=i->second.end(); ++j){
      // .. add undifferentiated functions to desired output
    
      for(int k=0; k<n_out; ++k){
        j->second.insert(oipair(k,-1));
      }
    }
  }

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
    void *fcn_ref=0, *eval_ref=0;
    if(it->mx->isEvaluation()){
      fcn_ref = it->mx->getFunction().get();
      eval_ref = it->mx.get();
    } else if(it->mx->isEvaluationOutput()){
      fcn_ref = it->mx->getFunction().get();
      eval_ref = it->mx->dep(0).get();
    } else if(it->mx->isJacobian()){
      fcn_ref = it->mx->getFunction().get();
      eval_ref = it->mx->dep(0)->dep(0).get();
    }
    
    // Update evaluation/evaluation output/jacobian node
    if(eval_ref || fcn_ref){
    
      // Check if the function should be replaced
      evmap::const_iterator ii=m.find(fcn_ref);
      if(ii!=m.end()){
        
        // Check if the function evaluation is to be replaced
        jbmap::const_iterator jj=ii->second.find(eval_ref);
        if(jj!=ii->second.end()){
          
          // Make the pointer unique so that other expressions won't be affected
          it->mx.makeUnique();
            
          if(it->mx->isEvaluation()){
            // Replace the function if its an evaluation node
            it->mx->getFunction() = jac[fcn_ref][jj->second];
              
          } else {
            // Jacobian reference or evaluation output
            casadi_assert(it->mx->isJacobian() || it->mx->isEvaluationOutput());
              
            // Get current input and output indices
            int iind = it->mx->getFunctionInput();
            int oind = it->mx->getFunctionOutput();
              
            // Update the output index if we have an evaluation node
            int oind_new = 0;
            for(jblocks::const_iterator kk=jj->second.begin(); kk!=jj->second.end(); ++kk){
              // Check if the input and output indices matches
              if(kk->first == oind && kk->second==iind){
                // Update the dependency index if its a Jacobian reference
                if(it->mx->isJacobian()){
                  it->ch[0] = alg[it->ch[0]].ch[0];
                }
                
                // Create a new evaluation output
                it->mx = MX::create(new EvaluationOutput(alg[it->ch[0]].mx,oind_new));
                break;
              }
              
              // Update index
              oind_new++;
            }
            casadi_assert(!it->mx->isJacobian());
          }
        }
      } // if(ii!=m.end()
    }
    casadi_assert(!it->mx->isJacobian());
    
    // Check if any of the dependencies have changed
    for(int i=0; i<it->ch.size(); ++i){
      int ch = it->ch[i];
      if(ch>=0){
        if(it->mx->dep(i).get() != alg[ch].mx.get()){
          // Make sure that the change does not affect other nodes
          it->mx.makeUnique();
            
          // Update dependencies
          it->mx->dep(i) = alg[ch].mx;
        }
      }
    }
  }
}

void MXFunctionInternal::setLiftingFunction(LiftingFunction liftfun, void* user_data){
  liftfun_ = liftfun;
  liftfun_ud_ = user_data;
}

void MXFunctionInternal::evaluate(int nfdir, int nadir){
  log("MXFunctionInternal::evaluate begin");
  
  // Pass the inputs
  for(int ind=0; ind<input_.size(); ++ind)
    alg[inputv_ind[ind]].val.data.set(input(ind));

  // Pass the forward seeds
  for(int dir=0; dir<nfdir; ++dir)
    for(int ind=0; ind<input_.size(); ++ind)
      alg[inputv_ind[ind]].val.dataF.at(dir).set(fwdSeed(ind,dir));
  
  // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
  for(vector<AlgEl>::iterator it=alg.begin(); it!=alg.end(); it++){
    it->mx->evaluate(it->input, it->val.data, it->fwdSeed, it->fwdSens, it->adjSeed, it->adjSens, nfdir, 0);
    // Lifting
    if(liftfun_ && it->mx->isNonLinear()){
      liftfun_(&it->val.data.front(),it->mx.size(),liftfun_ud_);
    }
  }
  
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
      it->mx->evaluate(it->input, it->val.data, it->fwdSeed, it->fwdSens, it->adjSeed, it->adjSens, 0, nadir);
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
  MXFunctionInternal* node = new MXFunctionInternal(inputv,outputv);
  node->setOption(dictionary());
  if(isInit()) node->init();
  return node;
}

/*CRSSparsity MXFunctionInternal::getJacSparsity(int iind, int oind){
  std::vector<MX> ret = jac(iind);
  return ret.at(oind).sparsity();
}*/
  
std::vector<MX> MXFunctionInternal::jac(int iind){
  casadi_assert(isInit());
  
  // Number of columns in the Jacobian
  int ncol = alg[inputv_ind[iind]].mx.size1();
  
  // Get theseed matrices
  vector<MX> fseed(input_.size());
  for(int ind=0; ind<input_.size(); ++ind){
    // Number of rows in the seed matrix
    int nrow = alg[inputv_ind[ind]].mx.size1();
    
    // Create seed matrix
    fseed[ind] = ind==iind ? MX::eye(ncol) : MX::zeros(nrow,ncol);
  }
  
  // Forward mode automatic differentiation, symbolically
  return adFwd(fseed);
}

std::vector<MX> MXFunctionInternal::adFwd(const std::vector<MX>& fseed){
  casadi_assert(isInit());
  
  // Get the number of columns
  const int ncol = fseed.front().size2();
  
  // Make sure that the number of columns is consistent
  for(vector<MX>::const_iterator it=fseed.begin(); it!=fseed.end(); ++it){
    casadi_assert_message(ncol==it->size2(),"Number of columns in seed matrices not consistent.");
  }
  
  // Directional derivative for each node
  std::vector<MX> derwork(alg.size());
  
  // Pass the seed matrices for the symbolic variables
  for(int ind=0; ind<input_.size(); ++ind){
    derwork[inputv_ind[ind]] = fseed[ind];
  }
  
  // Evaluate all the seed matrices of the algorithm sequentially
  for(int el=0; el<derwork.size(); ++el){
    
    // Skip the node if it has already been calculated (i.e. is symbolic or constant)
    if(!derwork[el].isNull()) continue;
    
    // Zero seed matrix if constant
    if(alg[el].mx->isConstant()){
      
      // Number of rows
      int nrow = alg[el].mx.numel();
      
      // Save zero matrix
      derwork[el] = MX::zeros(nrow,ncol);
      continue;
    }
    
    // Collect the seed matrices
    vector<MX> seed(alg[el].ch.size());
    casadi_assert(!seed.empty());
    for(int i=0; i<seed.size(); ++i){
      int ind = alg[el].ch[i];
      if(ind>=0)
        seed[i] = derwork[ind];
    }

    // Get the sensitivity matrix
    MX sens = alg[el].mx->adFwd(seed);
    
    // Save to the memory vector
    derwork[el] = sens;
  }

  // Collect the symbolic forward sensitivities
  vector<MX> ret(outputv.size());
  for(int ind=0; ind<outputv.size(); ++ind){
    ret[ind] = derwork[outputv_ind[ind]];
  }
  
  return ret;
}

FX MXFunctionInternal::hessian(int iind, int oind) {
  // Assert initialized
  casadi_assert_message(isInit(),"Function not initialized.");
  
  // Make sure that the function is scalar valued
  //casadi_assert_message(output(oind).numel()==1,"Can only create hessians for scalar valued functions");
  
  // Get the Jacobian of the function
  MX J = jac(iind).at(oind);

  // Construct the gradient function
  MXFunction gfcn(inputv,trans(J));
  gfcn.init();
  
  // And return the jacobian of that gradient function
  return gfcn.jacobian(iind,0);
}

void MXFunctionInternal::evaluateSX(const std::vector<Matrix<SX> >& input_s, std::vector<Matrix<SX> >& output_s, bool eliminate_constants){
  // Create a work array
  vector<SXMatrix> work(alg.size());
  for(int i=0; i<work.size(); ++i){
    work[i] = SXMatrix(alg[i].mx.sparsity());
  }
  
  // Pass the inputs
  for(int ind=0; ind<input_.size(); ++ind)
    work[inputv_ind[ind]].set(input_s[ind]);

  // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
  vector<SXMatrix*> d;
  for(int i=0; i<alg.size(); ++i){
    d.resize((alg[i].ch.size()));
    for(int c=0; c<d.size(); ++c){
      d[c] = &work[alg[i].ch[c]];
    }
    alg[i].mx->evaluateSX(d,work[i]);
  }
  
  // Get the outputs
  for(int ind=0; ind<outputv.size(); ++ind)
    work[outputv_ind[ind]].get(output_s[ind]);
}

SXFunction MXFunctionInternal::expand(const std::vector<SXMatrix>& inputv_sx ){
  casadi_assert(isInit());
  
  // Create inputs with the same name and sparsity as the matrix valued symbolic inputs
  vector<SXMatrix> arg(inputv.size());
  if(inputv_sx.empty()){ // No symbolic input provided
    for(int i=0; i<arg.size(); ++i){
      arg[i] = symbolic(inputv[i]->getName(),inputv[i].sparsity());
    }
  } else { // Use provided symbolic input
    // Make sure number of inputs matches
    casadi_assert(inputv_sx.size()==inputv.size());
    
    // Make sure that sparsity matches
    for(int i=0; i<inputv_sx.size(); ++i){
      casadi_assert(inputv_sx[i].sparsity() == inputv[i].sparsity());
    }

    // Copy to argument vector
    copy(inputv_sx.begin(),inputv_sx.end(),arg.begin());
  }

  // Create output vector with correct sparsity
  vector<SXMatrix> res(outputv.size());
  for(int i=0; i<res.size(); ++i){
    res[i] = SXMatrix(outputv[i].sparsity());
  }
  
  // Evaluate symbolically
  evaluateSX(arg,res);
  
  // Create function
  SXFunction f(arg,res);
  return f;
}



} // namespace CasADi

