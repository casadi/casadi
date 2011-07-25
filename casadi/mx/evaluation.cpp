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

#include "evaluation.hpp"
#include "../fx/fx_internal.hpp"
#include "../stl_vector_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../fx/x_function.hpp"
#include "jacobian_reference.hpp"

using namespace std;

namespace CasADi{

Evaluation::Evaluation(const FX& fcn, const vector<MX>& dep) : fcn_(fcn) {
  // Argument checking
  if (dep.size()!=fcn.getNumInputs()) {
    std::stringstream s;
    s << "Evaluation::Evaluation: number of passed-in dependencies (" << dep.size() << ") should match number of inputs of function (" << fcn.getNumInputs() << ").";
    throw CasadiException(s.str());
  }
  // Assumes initialised
  for (int i=0;i<fcn.getNumInputs();i++) {
     if (dep[i].isNull())
       continue;
      if (dep[i].size1()!=fcn.input(i).size1() || dep[i].size2()!=fcn.input(i).size2()) {
        std::stringstream s;
        s << "Evaluation::shapes of passed-in dependencies should match shapes of inputs of function." << std::endl;
        s << "Input argument " << i << " has shape (" << fcn.input(i).size1() << "," << fcn.input(i).size2() << ") while a shape (" << dep[i].size1() << "," << dep[i].size2() << ") was supplied.";
        throw CasadiException(s.str());
      }     
   }
  setDependencies(dep);
  setSparsity(CRSSparsity(1,1,true));
}

Evaluation* Evaluation::clone() const{
  return new Evaluation(*this);
}

void Evaluation::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << fcn_ << ".call(" << args << ")";
}

void Evaluation::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();
  
  // Pass the input and forward seeds to the function
  for(int i=0; i<ndep(); ++i){
    if(input[i] != 0 && input[i]->size() !=0 ){
      fcn_.setInput(input[i]->data(),i);
      for(int d=0; d<nfwd; ++d){
        fcn_.setFwdSeed(fwdSeed[d][i]->data(),i,d);
      }
    }
  }
  
  // Pass the adjoint seed to the function
  for(int i=0; i<getNumOutputs(); ++i){
    for(int d=0; d<nadj; ++d){
      if(adjSeed[d][0]!=0 && adjSeed[d][0]->size() != 0){
        fcn_.setAdjSeed(adjSeed[d][0]->data(),i,d);
      }
    }
  }

  // Evaluate
  fcn_.evaluate(nfwd, nadj);
  
  // Get the outputs and forward sensitivities
  for(int i=0; i<getNumOutputs(); ++i){
    if(output[i] != 0 && output[i]->size() !=0 ){
      fcn_.getOutput(output[i]->data(),i);
      for(int d=0; d<nfwd; ++d){
        fcn_.getFwdSens(fwdSens[d][i]->data(),i,d);
      }
    }
  }
  
  // Get the adjoint sensitivities
  for(int i=0; i<ndep(); ++i){
    for(int d=0; d<nadj; ++d){
      if(adjSens[d][i] != 0 && adjSens[d][i]->size() != 0){
        const vector<double>& asens = fcn_.adjSens(i,d).data();
        for(int j=0; j<asens.size(); ++j)
          adjSens[d][i]->data()[j] += asens[j];
      }
    }
  }
}


EvaluationOutput::EvaluationOutput(const MX& parent, int oind) : OutputNode(parent,oind){  
  // Save the sparsity pattern
  setSparsity(getFunction().output(oind).sparsity());
}

EvaluationOutput* EvaluationOutput::clone() const{
  return new EvaluationOutput(*this);
}

void EvaluationOutput::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << args[0] << "[" << oind_ <<  "]";
}

FX& Evaluation::getFunction(){ 
  return fcn_;
}

FX& EvaluationOutput::getFunction(){ 
  return dep(0)->getFunction();
}

void Evaluation::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  // Make sure that the function is an X-function
  XFunction fcn = shared_cast<XFunction>(fcn_);
  casadi_assert_message(!fcn.isNull(),"Function not an SXFunction or MXFunction");
  vector<SXMatrix> arg(input.size());
  for(int i=0; i<arg.size(); ++i){
    arg[i] = *input[i];
  }
  xs_ = fcn.eval(arg);
}

void Evaluation::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();
  
  if(!output_given){
    // Evaluate the function symbolically
    vector<MX> arg(input.size());
    for(int i=0; i<arg.size(); ++i){
      if(input[i])
        arg[i] = *input[i];
    }
    vector<MX> res = fcn_.call(arg);
    for(int i=0; i<res.size(); ++i){
      if(output[i])
        *output[i] = res[i];
    }
  }
  
  if(nfwd>0){
    for(int oind=0; oind<output.size(); ++oind){
      // Skip of not used
      if(output[oind]==0) continue;
          
      // Return matrix
      MX ret2 = MX::zeros(output[oind]->numel(),1);
      vector<MX> ret3(nfwd,ret2);
      
      for(int iind=0; iind<fwdSeed.front().size(); ++iind){
        if(fwdSeed[0][iind]!=0){
          MX J = MX::create(new JacobianReference(*output[oind],iind));
          for(int d=0; d<nfwd; ++d){
            ret3[d] += prod(J,vec(*fwdSeed[d][iind]));
          }
        }
      }
      
      for(int d=0; d<nfwd; ++d){
        *fwdSens[d][oind] = reshape(ret3[d],output[oind]->size1(),output[oind]->size2());
      }
    }
  }
  
/*  for(int d=0; d<nfwd; ++d){
    for(int oind=0; oind<output.size(); ++oind){
      if(fwdSens[d][oind]!=0 && !output[oind]->isNull()){
        continue;
        *fwdSens[d][oind] = MX::zeros(output[oind]->size1(),output[oind]->size2());
        for(int iind=0; iind<input.size(); ++iind){
          if(fwdSeed[d][iind]!=0){
            MX J = MX::create(new JacobianReference(*output[oind],iind));
            *fwdSens[d][oind] += prod(J,*fwdSeed[d][iind]);
          }
        }
      }
    }
  }*/
}

void EvaluationOutput::evaluateSX(const std::vector<SXMatrix*> &input, SXMatrix& output){
  // Get a reference the arguments
  const vector<SXMatrix>& xs = dynamic_cast<Evaluation*>(dep(0).get())->xs_;
  
  // Copy to output
  output.set(xs[oind_]);
}

void Evaluation::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
  MXNode::deepCopyMembers(already_copied);
  fcn_ = deepcopy(fcn_,already_copied);
}

} // namespace CasADi
