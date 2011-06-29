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
#include "jacobian_reference.hpp"
#include "../fx/fx_internal.hpp"
#include "../stl_vector_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../fx/x_function.hpp"

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

void Evaluation::evaluate(const std::vector<DMatrix*>& input, DMatrix& output, const VVDptr& fwdSeed, std::vector<DMatrix*>& fwdSens, const std::vector<DMatrix*>& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  
  // Pass the input and forward seeds to the function
  for(int i=0; i<ndep(); ++i){
    if(input[i] != 0 && input[i]->size() !=0 ){
      fcn_.setInput(input[i]->data(),i);
      for(int d=0; d<nfwd; ++d){
        fcn_.setFwdSeed(fwdSeed[i][d],i,d);
      }
    }
  }

  // Evaluate
  fcn_.evaluate(nfwd, nadj);
  
  // Get the adjoint sensitivities
  for(int i=0; i<ndep(); ++i){
    for(int d=0; d<nadj; ++d){
      if(adjSens[i][d] != 0){
        const vector<double>& asens = fcn_.adjSens(i,d).data();
        for(int j=0; j<asens.size(); ++j)
          adjSens[i][d][j] += asens[j];
      }
    }
  }
}


EvaluationOutput::EvaluationOutput(const MX& parent, int oind) : OutputNode(parent), oind_(oind){
  setDependencies(parent);
  
  // Save the sparsity pattern
  setSparsity(getFunction().output(oind).sparsity());
}

EvaluationOutput* EvaluationOutput::clone() const{
  return new EvaluationOutput(*this);
}

void EvaluationOutput::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << args[0] << "[" << oind_ <<  "]";
}

void EvaluationOutput::evaluate(const std::vector<DMatrix*>& input, DMatrix& output, const VVDptr& fwdSeed, std::vector<DMatrix*>& fwdSens, const std::vector<DMatrix*>& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  vector<double> &outputd = output.data();

  // Pass the adjoint seed to the function
  for(int d=0; d<nadj; ++d)
    if(adjSeed[d]!=0 && adjSeed[d]->size() != 0)
      getFunction().setAdjSeed(adjSeed[d]->data(),oind_,d);

    // Get the results
  getFunction().getOutput(outputd,oind_);

  // Get the fwd sensitivities
  for(int d=0; d<nfwd; ++d)
    if(fwdSens[d]!=0 && fwdSens[d]->size() != 0)
      getFunction().getFwdSens(fwdSens[d]->data(),oind_,d);
}

MX Evaluation::adFwd(const std::vector<MX>& jx){ 
  // Save the forward derivative
  x_ = jx;
  
  // Return null
  return MX();
}

MX EvaluationOutput::adFwd(const std::vector<MX>& jx){
  
  // Get a reference the arguments
  vector<MX>& x = dynamic_cast<Evaluation*>(dep(0).get())->x_;
  
  // Find the number of columns
  int ncol = -1;
  for(int i=0; i<x.size(); ++i){
    if(!x[i].isNull())
      ncol = x[i].size2();
  }
  casadi_assert(ncol>=0);
  
  // Return matrix
  MX ret = MX::zeros(size(),ncol);
  for(int i=0; i<x.size(); ++i){
    if(!x[i].isNull()){
      ret += prod(jac(i),x[i]);
    }
  }
  return ret;
}


MX EvaluationOutput::jac(int iind){
  return MX::create(new JacobianReference(MX::create(this),iind));
}

FX& Evaluation::getFunction(){ 
  return fcn_;
}

FX& EvaluationOutput::getFunction(){ 
  return dep(0)->getFunction();
}

void Evaluation::evaluateSX(const std::vector<SXMatrix*> &input, SXMatrix& output){
  // Make sure that the function is an X-function
  XFunction fcn = shared_cast<XFunction>(fcn_);
  casadi_assert_message(!fcn.isNull(),"Function not an SXFunction or MXFunction");
  vector<SXMatrix> arg(input.size());
  for(int i=0; i<arg.size(); ++i){
    arg[i] = *input[i];
  }
  xs_ = fcn.eval(arg);
}

void EvaluationOutput::evaluateSX(const std::vector<SXMatrix*> &input, SXMatrix& output){
  // Get a reference the arguments
  const vector<SXMatrix>& xs = dynamic_cast<Evaluation*>(dep(0).get())->xs_;
  
  // Copy to output
  output.set(xs[oind_]);
}


} // namespace CasADi
