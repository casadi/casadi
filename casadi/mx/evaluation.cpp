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

void Evaluation::evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  // Pass the input and forward seeds to the function
  for(int i=0; i<ndep(); ++i){
    if(input[i] != 0){
      fcn_.setInput(input[i],i);
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
  
  // Get the function
  const Evaluation* p = dynamic_cast<const Evaluation*>(parent.get());
  casadi_assert(p!=0);
  fcn_ = p->fcn_;

  // Save the sparsity pattern
  setSparsity(fcn_.output(oind).sparsity());
}

EvaluationOutput* EvaluationOutput::clone() const{
  return new EvaluationOutput(*this);
}

void EvaluationOutput::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << args[0] << "[" << oind_ <<  "]";
}

void EvaluationOutput::evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  // Pass the adjoint seed to the function
  for(int d=0; d<nadj; ++d)
    if(adjSeed[d]!=0)
      fcn_.setAdjSeed(adjSeed[d],oind_,d);

    // Get the results
  fcn_.getOutput(output,oind_);

  // Get the fwd sensitivities
  for(int d=0; d<nfwd; ++d)
    if(fwdSens[d]!=0)
      fcn_.getFwdSens(fwdSens[d],oind_,d);
}

MX Evaluation::adFwd(const std::vector<MX>& jx){ 
  // Save the forward derivative
  x_ = jx;
  
  // Return null
  return MX();
}

MX EvaluationOutput::adFwd(const std::vector<MX>& jx){
  cout << "called EvaluationOutput::adFwd " << endl;
  
  
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
    if(0){
      if(!x[i].isNull()){
        ret += prod(jac(i),x[i]);
      }
    } else {
      
      // Get the Jacobian (this is inefficient, unless jacobian gets a bit smarter!)
      FX J = fcn_.jacobian(i,oind_);
      
      // If the jacobian is not zero
      if(!J.isNull()){
        J.init();

        // Get a reference to the argument
        vector<MX> &Jarg = dep(0)->dep_;

        // Create an evaluation node
        MX Ji = J.call(Jarg).at(0);
        
        // Assemble the return matrix
        if(!x[i].isNull()){
          ret += prod(Ji,x[i]);
        }
      }
    }
  }
  return ret;
}


MX EvaluationOutput::jac(int iind){
  return MX::create(new JacobianReference(MX::create(this),iind));
}



} // namespace CasADi
