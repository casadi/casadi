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

using namespace std;

namespace CasADi{

Evaluation::Evaluation(const FX& fcn, const vector<MX>& dep) : fcn_(fcn) {
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
  fcn_.evaluate(nfwd>0, nadj>0);
  
  // Get the adjoint sensitivities
  for(int i=0; i<ndep(); ++i){
    for(int d=0; d<nadj; ++d){
      if(adjSens[i][d] != 0){
        const vector<double>& asens = fcn_.adjSens(i,d);
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

} // namespace CasADi
