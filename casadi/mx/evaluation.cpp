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

// Constructor
Evaluation::Evaluation(const FX& fcn, const vector<MX>& dep, int oind_) : fcn_(fcn), oind(oind_) {
  setDependencies(dep);
  setSparsity(fcn_.output(oind).sparsity());
}

Evaluation* Evaluation::clone() const{
  return new Evaluation(*this);
}

void Evaluation::print(ostream &stream) const{
  stream << fcn_ << "[" << dep() << "]";
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

  // Pass the adjoint seed to the function
  for(int d=0; d<nadj; ++d)
    fcn_.setAdjSeed(adjSeed[d],oind,d);
    
  // Set adjoint seed to zero for all other outputs in all other directions
  for(int ind=0; ind<fcn_.getNumOutputs(); ++ind){
    if(ind!=oind){
      for(int d=0; d<nadj; ++d){
        fill(fcn_.adjSeed(ind,d).begin(),fcn_.adjSeed(ind,d).end(),0);
      }
    }
  }

  // Evaluate
  fcn_.evaluate(nfwd>0, nadj>0);
  
  // Get the results
  fcn_.getOutput(output,oind);

  // Get the fwd sensitivities
  for(int d=0; d<nadj; ++d)
    fcn_.getFwdSens(fwdSens[d],oind,d);
  
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

} // namespace CasADi
