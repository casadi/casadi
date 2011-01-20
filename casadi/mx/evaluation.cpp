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
  setSize(fcn_->output(oind).size1(),fcn_->output(oind).size2());
}

Evaluation* Evaluation::clone() const{
  return new Evaluation(*this);
}

void Evaluation::print(ostream &stream) const{
  stream << fcn_ << "[" << dep() << "]";
}

void Evaluation::evaluate(int fsens_order, int asens_order){
  // Pass the input to the function
  for(int i=0; i<ndep(); ++i){
    if(!dep(i).isNull()){
      fcn_.setInput(input(i),i);
    }
  }

  // Give the forward seed to the function
  if(fsens_order>0){
    for(int i=0; i<ndep(); ++i){
      if(!dep(i).isNull()){
        fcn_.setFwdSeed(fwdSeed(i),i);
      }
    }
  }

  // Pass the adjoint seed to the function
  if(asens_order>0){
    fcn_.setAdjSeed(adjSeed(),oind);
  }

  // Evaluate
  fcn_.evaluate(fsens_order, asens_order);
  
  // Get the results
  fcn_.getOutput(output(),oind);

  // Fwd sens
  if(fsens_order>0){
    fcn_.getFwdSens(fwdSens(),oind);
  }

  // Adjoint sens
  if(asens_order>0)
    for(int i=0; i<ndep(); ++i)
      if(!dep(i).isNull())
        fcn_.getAdjSens(adjSens(i),i);
}

} // namespace CasADi
