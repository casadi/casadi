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

#include "inverse_mapping.hpp"
#include "mapping.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../fx/sx_function.hpp"

using namespace std;

namespace CasADi{

InverseMapping::InverseMapping(const MX& dep, const std::vector<CRSSparsity>& sp, const std::vector<int>& nzind, const std::vector<int>& depind) : sp_(sp), nzind_(nzind), depind_(depind){
  setDependencies(dep);
}

InverseMapping* InverseMapping::clone() const{
  return new InverseMapping(*this);
}

void InverseMapping::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  int nadj = adjSeed.size();
  int nfwd = fwdSens.size();
  const vector<double> &inputd = input[0]->data();
  
  // Set outputs and forward sensitivities to zero
  for(int oind=0; oind<output.size(); ++oind){
    fill(output[oind]->begin(),output[oind]->end(),0);
    for(int d=0; d<nfwd; ++d){
      fill(fwdSens[d][oind]->begin(),fwdSens[d][oind]->end(),0);
    }
  }
  
  for(int k=0; k<depind_.size(); ++k){
    output[depind_[k]]->data()[nzind_[k]] += inputd[k];
    
    for(int d=0; d<nfwd; ++d)
      fwdSens[d][depind_[k]]->data()[nzind_[k]] += fwdSeed[d][0]->data()[k];
    
    for(int d=0; d<nadj; ++d)
      adjSens[d][0]->data()[k] += adjSeed[d][depind_[k]]->data()[nzind_[k]];
  }
}

void InverseMapping::propagateSparsity(const DMatrixPtrV& input, DMatrixPtrV& output){
  casadi_assert_message(0,"not implemented");
}

void InverseMapping::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << "inverse_mapping";
}

void InverseMapping::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
  casadi_assert_message(0,"not implemented");
}

void InverseMapping::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  casadi_assert_message(0,"not implemented");
}

int InverseMapping::getNumOutputs() const{
  return sp_.size();
}
    
const CRSSparsity& InverseMapping::sparsity(int oind){
  return sp_.at(oind);
}


} // namespace CasADi
