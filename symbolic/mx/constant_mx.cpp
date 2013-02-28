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

#include "constant_mx.hpp"
#include <cassert>
#include <vector>
#include <algorithm>
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"

using namespace std;

namespace CasADi{

  ConstantMX::ConstantMX(const CRSSparsity& sp){
    setSparsity(sp);
  }

  ConstantMX::~ConstantMX(){
  }

  void ConstantMX::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    int nfwd = fwdSens.size();
    for(int d=0; d<nfwd; ++d){
      fwdSens[d][0]->setZero();
    }
  }

  void ConstantMX::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  }

  void ConstantMX::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    // Evaluate nondifferentiated
    if(!output_given){
      if(output[0]){
	*output[0] = shared_from_this<MX>();
      }
    }
  
    // Number of derivative directions
    int nfwd = fwdSens.size();
    if(nfwd==0) return; // Quick return
  
    // Derivatives
    MX zero_sens = MX::sparse(size1(),size2());
    for(int d=0; d<nfwd; ++d){
      if(fwdSens[d][0]){
	*fwdSens[d][0] = zero_sens;
      }
    }
  }

  void ConstantMX::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    if(fwd){
      bvec_t *outputd = get_bvec_t(output[0]->data());
      fill_n(outputd,output[0]->size(),0);
    }
  }

  void ConstantMX::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    // Print the constant
    int ind = gen.addConstant(shared_from_this<MX>());

    // Copy the constant to the work vector
    stream << "  for(i=0; i<" << sparsity().size() << "; ++i) ";
    stream << res.at(0) << "[i]=";
    stream << "c" << ind << "[i];" << endl;
  }

  bool ConstantMX::__nonzero__() const{
    if (numel()!=1) casadi_error("Can only determine truth value of scalar MX.");
    if (size()!=1) casadi_error("Can only determine truth value of dense scalar MX.");
    return !isZero();
  }

  bool ConstantDMatrix::isZero() const{
    return CasADi::isZero(x_);
  }

  bool ConstantDMatrix::isOne() const{
    return CasADi::isOne(x_);
  }

  bool ConstantDMatrix::isMinusOne() const{
    return CasADi::isMinusOne(x_);
  }

  bool ConstantDMatrix::isIdentity() const{
    return CasADi::isIdentity(x_);
  }

} // namespace CasADi

