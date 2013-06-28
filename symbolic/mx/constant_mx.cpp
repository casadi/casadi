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
    for(int d=0; d<fwdSens.size(); ++d){
      fwdSens[d][0]->setZero();
    }
    for(int d=0; d<adjSeed.size(); ++d){
      adjSeed[d][0]->setZero();
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
  
    // Forward sensitivities
    if(!fwdSens.empty()){
      MX zero_sens = MX::sparse(size1(),size2());
      for(int d=0; d<fwdSens.size(); ++d){
        if(fwdSens[d][0]!=0){
          *fwdSens[d][0] = zero_sens;
        }
      }
    }
    
    // Clear adjoint seeds
    for(int d=0; d<adjSeed.size(); ++d){
      if(adjSeed[d][0]!=0){
        *adjSeed[d][0] = MX();
      }      
    }
  }

  void ConstantMX::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    bvec_t *outputd = get_bvec_t(output[0]->data());
    fill_n(outputd,output[0]->size(),0);
  }

  void ConstantDMatrix::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    // Print the constant
    int ind = gen.getConstant(x_.data(),true);

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

  ConstantMX* ConstantMX::create(const CRSSparsity& sp, int val){
    switch(val){
    case 0: return new Constant<CompiletimeConst<0> >(sp);
    case 1: return new Constant<CompiletimeConst<1> >(sp);
    case -1: return new Constant<CompiletimeConst<(-1)> >(sp);
    default: return new Constant<RuntimeConst<int> >(sp,val);        
    }
  }

  ConstantMX* ConstantMX::create(const CRSSparsity& sp, double val){
    int intval(val);
    if(intval-val==0){
      return create(sp,intval);
    } else {
      return new Constant<RuntimeConst<double> >(sp,val);
    }
  }

  ConstantMX* ConstantMX::create(const Matrix<double>& val){
    if(val.size()==0){
      return create(val.sparsity(),0);
    } else if(val.scalar()){
      return create(val.sparsity(),val.toScalar());
    } else {
      // Check if all values are the same
      const vector<double> vdata = val.data();
      double v = vdata[0];
      for(vector<double>::const_iterator i=vdata.begin(); i!=vdata.end(); ++i){
        if(*i!=v){
          // Values not all the same
          return new ConstantDMatrix(val);
        }
      }

      // All values identical if reached this point
      return create(val.sparsity(),v);
    }
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

  MX ConstantMX::getMultiplication(const MX& y) const{
    if(y.isConstant()){
      // Constant folding
      DMatrix xv = getMatrixValue();
      DMatrix yv = y->getMatrixValue();
      return mul(xv,yv);
    } else {
      return MXNode::getMultiplication(y);
    }
  }

  MX ConstantMX::getInnerProd(const MX& y) const{
    if(y.isConstant()){
      // Constant folding
      DMatrix xv = getMatrixValue();
      DMatrix yv = y->getMatrixValue();
      return inner_prod(xv,yv);
    } else {
      return MXNode::getInnerProd(y);
    }
  }

  bool ConstantDMatrix::isEqual(const MXNode* node, int depth) const{
    // Check if same node
    const ConstantDMatrix* n = dynamic_cast<const ConstantDMatrix*>(node);
    if(n==0) return false;
    
    // Check sparsity
    if(this->sparsity()!=node->sparsity()) return false;
    
    // Check indices
    if(!std::equal(x_.begin(),x_.end(),n->x_.begin())) return false;
    
    return true;
  }

} // namespace CasADi

