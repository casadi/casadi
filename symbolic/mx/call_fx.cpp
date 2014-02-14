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

#include "call_fx.hpp"
#include "../fx/fx_internal.hpp"
#include "../stl_vector_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../matrix/matrix_tools.hpp"

using namespace std;

namespace CasADi {

  CallFX::CallFX(const FX& fcn, std::vector<MX> arg) : fcn_(fcn) {

    // Number inputs and outputs
    int num_in = fcn.getNumInputs();
    casadi_assert_message(arg.size()<=num_in,"Argument list length (" << arg.size() << ") exceeds number of inputs (" << num_in << ")");

    // Add arguments if needed
    arg.resize(num_in);

    // Replace nulls with zeros of the right dimension
    for(int i=0; i<arg.size(); ++i){
      if(arg[i].isNull()) arg[i] = MX::zeros(fcn_.input(i).sparsity());
    }

    setDependencies(arg);
    setSparsity(CRSSparsity(1, 1, true));
  }

  CallFX* CallFX::clone() const {
    return new CallFX(*this);
  }

  void CallFX::printPart(std::ostream &stream, int part) const {
    fcn_->printPart(this,stream,part);
  }

  void CallFX::evaluateD(const DMatrixPtrV& arg, DMatrixPtrV& res, std::vector<int>& itmp, std::vector<double>& rtmp) {
    fcn_->evaluateD(this,arg,res,itmp,rtmp);
  }

  int CallFX::getNumOutputs() const {
    return fcn_.getNumOutputs();
  }

  const CRSSparsity& CallFX::sparsity(int oind) const{
    return fcn_.output(oind).sparsity();
  }

  FX& CallFX::getFunction() {
    return fcn_;
  }

  void CallFX::evaluateSX(const SXMatrixPtrV& arg, SXMatrixPtrV& res, std::vector<int>& itmp, std::vector<SX>& rtmp) {
    fcn_->evaluateSX(this,arg,res,itmp,rtmp);
  }

  void CallFX::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given) {
    fcn_->evaluateMX(this,input,output,fwdSeed,fwdSens,adjSeed,adjSens,output_given);
  }

  void CallFX::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    MXNode::deepCopyMembers(already_copied);
    fcn_ = deepcopy(fcn_, already_copied);
  }

  void CallFX::propagateSparsity(DMatrixPtrV& arg, DMatrixPtrV& res, std::vector<int>& itmp, std::vector<double>& rtmp, bool use_fwd) {
    fcn_->propagateSparsity(this,arg,res,itmp,rtmp,use_fwd);
  }

  void CallFX::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    fcn_->generateOperation(this,stream,arg,res,gen);
  }
  
  void CallFX::nTmp(size_t& ni, size_t& nr){
    fcn_->nTmp(this,ni,nr);
  }

} // namespace CasADi
