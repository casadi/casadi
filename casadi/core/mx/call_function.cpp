/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#include "call_function.hpp"
#include "../function/function_internal.hpp"
#include "../std_vector_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../matrix/matrix_tools.hpp"

using namespace std;

namespace casadi {

  CallFunction::CallFunction(const Function& fcn, std::vector<MX> arg) : fcn_(fcn) {

    // Number inputs and outputs
    int num_in = fcn.getNumInputs();
    casadi_assert_message(arg.size()<=num_in, "Argument list length (" << arg.size()
                          << ") exceeds number of inputs (" << num_in << ")");

    // Add arguments if needed
    arg.resize(num_in);

    // Replace nulls with zeros of the right dimension
    for (int i=0; i<arg.size(); ++i) {
      if (arg[i].isEmpty(true)) arg[i] = MX::zeros(fcn_.input(i).sparsity());
    }

    setDependencies(arg);
    setSparsity(Sparsity::scalar());
  }

  CallFunction* CallFunction::clone() const {
    return new CallFunction(*this);
  }

  void CallFunction::printPart(std::ostream &stream, int part) const {
    fcn_->printPart(this, stream, part);
  }

  void CallFunction::evaluateD(const DMatrixPtrV& arg, DMatrixPtrV& res, std::vector<int>& itmp,
                               std::vector<double>& rtmp) {
    fcn_->evaluateD(this, arg, res, itmp, rtmp);
  }

  int CallFunction::getNumOutputs() const {
    return fcn_.getNumOutputs();
  }

  const Sparsity& CallFunction::sparsity(int oind) const {
    return fcn_.output(oind).sparsity();
  }

  Function& CallFunction::getFunction() {
    return fcn_;
  }

  void CallFunction::evaluateSX(const SXPtrV& arg, SXPtrV& res, std::vector<int>& itmp,
                                std::vector<SXElement>& rtmp) {
    fcn_->evaluateSX(this, arg, res, itmp, rtmp);
  }

  void CallFunction::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                                MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                                bool output_given) {
    fcn_->evaluateMX(this, input, output, fwdSeed, fwdSens, adjSeed, adjSens, output_given);
  }

  void CallFunction::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    MXNode::deepCopyMembers(already_copied);
    fcn_ = deepcopy(fcn_, already_copied);
  }

  void CallFunction::propagateSparsity(DMatrixPtrV& arg, DMatrixPtrV& res, std::vector<int>& itmp,
                                       std::vector<double>& rtmp, bool use_fwd) {
    fcn_->propagateSparsity(this, arg, res, itmp, rtmp, use_fwd);
  }

  void CallFunction::generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                       const std::vector<std::string>& res,
                                       CodeGenerator& gen) const {
    fcn_->generateOperation(this, stream, arg, res, gen);
  }

  void CallFunction::nTmp(size_t& ni, size_t& nr) {
    fcn_->nTmp(this, ni, nr);
  }

} // namespace casadi
