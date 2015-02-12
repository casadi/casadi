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

    // Assert argument of correct dimension and sparsity
    for (int i=0; i<arg.size(); ++i) {
      if (arg[i].shape()==fcn_.input(i).shape()) {
        // Insert sparsity projection nodes if needed
        arg[i] = arg[i].setSparse(fcn_.input(i).sparsity());
      } else {
        // Different dimensions
        if (arg[i].isEmpty() || fcn_.input(i).isEmpty()) { // NOTE: To permissive?
          // Replace nulls with zeros of the right dimension
          arg[i] = MX::zeros(fcn_.input(i).sparsity());
        } else if (arg[i].isScalar()) {
          // Scalar argument means set all
          arg[i] = MX(fcn_.input(i).sparsity(), arg[i]);
        } else {
          // Mismatching dimensions
          casadi_error("Cannot create function call node: Dimension mismatch for argument "
                       << i << ". Argument has shape " << arg[i].shape()
                       << " but function input is " << fcn_.input(i).shape());
        }
      }
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

  void CallFunction::evaluateD(const double* const* input, double** output,
                               int* itmp, double* rtmp) {
    // Set up timers for profiling
    double time_zero=0;
    double time_start=0;
    double time_stop=0;
    double time_offset=0;

    // Number of inputs and outputs
    int num_in = fcn_.getNumInputs();
    int num_out = fcn_.getNumOutputs();

    // Pass the inputs to the function
    for (int i = 0; i < num_in; ++i) {
      if (input[i] != 0) {
        fcn_.setInput(input[i], i);
      } else {
        fcn_.setInput(0., i);
      }
    }

    // Evaluate
    fcn_.evaluate();

    // Get the outputs
    for (int i=0; i<num_out; ++i) {
      if (output[i] != 0) fcn_.getOutput(output[i], i);
    }
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

  void CallFunction::evaluateSX(const SXElement* const* input, SXElement** output,
                                int* itmp, SXElement* rtmp) {
    // Number of inputs and outputs
    int num_in = fcn_.getNumInputs();
    int num_out = fcn_.getNumOutputs();

    // Create input arguments
    vector<SX> argv(num_in);
    for (int i=0; i<num_in; ++i) {
      argv[i] = SX::zeros(fcn_.input(i).sparsity());
      if (input[i] != 0) argv[i].set(input[i]);
    }

    // Evaluate symbolically
    vector<SX> resv;
    fcn_->evalSX(argv, resv);

    // Collect the result
    for (int i = 0; i < num_out; ++i) {
      if (output[i] != 0) resv[i].get(output[i]);
    }
  }

  void CallFunction::eval(const MXPtrV& input, MXPtrV& output) {
    vector<MX> arg = getVector(input);
    vector<MX> res = fcn_->createCall(arg);
    for (int i=0; i<res.size(); ++i) {
      if (output[i]!=0) {
        *output[i] = res[i];
      }
    }
  }

  void CallFunction::evalFwd(const MXPtrVV& fwdSeed, MXPtrVV& fwdSens) {
    // Nondifferentiated inputs and outputs
    vector<MX> arg(ndep());
    for (int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(getNumOutputs());
    for (int i=0; i<res.size(); ++i) res[i] = getOutput(i);

    // Collect seeds
    vector<vector<MX> > fseed(getVector(fwdSeed)), fsens;

    // Call the cached functions
    fcn_.callFwd(arg, res, fseed, fsens);

    // Store the forward sensitivities
    for (int d=0; d<fwdSens.size(); ++d) {
      for (int i=0; i<fwdSens[d].size(); ++i) {
        if (fwdSens[d][i]!=0) {
          *fwdSens[d][i] = fsens[d][i];
        }
      }
    }
  }

  void CallFunction::evalAdj(MXPtrVV& adjSeed, MXPtrVV& adjSens) {
    // Nondifferentiated inputs and outputs
    vector<MX> arg(ndep());
    for (int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(getNumOutputs());
    for (int i=0; i<res.size(); ++i) res[i] = getOutput(i);

    // Collect seeds
    vector<vector<MX> > aseed(getVector(adjSeed)), asens;

    // Call the cached functions
    fcn_.callAdj(arg, res, aseed, asens);

    // Free adjoint seeds
    clearVector(adjSeed);

    // Store the adjoint sensitivities
    for (int d=0; d<adjSens.size(); ++d) {
      for (int i=0; i<adjSens[d].size(); ++i) {
        if (adjSens[d][i]!=0) {
          adjSens[d][i]->addToSum(asens[d][i]);
        }
      }
    }
  }

  void CallFunction::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    MXNode::deepCopyMembers(already_copied);
    fcn_ = deepcopy(fcn_, already_copied);
  }

  void CallFunction::propagateSparsity(double** input, double** output,
                                       int* itmp, bvec_t* rtmp, bool use_fwd) {
    fcn_->propagateSparsity(this, input, output, itmp, rtmp, use_fwd);
  }

  void CallFunction::generateOperation(std::ostream &stream, const std::vector<int>& arg,
                                       const std::vector<int>& res,
                                       CodeGenerator& gen) const {
    fcn_->generateOperation(this, stream, arg, res, gen);
  }

  void CallFunction::nTmp(size_t& ni, size_t& nr) {
    fcn_->nTmp(this, ni, nr);
  }

} // namespace casadi
