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


#include "casadi_call.hpp"
#include "../function/function_internal.hpp"
#include "../std_vector_tools.hpp"

using namespace std;

namespace casadi {

  MX GenericCall::projectArg(const MX& x, const Sparsity& sp, int i) {
    if (x.size()==sp.size()) {
      // Insert sparsity projection nodes if needed
      return project(x, sp);
    } else {
      // Different dimensions
      if (x.isempty() || sp.isempty()) { // NOTE: To permissive?
        // Replace nulls with zeros of the right dimension
        return MX::zeros(sp);
      } else if (x.isscalar()) {
        // Scalar argument means set all
        return MX(sp, x);
      } else if (x.size1()==sp.size2() && x.size2()==sp.size1() && sp.isvector()) {
        // Transposed vector
        return projectArg(x.T(), sp, i);
      } else {
        // Mismatching dimensions
        casadi_error("Cannot create function call node: Dimension mismatch for argument "
                     << i << ". Argument has shape " << x.size()
                     << " but function input has shape " << sp.size());
      }
    }
  }

  Call::Call(const Function& fcn, const vector<MX>& arg) : fcn_(fcn) {

    // Number inputs and outputs
    int num_in = fcn.nIn();
    casadi_assert_message(arg.size()==num_in, "Argument list length (" << arg.size()
                          << ") does not match number of inputs (" << num_in << ")"
                          << " for function " << fcn.name());

    // Create arguments of the right dimensions and sparsity
    vector<MX> arg1(num_in);
    for (int i=0; i<num_in; ++i) {
      arg1[i] = projectArg(arg[i], fcn_.input(i).sparsity(), i);
    }
    setDependencies(arg1);
    setSparsity(Sparsity::scalar());
  }

  std::string Call::print(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << fcn_.name() << "(";
    for (int i=0; i<ndep(); ++i) {
      if (i!=0) ss << ", ";
      ss << arg.at(i);
    }
    ss << ")";
    return ss.str();
  }

  void Call::evalD(const double** arg, double** res,
                           int* iw, double* w) {
    fcn_->eval(arg, res, iw, w);
  }

  int Call::nout() const {
    return fcn_.nOut();
  }

  const Sparsity& Call::sparsity(int oind) const {
    return fcn_.output(oind).sparsity();
  }

  void Call::evalSX(const SXElement** arg, SXElement** res, int* iw, SXElement* w) {
    fcn_->evalSX(arg, res, iw, w);
  }

  void Call::evalMX(const vector<MX>& arg, vector<MX>& res) {
    res = create(fcn_, arg);
  }

  void Call::evalFwd(const vector<vector<MX> >& fseed,
                     vector<vector<MX> >& fsens) {
    // Nondifferentiated inputs and outputs
    vector<MX> arg(ndep());
    for (int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(nout());
    for (int i=0; i<res.size(); ++i) res[i] = getOutput(i);

    // Call the cached functions
    fcn_.callForward(arg, res, fseed, fsens);
  }

  void Call::evalAdj(const vector<vector<MX> >& aseed,
                     vector<vector<MX> >& asens) {
    // Nondifferentiated inputs and outputs
    vector<MX> arg(ndep());
    for (int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(nout());
    for (int i=0; i<res.size(); ++i) res[i] = getOutput(i);

    // Call the cached functions
    vector<vector<MX> > v;
    fcn_.callReverse(arg, res, aseed, v);
    for (int i=0; i<v.size(); ++i) {
      for (int j=0; j<v[i].size(); ++j) {
        if (!v[i][j].isempty()) { // TODO(@jaeandersson): Hack
          asens[i][j] += v[i][j];
        }
      }
    }
  }

  void Call::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    fcn_.spFwd(arg, res, iw, w);
  }

  void Call::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    fcn_.spAdj(arg, res, iw, w);
  }

  void Call::addDependency(CodeGenerator& g) const {
    fcn_->addDependency(g);
  }

  void Call::generate(const vector<int>& arg, const vector<int>& res,
                      CodeGenerator& g) const {
    fcn_->generate(g, arg, res);
  }

  size_t Call::sz_arg() const {
    return fcn_.sz_arg();
  }

  size_t Call::sz_res() const {
    return fcn_.sz_res();
  }

  size_t Call::sz_iw() const {
    return fcn_.sz_iw();
  }

  size_t Call::sz_w() const {
    return fcn_.sz_w();
  }

  std::vector<MX> Call::create(const Function& fcn, const std::vector<MX>& arg) {
    return MX::createMultipleOutput(new Call(fcn, arg));
  }

} // namespace casadi
