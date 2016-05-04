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


#include "oracle_function.hpp"

using namespace std;

namespace casadi {

  OracleFunction::OracleFunction(const std::string& name, const Function& oracle)
  : FunctionInternal(name), oracle_(oracle) {
  }

  OracleFunction::~OracleFunction() {
  }

  Function OracleFunction::create_function(const std::string& fname,
                                   const std::vector<std::string>& s_in,
                                   const std::vector<std::string>& s_out,
                                   const Function::AuxOut& aux,
                                   const Dict& opts, bool reg) {
    // Generate the function
    Function ret = oracle_.factory(fname, s_in, s_out, aux, opts);
    if (reg) register_function(ret);
    return ret;
  }

  void OracleFunction::register_function(const Function& fcn) {
    // Make sure that the function names are unique
    for (const Function& f : all_functions_) casadi_assert(fcn.name()!=f.name());

    // Add to the list
    all_functions_.push_back(fcn);
    alloc(fcn);
  }

  int OracleFunction::calc_function(OracleMemory* m, const Function& fcn,
                            std::initializer_list<const double*> arg,
                            std::initializer_list<double*> res) const {
    // Respond to a possible Crl+C signals
    InterruptHandler::check();

    // Get statistics structure
    FStats& fstats = m->fstats.at(fcn.name());

    // Prepare stats, start timer
    fstats.tic();

    // Number of inputs and outputs
    int n_in = fcn.n_in(), n_out = fcn.n_out();

    // Input buffers
    fill_n(m->arg, n_in, nullptr);
    auto arg_it = arg.begin();
    for (int i=0; i<n_in; ++i) m->arg[i] = *arg_it++;
    casadi_assert(arg_it==arg.end());

    // Output buffers
    fill_n(m->res, n_out, nullptr);
    auto res_it = res.begin();
    for (int i=0; i<n_out; ++i) m->res[i] = *res_it++;
    casadi_assert(res_it==res.end());

    // Evaluate memory-less
    try {
      fcn(m->arg, m->res, m->iw, m->w, 0);
    } catch(exception& ex) {
      // Fatal error
      userOut<true, PL_WARN>()
        << name() << ":" << fcn.name() << " failed:" << ex.what() << endl;
      return 1;
    }

    // Make sure not NaN or Inf
    for (int i=0; i<n_out; ++i) {
      if (!m->res[i]) continue;
      if (!all_of(m->res[i], m->res[i]+fcn.nnz_out(i), [](double v) { return isfinite(v);})) {
        userOut<true, PL_WARN>()
          << name() << ":" << fcn.name() << " failed: NaN or Inf detected for output "
          << fcn.name_out(i) << endl;
        return -1;
      }
    }

    // Update stats
    fstats.toc();

    // Success
    return 0;
  }

  void OracleFunction::generate_dependencies(const std::string& fname, const Dict& opts) {
    CodeGenerator gen(fname, opts);
    gen.add(oracle_);
    for (const Function& f : all_functions_) gen.add(f);
    gen.generate();
  }

  void OracleFunction::expand() {
    oracle_ = oracle_.expand();
  }

} // namespace casadi
