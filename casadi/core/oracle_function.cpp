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
#include "external.hpp"
#include "serializing_stream.hpp"

#include <iomanip>
#include <iostream>

using namespace std;

namespace casadi {

  OracleFunction::OracleFunction(const std::string& name, const Function& oracle)
  : FunctionInternal(name), oracle_(oracle) {
  }

  OracleFunction::~OracleFunction() {
  }

  const Options OracleFunction::options_
  = {{&FunctionInternal::options_},
     {{"expand",
       {OT_BOOL,
        "Replace MX with SX expressions in problem formulation [false]"}},
      {"monitor",
       {OT_STRINGVECTOR,
        "Set of user problem functions to be monitored"}},
      {"show_eval_warnings",
       {OT_BOOL,
        "Show warnings generated from function evaluations [true]"}},
      {"common_options",
       {OT_DICT,
        "Options for auto-generated functions"}},
      {"specific_options",
       {OT_DICT,
        "Options for specific auto-generated functions,"
        " overwriting the defaults from common_options. Nested dictionary."}}
    }
  };

  void OracleFunction::init(const Dict& opts) {

    FunctionInternal::init(opts);

    // Default options
    bool expand = false;

    show_eval_warnings_ = true;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="expand") {
        expand = op.second;
      } else if (op.first=="common_options") {
        common_options_ = op.second;
      } else if (op.first=="specific_options") {
        specific_options_ = op.second;
        for (auto&& i : specific_options_) {
          casadi_assert(i.second.is_dict(),
            "specific_option must be a nested dictionary."
            " Type mismatch for entry '" + i.first+ "': "
            " got type " + i.second.get_description() + ".");
        }
      } else if (op.first=="monitor") {
        monitor_ = op.second;
      } else if (op.first=="show_eval_warnings") {
        show_eval_warnings_ = op.second;
      }
    }

    // Replace MX oracle with SX oracle?
    if (expand) oracle_ = oracle_.expand();

  }

  void OracleFunction::finalize() {


    // Set corresponding monitors
    for (const string& fname : monitor_) {
      auto it = all_functions_.find(fname);
      if (it==all_functions_.end()) {
        casadi_warning("Ignoring monitor '" + fname + "'."
                       " Available functions: " + join(get_function()) + ".");
      } else {
        if (it->second.monitored) casadi_warning("Duplicate monitor " + fname);
        it->second.monitored = true;
      }
    }

    // Check specific options
    for (auto&& i : specific_options_) {
      if (all_functions_.find(i.first)==all_functions_.end())
        casadi_warning("Ignoring specific_options entry '" + i.first+"'."
                       " Available functions: " + join(get_function()) + ".");
    }

    // Recursive call
    FunctionInternal::finalize();
  }

  Function OracleFunction::create_function(const std::string& fname,
                                   const std::vector<std::string>& s_in,
                                   const std::vector<std::string>& s_out,
                                   const Function::AuxOut& aux) {
    // Print progress
    if (verbose_) {
      casadi_message(name_ + "::create_function " + fname + ":" + str(s_in) + "->" + str(s_out));
    }

    // Retrieve specific set of options if available
    Dict specific_options;
    auto it = specific_options_.find(fname);
    if (it!=specific_options_.end()) specific_options = it->second;

    // Combine specific and common options
    Dict opt = combine(specific_options, common_options_);

    // Generate the function
    Function ret = oracle_.factory(fname, s_in, s_out, aux, opt);

    // Make sure that it's sound
    if (ret.has_free()) {
      casadi_error("Cannot create '" + fname + "' since " + str(ret.get_free()) + " are free.");
    }

    // Save and return
    set_function(ret, fname, true);
    return ret;
  }

  void OracleFunction::
  set_function(const Function& fcn, const std::string& fname, bool jit) {
    casadi_assert(!has_function(fname), "Duplicate function " + fname);
    RegFun& r = all_functions_[fname];
    r.f = fcn;
    r.jit = jit;
    alloc(fcn);
  }

  int OracleFunction::
  calc_function(OracleMemory* m, const std::string& fcn,
                const double* const* arg) const {
    // Is the function monitored?
    bool monitored = this->monitored(fcn);

    // Print progress
    if (monitored) casadi_message("Calling \"" + fcn + "\"");

    // Respond to a possible Crl+C signals
    InterruptHandler::check();

    // Get function
    const Function& f = get_function(fcn);

    // Get statistics structure
    FStats& fstats = m->fstats.at(fcn);

    // Number of inputs and outputs
    casadi_int n_in = f.n_in(), n_out = f.n_out();

    // Prepare stats, start timer
    ScopedTiming tic(fstats);

    // Input buffers
    if (arg) {
      fill_n(m->arg, n_in, nullptr);
      for (casadi_int i=0; i<n_in; ++i) m->arg[i] = *arg++;
    }

    // Print inputs nonzeros
    if (monitored) {
      std::stringstream s;
      s << fcn << " input nonzeros:\n";
      for (casadi_int i=0; i<n_in; ++i) {
        s << " " << i << " (" << f.name_in(i) << "): ";
        if (m->arg[i]) {
          // Print nonzeros
          s << "[";
          for (casadi_int k=0; k<f.nnz_in(i); ++k) {
            if (k!=0) s << ", ";
            DM::print_scalar(s, m->arg[i][k]);
          }
          s << "]\n";
        } else {
          // All zero input
          s << "0\n";
        }
      }
      casadi_message(s.str());
    }

    // Evaluate memory-less
    try {
      f(m->arg, m->res, m->iw, m->w);
    } catch(exception& ex) {
      // Fatal error
      if (show_eval_warnings_) {
        casadi_warning(name_ + ":" + fcn + " failed:" + std::string(ex.what()));
      }
      return 1;
    }

    // Print output nonzeros
    if (monitored) {
      std::stringstream s;
      s << fcn << " output nonzeros:\n";
      for (casadi_int i=0; i<n_out; ++i) {
        s << " " << i << " (" << f.name_out(i) << "): ";
        if (m->res[i]) {
          // Print nonzeros
          s << "[";
          for (casadi_int k=0; k<f.nnz_out(i); ++k) {
            if (k!=0) s << ", ";
            DM::print_scalar(s, m->res[i][k]);
          }
          s << "]\n";
        } else {
          // Ignored output
          s << " N/A\n";
        }
      }
      casadi_message(s.str());
    }

    // Make sure not NaN or Inf
    for (casadi_int i=0; i<n_out; ++i) {
      if (!m->res[i]) continue;
      if (!all_of(m->res[i], m->res[i]+f.nnz_out(i), [](double v) { return isfinite(v);})) {
        std::stringstream ss;

        auto it = find_if(m->res[i], m->res[i]+f.nnz_out(i), [](double v) { return !isfinite(v);});
        casadi_int k = distance(m->res[i], it);
        bool is_nan = isnan(m->res[i][k]);
        ss << name_ << ":" << fcn << " failed: " << (is_nan? "NaN" : "Inf") <<
        " detected for output " << f.name_out(i) << ", at " << f.sparsity_out(i).repr_el(k) << ".";

        if (regularity_check_) {
          casadi_error(ss.str());
        } else {
          if (show_eval_warnings_) casadi_warning(ss.str());
        }
        return -1;
      }
    }

    // Success
    return 0;
  }

  std::string OracleFunction::
  generate_dependencies(const std::string& fname, const Dict& opts) const {
    CodeGenerator gen(fname, opts);
    gen.add(oracle_);
    for (auto&& e : all_functions_) {
      if (e.second.jit) gen.add(e.second.f);
    }
    return gen.generate();
  }

  void OracleFunction::jit_dependencies(const std::string& fname) {
    if (verbose_)
      if (verbose_) casadi_message("compiling to "+ fname+"'.");
    // JIT dependent functions
    compiler_ = Importer(generate_dependencies(fname, Dict()),
                         compiler_plugin_, jit_options_);

    // Replace the Oracle functions with generated functions
    for (auto&& e : all_functions_) {
      if (verbose_)
        if (verbose_) casadi_message("loading '" + e.second.f.name() + "' from '" + fname + "'.");
      if (e.second.jit) e.second.f = external(e.second.f.name(), compiler_);
    }
  }

  void OracleFunction::expand() {
    oracle_ = oracle_.expand();
  }

  Dict OracleFunction::get_stats(void *mem) const {
    Dict stats = FunctionInternal::get_stats(mem);
    //auto m = static_cast<OracleMemory*>(mem);
    return stats;
  }

  int OracleFunction::init_mem(void* mem) const {
    if (ProtoFunction::init_mem(mem)) return 1;
    if (!mem) return 1;
    auto m = static_cast<OracleMemory*>(mem);

    // Create statistics
    for (auto&& e : all_functions_) {
      m->add_stat(e.first);
    }
    return 0;
  }

  void OracleFunction::set_temp(void* mem, const double** arg, double** res,
                            casadi_int* iw, double* w) const {
    auto m = static_cast<OracleMemory*>(mem);
    m->arg = arg;
    m->res = res;
    m->iw = iw;
    m->w = w;
  }

  std::vector<std::string> OracleFunction::get_function() const {
    std::vector<std::string> ret;
    ret.reserve(all_functions_.size());
    for (auto&& e : all_functions_) {
      ret.push_back(e.first);
    }
    return ret;
  }

  const Function& OracleFunction::get_function(const std::string &name) const {
    auto it = all_functions_.find(name);
    casadi_assert(it!=all_functions_.end(),
      "No function \"" + name + "\" in " + name_ + ". " +
      "Available functions: " + join(get_function()) + ".");
    return it->second.f;
  }

  bool OracleFunction::monitored(const std::string &name) const {
    auto it = all_functions_.find(name);
    casadi_assert(it!=all_functions_.end(),
      "No function \"" + name + "\" in " + name_+ ". " +
      "Available functions: " + join(get_function()) + ".");
    return it->second.monitored;
  }

  bool OracleFunction::has_function(const std::string& fname) const {
    return all_functions_.find(fname) != all_functions_.end();
  }


  void OracleFunction::serialize_body(SerializingStream &s) const {
    FunctionInternal::serialize_body(s);

    s.version("OracleFunction", 1);
    s.pack("OracleFunction::oracle", oracle_);
    s.pack("OracleFunction::common_options", common_options_);
    s.pack("OracleFunction::specific_options", specific_options_);
    s.pack("OracleFunction::show_eval_warnings", show_eval_warnings_);
    s.pack("OracleFunction::all_functions::size", all_functions_.size());
    for (auto &e : all_functions_) {
      s.pack("OracleFunction::all_functions::key", e.first);
      s.pack("OracleFunction::all_functions::value::f", e.second.f);
      s.pack("OracleFunction::all_functions::value::jit", e.second.jit);
      s.pack("OracleFunction::all_functions::value::monitored", e.second.monitored);
    }
    s.pack("OracleFunction::monitor", monitor_);
  }

  OracleFunction::OracleFunction(DeserializingStream& s) : FunctionInternal(s) {

    s.version("OracleFunction", 1);
    s.unpack("OracleFunction::oracle", oracle_);
    s.unpack("OracleFunction::common_options", common_options_);
    s.unpack("OracleFunction::specific_options", specific_options_);
    s.unpack("OracleFunction::show_eval_warnings", show_eval_warnings_);
    size_t size;

    s.unpack("OracleFunction::all_functions::size", size);
    for (casadi_int i=0;i<size;++i) {
      std::string key;
      s.unpack("OracleFunction::all_functions::key", key);
      RegFun r;
      s.unpack("OracleFunction::all_functions::value::f", r.f);
      s.unpack("OracleFunction::all_functions::value::jit", r.jit);
      s.unpack("OracleFunction::all_functions::value::monitored", r.monitored);
      all_functions_[key] = r;
    }
    s.unpack("OracleFunction::monitor", monitor_);
  }

} // namespace casadi
