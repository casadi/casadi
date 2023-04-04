/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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

  max_num_threads_ = 1;

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

  stride_arg_ = 0;
  stride_res_ = 0;
  stride_iw_ = 0;
  stride_w_ = 0;

}

void OracleFunction::finalize() {

  // Allocate space for (parallel) evaluations
  // Lifted from set_function as max_num_threads_ is not known yet in that method
  for (auto&& e : all_functions_) {
    Function& fcn = e.second.f;
    // Compute strides for multi threading
    size_t sz_arg, sz_res, sz_iw, sz_w;
    fcn.sz_work(sz_arg, sz_res, sz_iw, sz_w);
    stride_arg_ = std::max(stride_arg_, sz_arg);
    stride_res_ = std::max(stride_res_, sz_res);
    stride_iw_ = std::max(stride_iw_, sz_iw);
    stride_w_ = std::max(stride_w_, sz_w);
    bool persistent = false;
    alloc(fcn, persistent, max_num_threads_);
  }

  // Set corresponding monitors
  for (const std::string& fname : monitor_) {
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

void OracleFunction::join_results(OracleMemory* m) const {
  // Combine runtime statistics
  // Note: probably not correct to simply add wall times
  for (int i = 0; i < max_num_threads_; ++i) {
    auto* ml = m->thread_local_mem[i];
    for (auto&& s : ml->fstats) {
      m->fstats.at(s.first).join(s.second);
    }
  }
}

Function OracleFunction::create_function(const std::string& fname,
    const std::vector<std::string>& s_in,
    const std::vector<std::string>& s_out,
    const Function::AuxOut& aux) {
  return create_function(oracle_, fname, s_in, s_out, aux);
}

Function OracleFunction::create_function(const Function& oracle, const std::string& fname,
    const std::vector<std::string>& s_in,
    const std::vector<std::string>& s_out,
    const Function::AuxOut& aux) {
  // Print progress
  if (verbose_) {
    casadi_message(name_ + "::create_function " + fname + ":" + str(s_in) + "->" + str(s_out));
  }

  // Check if function is already in cache
  Function ret;
  if (incache(fname, ret)) {
    // Consistency checks
    casadi_assert(ret.n_in() == s_in.size(), fname + " has wrong number of inputs");
    casadi_assert(ret.n_out() == s_out.size(), fname + " has wrong number of outputs");
  } else {
    // Retrieve specific set of options if available
    Dict specific_options;
    auto it = specific_options_.find(fname);
    if (it!=specific_options_.end()) specific_options = it->second;

    // Combine specific and common options
    Dict opt = combine(specific_options, common_options_);

    // Generate the function
    ret = oracle.factory(fname, s_in, s_out, aux, opt);

    // Make sure that it's sound
    if (ret.has_free()) {
      casadi_error("Cannot create '" + fname + "' since " + str(ret.get_free()) + " are free.");
    }

    // Add to cache
    tocache(ret);
  }

  // Save and return
  set_function(ret, fname, true);
  return ret;
}

Function OracleFunction::create_forward(const std::string& fname, casadi_int nfwd) {
  // Create derivative
  Function ret = get_function(fname).forward(nfwd);
  if (!has_function(ret.name())) set_function(ret, ret.name(), true);
  return ret;
}

void OracleFunction::
set_function(const Function& fcn, const std::string& fname, bool jit) {
  casadi_assert(!has_function(fname), "Duplicate function " + fname);
  RegFun& r = all_functions_[fname];
  r.f = fcn;
  r.jit = jit;
}

int OracleFunction::
calc_function(OracleMemory* m, const std::string& fcn,
              const double* const* arg, int thread_id) const {
  auto ml = m->thread_local_mem.at(thread_id);
  // Is the function monitored?
  bool monitored = this->monitored(fcn);

  // Print progress
  if (monitored) casadi_message("Calling \"" + fcn + "\"");

  // Respond to a possible Crl+C signals
  // Python interrupt checker needs the GIL.
  // We may not have access to it in a multi-threaded context
  // See issue #2955
  if (max_num_threads_==1) InterruptHandler::check();

  // Get function
  const Function& f = get_function(fcn);

  // Get statistics structure
  FStats& fstats = ml->fstats.at(fcn);

  // Number of inputs and outputs
  casadi_int n_in = f.n_in(), n_out = f.n_out();

  // Prepare stats, start timer
  ScopedTiming tic(fstats);

  // Input buffers
  if (arg) {
    std::fill_n(ml->arg, n_in, nullptr);
    for (casadi_int i=0; i<n_in; ++i) ml->arg[i] = *arg++;
  }

  // Print inputs nonzeros
  if (monitored) {
    std::stringstream s;
    s << fcn << " input nonzeros:\n";
    for (casadi_int i=0; i<n_in; ++i) {
      s << " " << i << " (" << f.name_in(i) << "): ";
      if (ml->arg[i]) {
        // Print nonzeros
        s << "[";
        for (casadi_int k=0; k<f.nnz_in(i); ++k) {
          if (k!=0) s << ", ";
          DM::print_scalar(s, ml->arg[i][k]);
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
    if (f(ml->arg, ml->res, ml->iw, ml->w)) {
      // Recoverable error
      if (monitored) casadi_message(name_ + ":" + fcn + " failed");
      return 1;
    }
  } catch(std::exception& ex) {
    // Fatal error: Generate stack trace
    casadi_error("Error in " + name_ + ":" + fcn + ":" + std::string(ex.what()));
  }

  // Print output nonzeros
  if (monitored) {
    std::stringstream s;
    s << fcn << " output nonzeros:\n";
    for (casadi_int i=0; i<n_out; ++i) {
      s << " " << i << " (" << f.name_out(i) << "): ";
      if (ml->res[i]) {
        // Print nonzeros
        s << "[";
        for (casadi_int k=0; k<f.nnz_out(i); ++k) {
          if (k!=0) s << ", ";
          DM::print_scalar(s, ml->res[i][k]);
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
    if (!ml->res[i]) continue;
    if (!std::all_of(ml->res[i], ml->res[i]+f.nnz_out(i), [](double v) { return isfinite(v);})) {
      std::stringstream ss;

      auto it = std::find_if(ml->res[i], ml->res[i] + f.nnz_out(i),
        [](double v) { return !isfinite(v);});
      casadi_int k = std::distance(ml->res[i], it);
      bool is_nan = isnan(ml->res[i][k]);
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

int OracleFunction::calc_sp_forward(const std::string& fcn, const bvec_t** arg, bvec_t** res,
    casadi_int* iw, bvec_t* w) const {
  return get_function(fcn)(arg, res, iw, w);
}

int OracleFunction::calc_sp_reverse(const std::string& fcn, bvec_t** arg, bvec_t** res,
    casadi_int* iw, bvec_t* w) const {
  return get_function(fcn).rev(arg, res, iw, w);
}

std::string OracleFunction::generate_dependencies(const std::string& fname,
    const Dict& opts) const {
  CodeGenerator gen(fname, opts);
  gen.add(oracle_);
  for (auto&& e : all_functions_) {
    if (e.second.jit) gen.add(e.second.f);
  }
  return gen.generate();
}

void OracleFunction::jit_dependencies(const std::string& fname) {
  if (compiler_.is_null()) {
    if (verbose_) casadi_message("compiling to "+ fname+"'.");
    // JIT dependent functions
    compiler_ = Importer(generate_dependencies(fname, Dict()),
                        compiler_plugin_, jit_options_);
  }
  // Replace the Oracle functions with generated functions
  for (auto&& e : all_functions_) {
    if (verbose_) casadi_message("loading '" + e.second.f.name() + "' from '" + fname + "'.");
    if (e.second.jit) {
      e.second.f_original = e.second.f;
      e.second.f = external(e.second.f.name(), compiler_);
    }
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

int OracleFunction::local_init_mem(void* mem) const {
  if (ProtoFunction::init_mem(mem)) return 1;
  if (!mem) return 1;
  auto m = static_cast<LocalOracleMemory*>(mem);

  // Create statistics
  for (auto&& e : all_functions_) {
    m->add_stat(e.first);
  }

  return 0;
}

int OracleFunction::init_mem(void* mem) const {
  if (ProtoFunction::init_mem(mem)) return 1;
  if (!mem) return 1;
  auto m = static_cast<OracleMemory*>(mem);

  // Create statistics
  for (auto&& e : all_functions_) {
    m->add_stat(e.first);
  }

  casadi_assert_dev(m->thread_local_mem.empty());

  // Allocate and initialize local memory for threads
  for (int i = 0; i < max_num_threads_; ++i) {
    m->thread_local_mem.push_back(new LocalOracleMemory());
    if (OracleFunction::local_init_mem(m->thread_local_mem[i])) return 1;
  }

  return 0;
}

void OracleFunction::free_mem(void *mem) const {
  auto m = static_cast<OracleMemory*>(mem);

  for (int i = 0; i < max_num_threads_; ++i) {
    local_free_mem(m->thread_local_mem[i]);
  }

  delete m;
}

void OracleFunction::set_temp(void* mem, const double** arg, double** res,
                          casadi_int* iw, double* w) const {

  auto m = static_cast<OracleMemory*>(mem);
  m->arg = arg;
  m->res = res;
  m->iw = iw;
  m->w = w;
  for (int i = 0; i < max_num_threads_; ++i) {
    auto* ml = m->thread_local_mem[i];
    ml->arg = arg;
    ml->res = res;
    ml->iw = iw;
    ml->w = w;
    arg += stride_arg_;
    res += stride_res_;
    iw += stride_iw_;
    w += stride_w_;
  }
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

  s.version("OracleFunction", 3);
  s.pack("OracleFunction::oracle", oracle_);
  s.pack("OracleFunction::common_options", common_options_);
  s.pack("OracleFunction::specific_options", specific_options_);
  s.pack("OracleFunction::show_eval_warnings", show_eval_warnings_);
  s.pack("OracleFunction::max_num_threads", max_num_threads_);
  s.pack("OracleFunction::all_functions::size", all_functions_.size());
  for (auto &e : all_functions_) {
    s.pack("OracleFunction::all_functions::key", e.first);
    s.pack("OracleFunction::all_functions::value::jit", e.second.jit);
    if (jit_ && e.second.jit) {
      if (jit_serialize_=="source") {
        // Save original f, such that it can be built
        s.pack("OracleFunction::all_functions::value::f", e.second.f_original);
      } else {
        std::string f_name = e.second.f.name();
        s.pack("OracleFunction::all_functions::value::f_name", f_name);
        // FunctionInternal will set compiler_
      }
    } else {
      // Save f
      s.pack("OracleFunction::all_functions::value::f", e.second.f);
    }
    s.pack("OracleFunction::all_functions::value::monitored", e.second.monitored);
  }
  s.pack("OracleFunction::monitor", monitor_);
  s.pack("OracleFunction::stride_arg", stride_arg_);
  s.pack("OracleFunction::stride_res", stride_res_);
  s.pack("OracleFunction::stride_iw", stride_iw_);
  s.pack("OracleFunction::stride_w", stride_w_);

}

OracleFunction::OracleFunction(DeserializingStream& s) : FunctionInternal(s) {

  int version = s.version("OracleFunction", 1, 3);
  s.unpack("OracleFunction::oracle", oracle_);
  s.unpack("OracleFunction::common_options", common_options_);
  s.unpack("OracleFunction::specific_options", specific_options_);
  s.unpack("OracleFunction::show_eval_warnings", show_eval_warnings_);

  if (version>=3) {
    s.unpack("OracleFunction::max_num_threads", max_num_threads_);
  } else {
    max_num_threads_ = 1;
  }

  size_t size;

  s.unpack("OracleFunction::all_functions::size", size);
  for (casadi_int i=0;i<size;++i) {
    std::string key;
    s.unpack("OracleFunction::all_functions::key", key);
    RegFun r;
    if (version==1) {
      s.unpack("OracleFunction::all_functions::value::f", r.f);
      s.unpack("OracleFunction::all_functions::value::jit", r.jit);
    } else {
      s.unpack("OracleFunction::all_functions::value::jit", r.jit);
      if (jit_ && r.jit) {
        if (jit_serialize_=="source") {
          s.unpack("OracleFunction::all_functions::value::f", r.f);
        } else {
          std::string f_name;
          s.unpack("OracleFunction::all_functions::value::f_name", f_name);
          r.f = Function(f_name, std::vector<MX>{}, std::vector<MX>{});
          // FunctionInternal will set compiler_
        }
      } else {
        s.unpack("OracleFunction::all_functions::value::f", r.f);
      }
    }
    s.unpack("OracleFunction::all_functions::value::monitored", r.monitored);
    all_functions_[key] = r;
  }
  s.unpack("OracleFunction::monitor", monitor_);
  if (version>=3) {
    s.unpack("OracleFunction::stride_arg", stride_arg_);
    s.unpack("OracleFunction::stride_res", stride_res_);
    s.unpack("OracleFunction::stride_iw", stride_iw_);
    s.unpack("OracleFunction::stride_w", stride_w_);
  } else {
    stride_arg_ = 0;
    stride_res_ = 0;
    stride_iw_ = 0;
    stride_w_ = 0;
  }
}

} // namespace casadi
