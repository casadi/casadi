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

#include <iostream>
#include <iomanip>

using namespace std;

namespace casadi {

  OracleFunction::OracleFunction(const std::string& name, const Function& oracle)
  : FunctionInternal(name), oracle_(oracle) {
  }

  OracleFunction::~OracleFunction() {
  }

  Options OracleFunction::options_
  = {{&FunctionInternal::options_},
     {{"monitor",
       {OT_STRINGVECTOR,
        "Set of user problem functions to be monitored"}},
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

    // Read options
    for (auto&& op : opts) {
      if (op.first=="common_options") {
        common_options_ = op.second;
      } else if (op.first=="specific_options") {
        specific_options_ = op.second;
        for (auto&& i : specific_options_) {
          casadi_assert(i.second.is_dict(),
            "specific_option must be a nested dictionary."
            " Type mismatch for entry '" + i.first+ "': "
            " got type " + i.second.get_description() + ".");
        }
      }
    }
  }

  void OracleFunction::finalize(const Dict& opts) {
    // Default options
    vector<string> monitor;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="monitor") {
        monitor = op.second;
      }
    }

    // Set corresponding monitors
    for (const string& fname : monitor) {
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
    FunctionInternal::finalize(opts);
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
    int n_in = f.n_in(), n_out = f.n_out();

    // Prepare stats, start timer
    fstats.tic();

    // Input buffers
    if (arg) {
      fill_n(m->arg, n_in, nullptr);
      for (int i=0; i<n_in; ++i) m->arg[i] = *arg++;
    }

    // Print inputs nonzeros
    if (monitored) {
      std::stringstream s;
      s << fcn << " input nonzeros:\n";
      for (int i=0; i<n_in; ++i) {
        s << " " << i << " (" << f.name_in(i) << "): ";
        if (m->arg[i]) {
          // Print nonzeros
          s << "[";
          for (int k=0; k<f.nnz_in(i); ++k) {
            if (k!=0) s << ", ";
            s << m->arg[i][k];
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
      f(m->arg, m->res, m->iw, m->w, 0);
    } catch(exception& ex) {
      // Fatal error
      casadi_warning(name() + ":" + fcn + " failed:" + std::string(ex.what()));
      return 1;
    }

    // Print output nonzeros
    if (monitored) {
      std::stringstream s;
      s << fcn << " output nonzeros:\n";
      for (int i=0; i<n_out; ++i) {
        s << " " << i << " (" << f.name_out(i) << "): ";
        if (m->res[i]) {
          // Print nonzeros
          s << "[";
          for (int k=0; k<f.nnz_out(i); ++k) {
            if (k!=0) s << ", ";
            s << m->res[i][k];
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
    for (int i=0; i<n_out; ++i) {
      if (!m->res[i]) continue;
      if (!all_of(m->res[i], m->res[i]+f.nnz_out(i), [](double v) { return isfinite(v);})) {
        std::stringstream ss;

        auto it = find_if(m->res[i], m->res[i]+f.nnz_out(i), [](double v) { return !isfinite(v);});
        int k = distance(m->res[i], it);
        bool is_nan = isnan(m->res[i][k]);
        ss << name() << ":" << fcn << " failed: " << (is_nan? "NaN" : "Inf") <<
        " detected for output " << f.name_out(i) << ", at " << f.sparsity_out(i).repr_el(k) << ".";

        if (regularity_check_) {
          casadi_error(ss.str());
        } else {
          casadi_warning(ss.str());
        }
        return -1;
      }
    }

    // Update stats
    fstats.toc();

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
                         compilerplugin_, jit_options_);

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

  void OracleFunction::print_fstats(const OracleMemory* m) const {

    size_t maxNameLen=0;

    // Retrieve all nlp keys
    std::vector<std::string> keys;
    std::vector<std::string> keys_other;
    for (auto &&s : m->fstats) {
      maxNameLen = max(s.first.size(), maxNameLen);
      if (s.first.find("nlp")!=std::string::npos) {
        keys.push_back(s.first);
      } else if (s.first.find("mainloop")==std::string::npos) {
        keys_other.push_back(s.first);
      } else {
        continue;
      }
    }

    maxNameLen = max(std::string("all previous").size(), maxNameLen);
    maxNameLen = max(std::string("solver").size(), maxNameLen);

    // Print header
    std::stringstream s;
    std::string blankName(maxNameLen, ' ');
    s << blankName
      << "      proc           wall      num           mean             mean"
      << endl << blankName
      << "      time           time     evals       proc time        wall time"
      << endl;
    uout() << s.str();

    // Sort the keys according to order
    std::vector<std::string> keys_order0;
    std::vector<std::string> keys_order1;
    std::vector<std::string> keys_order2;
    for (auto k : keys) {
      if (k.find("hess")!=std::string::npos) {
        keys_order2.push_back(k);
        continue;
      }
      if (k.find("grad")!=std::string::npos ||
          k.find("jac")!=std::string::npos) {
        keys_order1.push_back(k);
        continue;
      }
      keys_order0.push_back(k);
    }

    // Print all NLP stats
    for (auto keys : {&keys_order0, &keys_order1, &keys_order2}) {
        std::sort(keys->begin(), keys->end());
        for (auto k : *keys) {
          const FStats& fs = m->fstats.at(k);
          if (fs.n_call!=0) {
            print(("%" + str(maxNameLen) + "s ").c_str(), k.c_str());
            print_stats_line(fs.n_call, fs.t_proc, fs.t_wall);
          }
        }
    }

    // Sum the previously printed stats
    double t_wall_all_previous = 0;
    double t_proc_all_previous = 0;
    for (auto k : keys) {
      const FStats& fs = m->fstats.at(k);
      t_proc_all_previous += fs.t_proc;
      t_wall_all_previous += fs.t_wall;
    }
    print(("%" + str(maxNameLen) + "s ").c_str(), "all previous");
    print_stats_line(-1, t_proc_all_previous, t_wall_all_previous);

    // Sort and show the remainder of keys
    std::sort(keys_other.begin(), keys_other.end());
    for (std::string k : keys_other) {
      const FStats& fs = m->fstats.at(k);
      if (fs.n_call!=0) {
        print(("%" + str(maxNameLen) + "s ").c_str(), k.c_str());
        print_stats_line(fs.n_call, fs.t_proc, fs.t_wall);
      }
      t_proc_all_previous += fs.t_proc;
      t_wall_all_previous += fs.t_wall;
    }

    // Show the mainloop stats
    const FStats& fs_mainloop = m->fstats.at("mainloop");
    if (fs_mainloop.n_call>0) {
      print(("%" + str(maxNameLen) + "s ").c_str(), "solver");
      print_stats_line(-1,
        fs_mainloop.t_proc-t_proc_all_previous, fs_mainloop.t_wall-t_wall_all_previous);
      print(("%" + str(maxNameLen) + "s ").c_str(), "mainloop");
      print_stats_line(-1, fs_mainloop.t_proc, fs_mainloop.t_wall);
    }
  }

  Dict OracleFunction::get_stats(void *mem) const {
    auto m = static_cast<OracleMemory*>(mem);

    // Add timing statistics
    Dict stats;
    for (auto&& s : m->fstats) {
      stats["n_call_" +s.first] = s.second.n_call;
      stats["t_wall_" +s.first] = s.second.t_wall;
      stats["t_proc_" +s.first] = s.second.t_proc;
    }
    return stats;
  }

  int OracleFunction::init_mem(void* mem) const {
    if (!mem) return 1;
    auto m = static_cast<OracleMemory*>(mem);

    // Create statistics
    for (auto&& e : all_functions_) {
      m->fstats[e.first] = FStats();
    }
    return 0;
  }

  void OracleFunction::set_temp(void* mem, const double** arg, double** res,
                            int* iw, double* w) const {
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
      "No function \"" + name + "\" in " + this->name() + ". " +
      "Available functions: " + join(get_function()) + ".");
    return it->second.f;
  }

  bool OracleFunction::monitored(const std::string &name) const {
    auto it = all_functions_.find(name);
    casadi_assert(it!=all_functions_.end(),
      "No function \"" + name + "\" in " + this->name()+ ". " +
      "Available functions: " + join(get_function()) + ".");
    return it->second.monitored;
  }

  bool OracleFunction::has_function(const std::string& fname) const {
    return all_functions_.find(fname) != all_functions_.end();
  }


} // namespace casadi
