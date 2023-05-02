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


#include "fmu_function.hpp"
#include "casadi_misc.hpp"
#include "serializing_stream.hpp"
#include "dae_builder_internal.hpp"
#include "fmu_impl.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

#ifdef WITH_OPENMP
#include <omp.h>
#endif // WITH_OPENMP

#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.thread.h>
#else // CASADI_WITH_THREAD_MINGW
#include <thread>
#endif // CASADI_WITH_THREAD_MINGW
#endif // CASADI_WITH_THREAD

namespace casadi {

#ifdef WITH_FMU

int FmuFunction::init_mem(void* mem) const {
  casadi_assert(mem != 0, "Memory is null");
  // Instantiate base classes
  if (FunctionInternal::init_mem(mem)) return 1;
  // Number of memory instances needed
  casadi_int n_mem = std::max(static_cast<casadi_int>(1),
    std::max(max_jac_tasks_, max_hess_tasks_));
  // Initialize master and all slaves
  FmuMemory* m = static_cast<FmuMemory*>(mem);
  for (casadi_int i = 0; i < n_mem; ++i) {
    // Initialize the memory object itself or a slave
    FmuMemory* m1 = i == 0 ? m : m->slaves.at(i - 1);
    if (get_fmu(fmu_)->init_mem(m1)) return 1;
  }
  return 0;
}

void* FmuFunction::alloc_mem() const {
  // Create (master) memory object
  FmuMemory* m = new FmuMemory(*this);
  // Attach additional (slave) memory objects
  for (casadi_int i = 1; i < max_jac_tasks_; ++i) {
    m->slaves.push_back(new FmuMemory(*this));
  }
  return m;
}

void FmuFunction::free_mem(void *mem) const {
  // Consistency check
  casadi_assert(mem != nullptr, "Memory is null");
  FmuMemory* m = static_cast<FmuMemory*>(mem);
  // Free slave memory
  for (FmuMemory*& s : m->slaves) {
    if (!s) continue;
    // Free FMU memory
    if (s->instance) {
      get_fmu(fmu_)->free_instance(s->instance);
      s->instance = nullptr;
    }
    // Free the slave
    delete s;
  }
  // Free FMI memory
  if (m->instance) {
    get_fmu(fmu_)->free_instance(m->instance);
    m->instance = nullptr;
  }
  // Free the memory object
  delete m;
}

FmuFunction::FmuFunction(const std::string& name, const Fmu& fmu,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out)
    : FunctionInternal(name), fmu_(fmu) {
  // Parse input IDs
  in_.resize(name_in.size());
  for (size_t k = 0; k < name_in.size(); ++k) {
    try {
      in_[k] = InputStruct::parse(name_in[k], get_fmu(fmu));
    } catch (std::exception& e) {
      casadi_error("Cannot process input " + name_in[k] + ": " + std::string(e.what()));
    }
  }
  // Parse output IDs
  out_.resize(name_out.size());
  for (size_t k = 0; k < name_out.size(); ++k) {
    try {
      out_[k] = OutputStruct::parse(name_out[k], get_fmu(fmu));
    } catch (std::exception& e) {
      casadi_error("Cannot process output " + name_out[k] + ": " + std::string(e.what()));
    }
  }
  // Set input/output names
  name_in_ = name_in;
  name_out_ = name_out;
  // Default options
  enable_ad_ = get_fmu(fmu)->get_directional_derivative_ != 0;
  validate_ad_ = false;
  validate_ad_file_ = "";
  make_symmetric_ = true;
  check_hessian_ = false;
  enable_fd_op_ = enable_ad_ && !all_regular();  // Use FD for second and higher order derivatives
  step_ = 1e-6;
  abstol_ = 1e-3;
  reltol_ = 1e-3;
  print_progress_ = false;
  new_jacobian_ = true;
  new_hessian_ = true;
  hessian_coloring_ = true;
  parallelization_ = Parallelization::SERIAL;
  // Number of parallel tasks, by default
  max_n_tasks_ = 1;
  max_jac_tasks_ = max_hess_tasks_ = 0;
}

FmuFunction::~FmuFunction() {
  // Free memory
  clear_mem();
}

const Options FmuFunction::options_
= {{&FunctionInternal::options_},
   {{"scheme_in",
     {OT_STRINGVECTOR,
      "Names of the inputs in the scheme"}},
    {"scheme_out",
     {OT_STRINGVECTOR,
      "Names of the outputs in the scheme"}},
    {"scheme",
     {OT_DICT,
      "Definitions of the scheme variables"}},
    {"aux",
     {OT_STRINGVECTOR,
      "Auxilliary variables"}},
    {"enable_ad",
     {OT_BOOL,
      "Calculate first order derivatives using FMU directional derivative support"}},
    {"validate_ad",
     {OT_BOOL,
      "Compare analytic derivatives with finite differences for validation"}},
    {"validate_ad_file",
     {OT_STRING,
      "Redirect results of Hessian validation to a file instead of generating a warning"}},
    {"check_hessian",
     {OT_BOOL,
      "Symmetry check for Hessian"}},
    {"make_symmetric",
     {OT_BOOL,
      "Ensure Hessian is symmetric"}},
    {"step",
     {OT_DOUBLE,
      "Step size, scaled by nominal value"}},
    {"abstol",
     {OT_DOUBLE,
      "Absolute error tolerance, scaled by nominal value"}},
    {"reltol",
     {OT_DOUBLE,
      "Relative error tolerance"}},
    {"parallelization",
     {OT_STRING,
      "Parallelization [SERIAL|openmp|thread]"}},
    {"print_progress",
     {OT_BOOL,
      "Print progress during Jacobian/Hessian evaluation"}},
    {"new_jacobian",
     {OT_BOOL,
      "Use Jacobian implementation in class"}},
    {"new_hessian",
     {OT_BOOL,
      "Use Hessian implementation in class"}},
    {"hessian_coloring",
     {OT_BOOL,
      "Enable the use of graph coloring (star coloring) for Hessian calculation. "
      "Note that disabling the coloring can improve symmetry check diagnostics."}}
   }
};

void FmuFunction::init(const Dict& opts) {
  // Read options
  for (auto&& op : opts) {
    if (op.first=="enable_ad") {
      enable_ad_ = op.second;
    } else if (op.first=="validate_ad") {
      validate_ad_ = op.second;
    } else if (op.first=="validate_ad_file") {
      validate_ad_file_ = op.second.to_string();
    } else if (op.first=="check_hessian") {
      check_hessian_ = op.second;
    } else if (op.first=="make_symmetric") {
      make_symmetric_ = op.second;
    } else if (op.first=="step") {
      step_ = op.second;
    } else if (op.first=="abstol") {
      abstol_ = op.second;
    } else if (op.first=="reltol") {
      reltol_ = op.second;
    } else if (op.first=="parallelization") {
      parallelization_ = to_enum<Parallelization>(op.second, "serial");
    } else if (op.first=="print_progress") {
      print_progress_ = op.second;
    } else if (op.first=="new_jacobian") {
      new_jacobian_ = op.second;
    } else if (op.first=="new_hessian") {
      new_hessian_ = op.second;
    } else if (op.first=="hessian_coloring") {
      hessian_coloring_ = op.second;
    }
  }

  // Call the initialization method of the base class
  FunctionInternal::init(opts);

  // Read FD mode
  fd_ = to_enum<FdMode>(fd_method_, "forward");

  // Consistency checks
  if (enable_ad_) casadi_assert(get_fmu(fmu_)->get_directional_derivative_ != nullptr,
    "FMU does not provide support for analytic derivatives");
  if (validate_ad_ && !enable_ad_) casadi_error("Inconsistent options");

  // New AD validation file, if any
  if (!validate_ad_file_.empty()) {
    std::ofstream valfile;
    valfile.open(validate_ad_file_);
    valfile << "Output Input Value Nominal Min Max AD FD Step Offset Stencil" << std::endl;
  }
  // Which inputs and outputs exist
  has_fwd_ = has_adj_ = has_jac_ = has_hess_ = false;
  for (auto&& i : out_) {
    switch (i.type) {
    case OutputType::JAC:
    case OutputType::JAC_TRANS:
      has_jac_ = true;
      break;
    case OutputType::FWD:
      has_fwd_ = true;
      break;
    case OutputType::ADJ:
      has_adj_ = true;
      break;
    case OutputType::HESS:
      has_adj_ = true;
      has_hess_ = true;
    default:
      break;
    }
  }

  // Forward derivatives not yet implemented
  if (has_fwd_) casadi_warning("Forward derivatives not implemented, ignored");


  // Quick return if no Jacobian calculation
  if (!has_jac_ && !has_adj_ && !has_hess_) return;

  // Parallelization
  switch (parallelization_) {
    case Parallelization::SERIAL:
      if (verbose_) casadi_message("Serial evaluation");
      break;
#ifdef WITH_OPENMP
    case Parallelization::OPENMP:
      max_n_tasks_ = omp_get_max_threads();
      if (verbose_) casadi_message("OpenMP using at most " + str(max_n_tasks_) + " threads");
      break;
#endif // WITH_OPENMP
#ifdef CASADI_WITH_THREAD
    case Parallelization::THREAD:
      max_n_tasks_ = std::thread::hardware_concurrency();
      if (verbose_) casadi_message("std::thread using at most " + str(max_n_tasks_) + " threads");
      break;
#endif // CASADI_WITH_THREAD
    default:
      casadi_warning("Parallelization " + to_string(parallelization_)
        + " not enabled during compilation. Falling back to serial evaluation");
      parallelization_ = Parallelization::SERIAL;
      break;
  }

  // Collect all inputs in any Jacobian, Hessian or adjoint block
  std::vector<size_t> in_jac(get_fmu(fmu_)->iind_.size(), 0);
  jac_in_.clear();
  jac_nom_in_.clear();
  for (auto&& i : out_) {
    if (i.type == OutputType::JAC || i.type == OutputType::JAC_TRANS
        || i.type == OutputType::ADJ || i.type == OutputType::HESS) {
      // Get input indices
      const std::vector<size_t>& iind = get_fmu(fmu_)->ired_.at(i.wrt);
      // Skip if no entries
      if (iind.empty()) continue;
      // Consistency check
      bool exists = in_jac[iind.front()] > 0;
      for (size_t j : iind) casadi_assert((in_jac[j] > 0) == exists, "Jacobian not a block");
      // Add selection
      if (!exists) {
        for (size_t j : iind) {
          jac_in_.push_back(j);
          jac_nom_in_.push_back(get_fmu(fmu_)->nominal_in_[j]);
          in_jac[j] = jac_in_.size();
        }
      }
      // Add column interval
      i.cbegin = in_jac[iind.front()] - 1;
      i.cend = i.cbegin + iind.size();
      // Also rows for Hessian blocks
      if (i.type == OutputType::HESS) {
        // Get input indices
        const std::vector<size_t>& iind = get_fmu(fmu_)->ired_.at(i.ind);
        // Skip if no entries
        if (iind.empty()) continue;
        // Consistency check
        bool exists = in_jac[iind.front()] > 0;
        for (size_t j : iind) casadi_assert((in_jac[j] > 0) == exists, "Hessian not a block");
        // Add selection
        if (!exists) {
          for (size_t j : iind) {
            jac_in_.push_back(j);
            jac_nom_in_.push_back(get_fmu(fmu_)->nominal_in_[j]);
            in_jac[j] = jac_in_.size();
          }
        }
        // Add column interval
        i.rbegin = in_jac[iind.front()] - 1;
        i.rend = i.rbegin + iind.size();
      }
    }
  }

  // Transpose of Jacobian sparsity
  sp_trans_map_.resize(out_.size(), -1);
  sp_trans_.clear();

  // Collect all outputs in any Jacobian or adjoint block
  in_jac.resize(get_fmu(fmu_)->oind_.size());
  std::fill(in_jac.begin(), in_jac.end(), 0);
  jac_out_.clear();
  for (size_t k = 0; k < out_.size(); ++k) {
    OutputStruct& i = out_[k];
    if (i.type == OutputType::JAC || i.type == OutputType::JAC_TRANS) {
      // Get output indices
      const std::vector<size_t>& oind = get_fmu(fmu_)->ored_.at(i.ind);
      // Skip if no entries
      if (oind.empty()) continue;
      // Consistency check
      bool exists = in_jac[oind.front()] > 0;
      for (size_t j : oind) casadi_assert((in_jac[j] > 0) == exists, "Jacobian not a block");
      // Add selection
      if (!exists) {
        for (size_t j : oind) {
          jac_out_.push_back(j);
          in_jac[j] = jac_out_.size();
        }
      }
      // Add row interval
      i.rbegin = in_jac[oind.front()] - 1;
      i.rend = i.rbegin + oind.size();
      // Additional memory for transpose
      if (i.type == OutputType::JAC_TRANS) {
        // Retrieve the sparsity pattern
        const Sparsity& sp = sparsity_out(k);
        // Save transpose of sparsity pattern
        sp_trans_map_.at(k) = sp_trans_.size();
        sp_trans_.push_back(sp.T());
        // Work vectors for casadi_trans
        alloc_w(sp.nnz());
        alloc_iw(sp.size2());
      }
    }
  }
  for (auto&& i : in_) {
    if (i.type == InputType::ADJ) {
      // Get output indices
      const std::vector<size_t>& oind = get_fmu(fmu_)->ored_.at(i.ind);
      // Skip if no entries
      if (oind.empty()) continue;
      // Consistency check
      bool exists = in_jac[oind.front()] > 0;
      for (size_t j : oind) casadi_assert((in_jac[j] > 0) == exists, "Jacobian not a block");
      // Add selection
      if (!exists) {
        for (size_t j : oind) {
          jac_out_.push_back(j);
          in_jac[j] = jac_out_.size();
        }
      }
    }
  }

  // Get sparsity pattern for extended Jacobian
  jac_sp_ = get_fmu(fmu_)->jac_sparsity(jac_out_, jac_in_);

  // Calculate graph coloring
  jac_colors_ = jac_sp_.uni_coloring();
  if (verbose_) casadi_message("Jacobian graph coloring: " + str(jac_sp_.size2())
    + " -> " + str(jac_colors_.size2()) + " directions");

  // Setup Jacobian memory
  casadi_jac_setup(&p_, jac_sp_, jac_colors_);
  p_.nom_in = get_ptr(jac_nom_in_);
  p_.map_out = get_ptr(jac_out_);
  p_.map_in = get_ptr(jac_in_);

  // Do not use more threads than there are colors in the Jacobian
  max_jac_tasks_ = std::min(max_n_tasks_, jac_colors_.size2());

  // Work vector for storing extended Jacobian, shared between threads
  if (has_jac_) {
    alloc_w(jac_sp_.nnz(), true);  // jac_nz
  }

  // Work vectors for adjoint derivative calculation, shared between threads
  if (has_adj_) {
    alloc_w(get_fmu(fmu_)->oind_.size(), true);  // aseed
    alloc_w(get_fmu(fmu_)->iind_.size(), true);  // asens
  }

  // If Hessian calculation is needed
  if (has_hess_) {
    // Get sparsity pattern for extended Hessian
    hess_sp_ = get_fmu(fmu_)->hess_sparsity(jac_in_, jac_in_);
    casadi_assert(hess_sp_.size1() == jac_in_.size(), "Inconsistent Hessian dimensions");
    casadi_assert(hess_sp_.size2() == jac_in_.size(), "Inconsistent Hessian dimensions");
    const casadi_int *hess_row = hess_sp_.row();
    casadi_int hess_nnz = hess_sp_.nnz();

    // Get nonlinearly entering variables
    std::vector<bool> is_nonlin(jac_in_.size(), false);
    for (casadi_int k = 0; k < hess_nnz; ++k) is_nonlin[hess_row[k]] = true;
    nonlin_.clear();
    for (casadi_int c = 0; c < jac_in_.size(); ++c) {
      if (is_nonlin[c]) nonlin_.push_back(c);
    }
    // Calculate graph coloring
    if (hessian_coloring_) {
      // Star coloring
      hess_colors_ = hess_sp_.star_coloring();
      if (verbose_) casadi_message("Hessian graph coloring: " + str(nonlin_.size())
        + " -> " + str(hess_colors_.size2()) + " directions");
      // Indices of non-nonlinearly entering variables
      std::vector<casadi_int> zind;
      for (casadi_int c = 0; c < jac_in_.size(); ++c) {
        if (!is_nonlin[c]) zind.push_back(c);
      }
      // Zero out corresponding rows
      hess_colors_.erase(zind, range(hess_colors_.size2()));
    } else {
      // One color for each nonlinear variable
      hess_colors_ = Sparsity(jac_in_.size(), nonlin_.size(), range(nonlin_.size() + 1), nonlin_);
      if (verbose_) casadi_message("Hessian calculation for " + str(nonlin_.size()) + " variables");
    }

    // Number of threads to be used for Hessian calculation
    max_hess_tasks_ = std::min(max_n_tasks_, hess_colors_.size2());

    // Work vector for storing extended Hessian, shared between threads
    alloc_w(hess_sp_.nnz(), true);  // hess_nz

    // Work vector for perturbed adjoint sensitivities
    alloc_w(max_hess_tasks_ * get_fmu(fmu_)->iind_.size(), true);  // pert_asens

    // Work vector for avoiding conflicting assignments for star coloring
    alloc_iw(max_hess_tasks_ * get_fmu(fmu_)->iind_.size(), true);  // star_iw

    // Work vector for making symmetric or checking symmetry
    alloc_iw(hess_sp_.size2());
  }

  // Total number of threads used for Jacobian/adjoint calculation
  // Note: Adjoint calculation also used for Hessian
  max_n_tasks_ = std::min(max_n_tasks_, std::max(max_jac_tasks_, max_hess_tasks_));
  if (verbose_) casadi_message("Jacobian calculation with " + str(max_n_tasks_) + " threads");

  // Work vectors for Jacobian/adjoint/Hessian calculation, for each thread
  casadi_int jac_iw, jac_w;
  casadi_jac_work(&p_, &jac_iw, &jac_w);
  alloc_iw(max_n_tasks_ * jac_iw, true);
  alloc_w(max_n_tasks_ * jac_w, true);
}

void FmuFunction::identify_io(
    std::vector<std::string>* scheme_in,
    std::vector<std::string>* scheme_out,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out) {
  // Clear returns
  if (scheme_in) scheme_in->clear();
  if (scheme_out) scheme_out->clear();
  // Parse FmuFunction inputs
  for (const std::string& n : name_in) {
    try {
      (void)InputStruct::parse(n, 0, scheme_in, scheme_out);
    } catch (std::exception& e) {
      casadi_error("Cannot process input " + n + ": " + std::string(e.what()));
    }
  }
  // Parse FmuFunction outputs
  for (const std::string& n : name_out) {
    try {
      (void)OutputStruct::parse(n, 0, scheme_in, scheme_out);
    } catch (std::exception& e) {
      casadi_error("Cannot process output " + n + ": " + std::string(e.what()));
    }
  }
  // Remove duplicates in scheme_in, also sorts alphabetically
  if (scheme_in) {
    std::set<std::string> s(scheme_in->begin(), scheme_in->end());
    scheme_in->assign(s.begin(), s.end());
  }
  // Remove duplicates in scheme_out, also sorts alphabetically
  if (scheme_out) {
    std::set<std::string> s(scheme_out->begin(), scheme_out->end());
    scheme_out->assign(s.begin(), s.end());
  }
}

InputStruct InputStruct::parse(const std::string& n, const Fmu2* fmu,
    std::vector<std::string>* name_in, std::vector<std::string>* name_out) {
  // Return value
  InputStruct s;
  // Look for a prefix
  if (has_prefix(n)) {
    // Get the prefix
    std::string pref, rem;
    pref = pop_prefix(n, &rem);
    if (pref == "out") {
      if (has_prefix(rem)) {
        // Second order function output (unused): Get the prefix
        pref = pop_prefix(rem, &rem);
        if (pref == "adj") {
          s.type = InputType::ADJ_OUT;
          s.ind = fmu ? fmu->index_in(rem) : -1;
          if (name_in) name_in->push_back(rem);
        } else {
          casadi_error("Cannot process: " + n);
        }
      } else {
        // Nondifferentiated function output (unused)
        s.type = InputType::OUT;
        s.ind = fmu ? fmu->index_out(rem) : -1;
        if (name_out) name_out->push_back(rem);
      }
    } else if (pref == "fwd") {
      // Forward seed
      s.type = InputType::FWD;
      s.ind = fmu ? fmu->index_in(rem) : 0;
      if (name_in) name_in->push_back(rem);
    } else if (pref == "adj") {
      // Adjoint seed
      s.type = InputType::ADJ;
      s.ind = fmu ? fmu->index_out(rem) : 0;
      if (name_out) name_out->push_back(rem);
    } else {
      // No such prefix
      casadi_error("No such prefix: " + pref);
    }
  } else {
    // No prefix - regular input
    s.type = InputType::REG;
    s.ind = fmu ? fmu->index_in(n) : 0;
    if (name_in) name_in->push_back(n);
  }
  // Return input struct
  return s;
}

OutputStruct OutputStruct::parse(const std::string& n, const Fmu2* fmu,
    std::vector<std::string>* name_in, std::vector<std::string>* name_out) {
  // Return value
  OutputStruct s;
  // Look for prefix
  if (has_prefix(n)) {
    // Get the prefix
    std::string pref, rem;
    pref = pop_prefix(n, &rem);
    if (pref == "jac") {
      // Jacobian block
      casadi_assert(has_prefix(rem), "Two arguments expected for Jacobian block");
      pref = pop_prefix(rem, &rem);
      if (pref == "adj") {
        // Jacobian of adjoint sensitivity
        casadi_assert(has_prefix(rem), "Two arguments expected for Jacobian block");
        pref = pop_prefix(rem, &rem);
        if (has_prefix(rem)) {
          // Jacobian with respect to a sensitivity seed
          std::string sens = pref;
          pref = pop_prefix(rem, &rem);
          if (pref == "adj") {
            // Jacobian of adjoint sensitivity w.r.t. adjoint seed -> Transpose of Jacobian
            s.type = OutputType::JAC_TRANS;
            s.ind = fmu ? fmu->index_out(rem) : -1;
            if (name_out) name_out->push_back(rem);
            s.wrt = fmu ? fmu->index_in(sens) : -1;
            if (name_in) name_in->push_back(sens);
          } else if (pref == "out") {
            // Jacobian w.r.t. to dummy output
            s.type = OutputType::JAC_ADJ_OUT;
            s.ind = fmu ? fmu->index_in(sens) : -1;
            if (name_in) name_in->push_back(sens);
            s.wrt = fmu ? fmu->index_out(rem) : -1;
            if (name_in) name_out->push_back(rem);
          } else {
            casadi_error("No such prefix: " + pref);
          }
        } else {
          // Hessian output
          s.type = OutputType::HESS;
          s.ind = fmu ? fmu->index_in(pref) : -1;
          if (name_in) name_in->push_back(pref);
          s.wrt = fmu ? fmu->index_in(rem) : -1;
          if (name_in) name_in->push_back(rem);
        }
      } else {
        if (has_prefix(rem)) {
          std::string out = pref;
          pref = pop_prefix(rem, &rem);
          if (pref == "adj") {
            // Jacobian of regular output w.r.t. adjoint sensitivity seed
            s.type = OutputType::JAC_REG_ADJ;
            s.ind = fmu ? fmu->index_out(out) : -1;
            if (name_out) name_out->push_back(out);
            s.wrt = fmu ? fmu->index_out(rem) : -1;
            if (name_out) name_out->push_back(rem);
          } else {
            casadi_error("No such prefix: " + pref);
          }
        } else {
          // Regular Jacobian
          s.type = OutputType::JAC;
          s.ind = fmu ? fmu->index_out(pref) : -1;
          if (name_out) name_out->push_back(pref);
          s.wrt = fmu ? fmu->index_in(rem) : -1;
          if (name_in) name_in->push_back(rem);
        }
      }
    } else if (pref == "fwd") {
      // Forward sensitivity
      s.type = OutputType::FWD;
      s.ind = fmu ? fmu->index_out(rem) : -1;
      if (name_out) name_out->push_back(rem);
    } else if (pref == "adj") {
      // Adjoint sensitivity
      s.type = OutputType::ADJ;
      s.wrt = fmu ? fmu->index_in(rem) : -1;
      if (name_in) name_in->push_back(rem);
    } else {
      // No such prefix
      casadi_error("No such prefix: " + pref);
    }
  } else {
    // No prefix - regular output
    s.type = OutputType::REG;
    s.ind = fmu ? fmu->index_out(n) : -1;
    if (name_out) name_out->push_back(n);
  }
  // Return output struct
  return s;
}

Sparsity FmuFunction::get_sparsity_in(casadi_int i) {
  switch (in_.at(i).type) {
    case InputType::REG:
      return Sparsity::dense(get_fmu(fmu_)->ired_.at(in_.at(i).ind).size(), 1);
    case InputType::FWD:
      return Sparsity::dense(get_fmu(fmu_)->ired_.at(in_.at(i).ind).size(), 1);
    case InputType::ADJ:
      return Sparsity::dense(get_fmu(fmu_)->ored_.at(in_.at(i).ind).size(), 1);
    case InputType::OUT:
      return Sparsity(get_fmu(fmu_)->ored_.at(in_.at(i).ind).size(), 1);
    case InputType::ADJ_OUT:
      return Sparsity(get_fmu(fmu_)->ired_.at(in_.at(i).ind).size(), 1);
  }
  return Sparsity();
}

Sparsity FmuFunction::get_sparsity_out(casadi_int i) {
  const OutputStruct& s = out_.at(i);
  switch (out_.at(i).type) {
    case OutputType::REG:
      return Sparsity::dense(get_fmu(fmu_)->ored_.at(s.ind).size(), 1);
    case OutputType::FWD:
      return Sparsity::dense(get_fmu(fmu_)->ored_.at(s.ind).size(), 1);
    case OutputType::ADJ:
      return Sparsity::dense(get_fmu(fmu_)->ired_.at(s.wrt).size(), 1);
    case OutputType::JAC:
      return get_fmu(fmu_)->jac_sparsity(s.ind, s.wrt);
    case OutputType::JAC_TRANS:
      return get_fmu(fmu_)->jac_sparsity(s.ind, s.wrt).T();
    case OutputType::JAC_ADJ_OUT:
      return Sparsity(get_fmu(fmu_)->ired_.at(s.ind).size(), get_fmu(fmu_)->ored_.at(s.wrt).size());
    case OutputType::JAC_REG_ADJ:
      return Sparsity(get_fmu(fmu_)->ored_.at(s.ind).size(), get_fmu(fmu_)->ored_.at(s.wrt).size());
    case OutputType::HESS:
      return get_fmu(fmu_)->hess_sparsity(s.ind, s.wrt);
  }
  return Sparsity();
}

std::vector<double> FmuFunction::get_nominal_in(casadi_int i) const {
  switch (in_.at(i).type) {
    case InputType::REG:
      return get_fmu(fmu_)->get_nominal_in(in_.at(i).ind);
    case InputType::FWD:
      break;
    case InputType::ADJ:
      break;
    case InputType::OUT:
      return get_fmu(fmu_)->get_nominal_out(in_.at(i).ind);
    case InputType::ADJ_OUT:
      break;
  }
  // Default: Base class
  return FunctionInternal::get_nominal_in(i);
}

std::vector<double> FmuFunction::get_nominal_out(casadi_int i) const {
  switch (out_.at(i).type) {
    case OutputType::REG:
      return get_fmu(fmu_)->get_nominal_out(out_.at(i).ind);
    case OutputType::FWD:
      break;
    case OutputType::ADJ:
      break;
    case OutputType::JAC:
      casadi_warning("FmuFunction::get_nominal_out not implemented for OutputType::JAC");
      break;
    case OutputType::JAC_TRANS:
      casadi_warning("FmuFunction::get_nominal_out not implemented for OutputType::JAC_TRANS");
      break;
    case OutputType::JAC_ADJ_OUT:
      casadi_warning("FmuFunction::get_nominal_out not implemented for OutputType::JAC_ADJ_OUT");
      break;
    case OutputType::JAC_REG_ADJ:
      casadi_warning("FmuFunction::get_nominal_out not implemented for OutputType::JAC_REG_ADJ");
      break;
    case OutputType::HESS:
      casadi_warning("FmuFunction::get_nominal_out not implemented for OutputType::HESS");
      break;
  }
  // Default: Base class
  return FunctionInternal::get_nominal_out(i);
}

int FmuFunction::eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // Get memory struct
  FmuMemory* m = static_cast<FmuMemory*>(mem);
  casadi_assert(m != 0, "Memory is null");
  // What blocks are there?
  bool need_jac = false, need_adj = false, need_hess = false;
  for (size_t k = 0; k < out_.size(); ++k) {
    if (res[k]) {
      switch (out_[k].type) {
        case OutputType::JAC:
        case OutputType::JAC_TRANS:
          need_jac = true;
          break;
        case OutputType::ADJ:
          need_adj = true;
          break;
        case OutputType::HESS:
          need_adj = true;
          need_hess = true;
          break;
        default:
          break;
      }
    }
  }
  // Work vectors, shared between threads
  double *aseed = 0, *asens = 0, *jac_nz = 0, *hess_nz = 0;
  if (need_jac) {
    // Jacobian nonzeros, initialize to NaN
    jac_nz = w; w += jac_sp_.nnz();
    std::fill(jac_nz, jac_nz + jac_sp_.nnz(), casadi::nan);
  }
  if (need_adj) {
    // Set up vectors
    aseed = w; w += get_fmu(fmu_)->oind_.size();
    asens = w; w += get_fmu(fmu_)->iind_.size();
    // Clear seed/sensitivity vectors
    std::fill(aseed, aseed + get_fmu(fmu_)->oind_.size(), 0);
    std::fill(asens, asens + get_fmu(fmu_)->iind_.size(), 0);
    // Copy adjoint seeds to aseed
    for (size_t i = 0; i < in_.size(); ++i) {
      if (arg[i] && in_[i].type == InputType::ADJ) {
        const std::vector<size_t>& oind = get_fmu(fmu_)->ored_[in_[i].ind];
        for (size_t k = 0; k < oind.size(); ++k) aseed[oind[k]] = arg[i][k];
      }
    }
  }
  if (need_hess) {
    // Hessian nonzeros, initialize to NaN
    hess_nz = w; w += hess_sp_.nnz();
    std::fill(hess_nz, hess_nz + hess_sp_.nnz(), casadi::nan);
  }
  // Setup memory for threads
  for (casadi_int task = 0; task < max_n_tasks_; ++task) {
    FmuMemory* s = task == 0 ? m : m->slaves.at(task - 1);
    // Shared memory
    s->arg = arg;
    s->res = res;
    s->aseed = aseed;
    s->asens = asens;
    s->jac_nz = jac_nz;
    s->hess_nz = hess_nz;
    // Thread specific memory
    casadi_jac_init(&p_, &s->d, &iw, &w);
    if (task < max_hess_tasks_) {
      // Perturbed adjoint sensitivities
      s->pert_asens = w;
      w += get_fmu(fmu_)->iind_.size();
      // Work vector for avoiding assignment of illegal nonzeros
      s->star_iw = iw;
      iw += get_fmu(fmu_)->iind_.size();
    }
  }
  // Evaluate everything except Hessian, possibly in parallel
  if (verbose_) casadi_message("Evaluating regular outputs, forward sens, extended Jacobian");
  if (eval_all(m, max_jac_tasks_, true, need_jac, need_adj, false)) return 1;
  // Evaluate Hessian
  if (need_hess) {
    if (verbose_) casadi_message("Evaluating extended Hessian");
    if (eval_all(m, max_hess_tasks_, false, false, false, true)) return 1;
    // Post-process Hessian
    remove_nans(hess_nz, iw);
    if (check_hessian_) check_hessian(m, hess_nz, iw);
    if (make_symmetric_) make_symmetric(hess_nz, iw);
  }
  // Fetch calculated blocks
  for (size_t k = 0; k < out_.size(); ++k) {
    // Get nonzeros, skip if not needed
    double* r = res[k];
    if (!r) continue;
    // Get by type
    switch (out_[k].type) {
      case OutputType::JAC:
        casadi_get_sub(r, jac_sp_, jac_nz,
          out_[k].rbegin, out_[k].rend, out_[k].cbegin, out_[k].cend);
        break;
      case OutputType::JAC_TRANS:
        casadi_get_sub(w, jac_sp_, jac_nz,
          out_[k].rbegin, out_[k].rend, out_[k].cbegin, out_[k].cend);
        casadi_trans(w, sp_trans_[sp_trans_map_[k]], r, sparsity_out(k), iw);
        break;
      case OutputType::ADJ:
        for (size_t id : get_fmu(fmu_)->ired_[out_[k].wrt]) *r++ = asens[id];
        break;
      case OutputType::HESS:
        casadi_get_sub(r, hess_sp_, hess_nz,
          out_[k].rbegin, out_[k].rend, out_[k].cbegin, out_[k].cend);
        break;
      default:
        break;
    }
  }
  // Successful return
  return 0;
}

int FmuFunction::eval_all(FmuMemory* m, casadi_int n_task,
    bool need_nondiff, bool need_jac, bool need_adj, bool need_hess) const {
  // Return flag
  int flag = 0;
  // Evaluate, serially or in parallel
  if (parallelization_ == Parallelization::SERIAL || n_task == 1
      || (!need_jac && !need_adj && !need_hess)) {
    // Evaluate serially
    flag = eval_task(m, 0, 1, need_nondiff, need_jac, need_adj, need_hess);
  } else if (parallelization_ == Parallelization::OPENMP) {
    #ifdef WITH_OPENMP
    // Parallel region
    #pragma omp parallel reduction(||:flag)
    {
      // Get thread number
      casadi_int task = omp_get_thread_num();
      // Get number of threads in region
      casadi_int num_threads = omp_get_num_threads();
      // Number of threads that are actually used
      casadi_int num_used_threads = std::min(num_threads, n_task);
      // Evaluate in parallel
      if (task < num_used_threads) {
        FmuMemory* s = task == 0 ? m : m->slaves.at(task - 1);
        flag = eval_task(s, task, num_used_threads,
          need_nondiff && task == 0, need_jac, need_adj, need_hess);
      } else {
        // Nothing to do for thread
        flag = 0;
      }
    }
    #else   // WITH_OPENMP
    flag = 1;
    #endif  // WITH_OPENMP
  } else if (parallelization_ == Parallelization::THREAD) {
    #ifdef CASADI_WITH_THREAD
    // Return value for each thread
    std::vector<int> flag_task(n_task);
    // Spawn threads
    std::vector<std::thread> threads;
    for (casadi_int task = 0; task < n_task; ++task) {
      threads.emplace_back(
        [&, task](int* fl) {
          FmuMemory* s = task == 0 ? m : m->slaves.at(task - 1);
          *fl = eval_task(s, task, n_task,
            need_nondiff && task == 0, need_jac, need_adj, need_hess);
        }, &flag_task[task]);
    }
    // Join threads
    for (auto&& th : threads) th.join();
    // Join return flags
    for (int fl : flag_task) flag = flag || fl;
    #else   // CASADI_WITH_THREAD
    flag = 1;
    #endif  // CASADI_WITH_THREAD
  } else {
    casadi_error("Unknown parallelization: " + to_string(parallelization_));
  }
  // Return combined error flag
  return flag;
}

int FmuFunction::eval_task(FmuMemory* m, casadi_int task, casadi_int n_task,
    bool need_nondiff, bool need_jac, bool need_adj, bool need_hess) const {
  // Pass all regular inputs
  for (size_t k = 0; k < in_.size(); ++k) {
    if (in_[k].type == InputType::REG) {
      get_fmu(fmu_)->set(m, in_[k].ind, m->arg[k]);
    }
  }
  // Request all regular outputs to be evaluated
  for (size_t k = 0; k < out_.size(); ++k) {
    if (m->res[k] && out_[k].type == OutputType::REG) {
      get_fmu(fmu_)->request(m, out_[k].ind);
    }
  }
  // Evaluate
  if (get_fmu(fmu_)->eval(m)) return 1;
  // Get regular outputs (master thread only)
  if (need_nondiff) {
    for (size_t k = 0; k < out_.size(); ++k) {
      if (m->res[k] && out_[k].type == OutputType::REG) {
        get_fmu(fmu_)->get(m, out_[k].ind, m->res[k]);
      }
    }
  }
  // Evalute extended Jacobian
  if (need_jac || need_adj) {
    // Selection of colors to be evaluated for the thread
    casadi_int c_begin = (task * jac_colors_.size2()) / n_task;
    casadi_int c_end = ((task + 1) * jac_colors_.size2()) / n_task;
    // Loop over colors
    for (casadi_int c = c_begin; c < c_end; ++c) {
      // Print progress
      if (print_progress_) print("Jacobian calculation, thread %d/%d: Seeding variable %d/%d\n",
        task + 1, n_task, c - c_begin + 1, c_end - c_begin);
      // Get derivative directions
      casadi_jac_pre(&p_, &m->d, c);
      // Calculate derivatives
      get_fmu(fmu_)->set_seed(m, m->d.nseed, m->d.iseed, m->d.seed);
      get_fmu(fmu_)->request_sens(m, m->d.nsens, m->d.isens, m->d.wrt);
      if (get_fmu(fmu_)->eval_derivative(m, true)) return 1;
      get_fmu(fmu_)->get_sens(m, m->d.nsens, m->d.isens, m->d.sens);
      // Scale derivatives
      casadi_jac_scale(&p_, &m->d);
      // Collect Jacobian nonzeros
      if (need_jac) {
        for (casadi_int i = 0; i < m->d.nsens; ++i) {
          m->jac_nz[m->d.nzind[i]] = m->d.sens[i];
        }
      }
      // Propagate adjoint sensitivities
      if (need_adj) {
        for (casadi_int i = 0; i < m->d.nsens; ++i)
          m->asens[m->d.wrt[i]] += m->aseed[m->d.isens[i]] * m->d.sens[i];
      }
    }
  }
  // Evaluate extended Hessian
  if (need_hess) {
    // Hessian coloring
    casadi_int n_hc = hess_colors_.size2();
    const casadi_int *hc_colind = hess_colors_.colind(), *hc_row = hess_colors_.row();
    // Hessian sparsity
    const casadi_int *hess_colind = hess_sp_.colind(), *hess_row = hess_sp_.row();
    // Selection of colors to be evaluated for the thread
    casadi_int c_begin = (task * n_hc) / n_task;
    casadi_int c_end = ((task + 1) * n_hc) / n_task;
    // Unperturbed values, step size
    std::vector<double> x, h;
    // Loop over colors
    for (casadi_int c = c_begin; c < c_end; ++c) {
      // Print progress
      if (print_progress_) print("Hessian calculation, thread %d/%d: Seeding variable %d/%d\n",
        task + 1, n_task, c - c_begin + 1, c_end - c_begin);
      // Variables being seeded
      casadi_int v_begin = hc_colind[c];
      casadi_int v_end = hc_colind[c + 1];
      casadi_int nv = v_end - v_begin;
      // Loop over variables being seeded for color
      x.resize(nv);
      h.resize(nv);
      for (casadi_int v = 0; v < nv; ++v) {
        // Corresponding input in Fmu
        casadi_int ind1 = hc_row[v_begin + v];
        casadi_int id = jac_in_.at(ind1);
        // Get unperturbed value
        x[v] = m->ibuf_.at(id);
        // Step size
        h[v] = m->self.step_ * get_fmu(fmu_)->nominal_in_.at(id);
        // Make sure a a forward step remains in bounds
        if (x[v] + h[v] > get_fmu(fmu_)->max_in_.at(id)) {
          // Ensure a negative step is possible
          if (m->ibuf_.at(id) - h[v] < get_fmu(fmu_)->min_in_.at(id)) {
            std::stringstream ss;
            ss << "Cannot perturb " << get_fmu(fmu_)->vn_in_.at(id) << " at " << x[v]
              << " with step size "
              << m->self.step_ << ", nominal " << get_fmu(fmu_)->nominal_in_.at(id) << " min "
              << get_fmu(fmu_)->min_in_.at(id) << ", max " << get_fmu(fmu_)->max_in_.at(id);
            casadi_warning(ss.str());
            return 1;
          }
          // Take reverse step instead
          h[v] = -h[v];
        }
        // Perturb the input
        m->ibuf_.at(id) += h[v];
        m->changed_.at(id) = true;
        // Inverse of step size
        h[v] = 1. / h[v];
      }
      // Request all outputs
      for (size_t i : jac_out_) {
        m->requested_.at(i) = true;
        m->wrt_.at(i) = -1;
      }
      // Calculate perturbed inputs
      if (get_fmu(fmu_)->eval(m)) return 1;
      // Clear perturbed adjoint sensitivities
      std::fill(m->pert_asens, m->pert_asens + get_fmu(fmu_)->iind_.size(), 0);
      // Loop over colors of the Jacobian
      for (casadi_int c1 = 0; c1 < jac_colors_.size2(); ++c1) {
       // Get derivative directions
       casadi_jac_pre(&p_, &m->d, c1);
       // Calculate derivatives
       get_fmu(fmu_)->set_seed(m, m->d.nseed, m->d.iseed, m->d.seed);
       get_fmu(fmu_)->request_sens(m, m->d.nsens, m->d.isens, m->d.wrt);
       if (get_fmu(fmu_)->eval_derivative(m, true)) return 1;
       get_fmu(fmu_)->get_sens(m, m->d.nsens, m->d.isens, m->d.sens);
       // Scale derivatives
       casadi_jac_scale(&p_, &m->d);
       // Propagate adjoint sensitivities
       for (casadi_int i = 0; i < m->d.nsens; ++i)
         m->pert_asens[m->d.wrt[i]] += m->aseed[m->d.isens[i]] * m->d.sens[i];
      }
      // Count how many times each input is calculated
      std::fill(m->star_iw, m->star_iw + get_fmu(fmu_)->iind_.size(), 0);
      for (casadi_int v = 0; v < nv; ++v) {
        casadi_int ind1 = hc_row[v_begin + v];
        for (casadi_int k = hess_colind[ind1]; k < hess_colind[ind1 + 1]; ++k) {
          casadi_int id2 = jac_in_.at(hess_row[k]);
          m->star_iw[id2]++;
        }
      }
      // Loop over variables being seeded for color
      for (casadi_int v = 0; v < nv; ++v) {
        // Corresponding input in Fmu
        casadi_int ind1 = hc_row[v_begin + v];
        casadi_int id = jac_in_.at(ind1);
        // Restore input
        m->ibuf_.at(id) = x[v];
        m->changed_.at(id) = true;
        // Get column in Hessian
        for (casadi_int k = hess_colind[ind1]; k < hess_colind[ind1 + 1]; ++k) {
          casadi_int id2 = jac_in_.at(hess_row[k]);
          // Save
          if (m->star_iw[id2] > 1) {
            // Input gets perturbed by multiple variables for the same color:
            // Replace with NaN to use mirror element instead
            m->hess_nz[k] = casadi::nan;
          } else {
            // Finite difference approximation
            m->hess_nz[k] = h[v] * (m->pert_asens[id2] - m->asens[id2]);
          }
        }
      }
    }
  }
  // Successful return
  return 0;
}

void FmuFunction::check_hessian(FmuMemory* m, const double *hess_nz, casadi_int* iw) const {
  // Get Hessian sparsity pattern
  casadi_int n = hess_sp_.size1();
  const casadi_int *colind = hess_sp_.colind(), *row = hess_sp_.row();
  // Nonzero counters for transpose
  casadi_copy(colind, n, iw);
  // Loop over Hessian columns
  for (casadi_int c = 0; c < n; ++c) {
    // Loop over nonzeros for the column
    for (casadi_int k = colind[c]; k < colind[c + 1]; ++k) {
      // Get row of Hessian
      casadi_int r = row[k];
      // Get nonzero of transpose
      casadi_int k_tr = iw[r]++;
      // Get indices
      casadi_int id_c = jac_in_[c], id_r = jac_in_[r];
      // Nonzero
      double nz = hess_nz[k], nz_tr = hess_nz[k_tr];
      // Check if entry is NaN of inf
      if (std::isnan(nz) || std::isinf(nz)) {
        std::stringstream ss;
        ss << "Second derivative w.r.t. " << get_fmu(fmu_)->desc_in(m, id_r) << " and "
          << get_fmu(fmu_)->desc_in(m, id_r) << " is " << nz;
        casadi_warning(ss.str());
        // Further checks not needed for entry
        continue;
      }
      // Symmetry check (strict upper triangular part only)
      if (r < c) {
        // Normaliation factor to be used for relative tolerance
        double nz_max = std::fmax(std::fabs(nz), std::fabs(nz_tr));
        // Check if above absolute and relative tolerance bounds
        if (nz_max > abstol_ && std::fabs(nz - nz_tr) > nz_max * reltol_) {
          std::stringstream ss;
          ss << "Hessian appears nonsymmetric. Got " << nz << " vs. " << nz_tr
            << " for second derivative w.r.t. " << get_fmu(fmu_)->desc_in(m, id_r) << " and "
            << get_fmu(fmu_)->desc_in(m, id_c) << ", hess_nz = " << k << "/" <<  k_tr;
          casadi_warning(ss.str());
        }
      }
    }
  }
}

void FmuFunction::remove_nans(double *hess_nz, casadi_int* iw) const {
  // Get Hessian sparsity pattern
  casadi_int n = hess_sp_.size1();
  const casadi_int *colind = hess_sp_.colind(), *row = hess_sp_.row();
  // Mark variables that have been perturbed
  std::fill(iw, iw + n, 0);
  casadi_int hess_colors_nnz = hess_colors_.nnz();
  const casadi_int* hess_colors_row = hess_colors_.row();
  for (casadi_int k = 0; k < hess_colors_nnz; ++k)
    iw[hess_colors_row[k]] = 1;
  // Nonzero counters for transpose
  casadi_copy(colind, n, iw);
  // Loop over Hessian columns
  for (casadi_int c = 0; c < n; ++c) {
    // Loop over nonzeros for the column
    for (casadi_int k = colind[c]; k < colind[c + 1]; ++k) {
      // Get row of Hessian
      casadi_int r = row[k];
      // Get nonzero of transpose
      casadi_int k_tr = iw[r]++;
      // If NaN, use element from transpose
      if (std::isnan(hess_nz[k])) {
        hess_nz[k] = hess_nz[k_tr];
      }
    }
  }
}

void FmuFunction::make_symmetric(double *hess_nz, casadi_int* iw) const {
  // Get Hessian sparsity pattern
  casadi_int n = hess_sp_.size1();
  const casadi_int *colind = hess_sp_.colind(), *row = hess_sp_.row();
  // Nonzero counters for transpose
  casadi_copy(colind, n, iw);
  // Loop over Hessian columns
  for (casadi_int c = 0; c < n; ++c) {
    // Loop over nonzeros for the column
    for (casadi_int k = colind[c]; k < colind[c + 1]; ++k) {
      // Get row of Hessian
      casadi_int r = row[k];
      // Get nonzero of transpose
      casadi_int k_tr = iw[r]++;
      // Make symmetric
      if (r < c) hess_nz[k] = hess_nz[k_tr] = 0.5 * (hess_nz[k] + hess_nz[k_tr]);
    }
  }
}

std::string to_string(Parallelization v) {
  switch (v) {
  case Parallelization::SERIAL: return "serial";
  case Parallelization::OPENMP: return "openmp";
  case Parallelization::THREAD: return "thread";
  default: break;
  }
  return "";
}

bool has_prefix(const std::string& s) {
  return s.find('_') < s.size();
}

std::string pop_prefix(const std::string& s, std::string* rem) {
  // Get prefix
  casadi_assert_dev(!s.empty());
  size_t pos = s.find('_');
  casadi_assert(pos < s.size(), "Cannot process \"" + s + "\"");
  // Get prefix
  std::string r = s.substr(0, pos);
  // Remainder, if requested (note that rem == &s is possible)
  if (rem) *rem = s.substr(pos+1, std::string::npos);
  // Return prefix
  return r;
}

bool FmuFunction::all_regular() const {
  // Look for any non-regular input
  for (auto&& e : in_) if (e.type != InputType::REG) return false;
  // Look for any non-regular output
  for (auto&& e : out_) if (e.type != OutputType::REG) return false;
  // Only regular inputs and outputs
  return true;
}

bool FmuFunction::all_vectors() const {
  // Check inputs
  for (auto&& e : in_) {
    switch (e.type) {
      // Supported for derivative calculations
      case InputType::REG:
      case InputType::ADJ:
      case InputType::OUT:
        break;
      // Not supported
      default:
        return false;
    }
  }
  // Check outputs
  for (auto&& e : out_) {
    // Supported for derivative calculations
    switch (e.type) {
      case OutputType::REG:
      case OutputType::ADJ:
        break;
      // Not supported
      default:
        return false;
    }
  }
  // OK if reached this point
  return true;
}

Function FmuFunction::factory(const std::string& name,
    const std::vector<std::string>& s_in,
    const std::vector<std::string>& s_out,
    const Function::AuxOut& aux,
    const Dict& opts) const {
  // Assume we can call constructor directly
  try {
    // Hack: Inherit parallelization, verbosity option
    Dict opts1 = opts;
    opts1["parallelization"] = to_string(parallelization_);
    opts1["verbose"] = verbose_;
    opts1["print_progress"] = print_progress_;
    // Replace ':' with '_' in s_in and s_out
    std::vector<std::string> s_in_mod = s_in, s_out_mod = s_out;
    for (std::string& s : s_in_mod) std::replace(s.begin(), s.end(), ':', '_');
    for (std::string& s : s_out_mod) std::replace(s.begin(), s.end(), ':', '_');
    // New instance of the same class, using the same Fmu instance
    Function ret;
    ret.own(new FmuFunction(name, fmu_, s_in_mod, s_out_mod));
    ret->construct(opts1);
    return ret;
  } catch (std::exception& e) {
    casadi_warning("FmuFunction::factory call for constructing " + name + " from " + name_
      + " failed:\n" + std::string(e.what()) + "\nFalling back to base class implementation");
  }
  // Fall back to base class
  return FunctionInternal::factory(name, s_in, s_out, aux, opts);
}

bool FmuFunction::has_jacobian() const {
  // Calculation of Hessian inside FmuFunction (in development)
  if (new_jacobian_ && all_vectors()) return true;
  // Only first order
  return all_regular();
}

Function FmuFunction::get_jacobian(const std::string& name, const std::vector<std::string>& inames,
    const std::vector<std::string>& onames, const Dict& opts) const {
  // Hack: Inherit parallelization, verbosity option
  Dict opts1 = opts;
  opts1["parallelization"] = to_string(parallelization_);
  opts1["verbose"] = verbose_;
  opts1["print_progress"] = print_progress_;
  // Return new instance of class
  Function ret;
  ret.own(new FmuFunction(name, fmu_, inames, onames));
  ret->construct(opts1);
  return ret;
}

bool FmuFunction::has_reverse(casadi_int nadj) const {
  // Only first order analytic derivative possible
  if (!all_regular()) return false;
  // Otherwise: Only 1 direction implemented
  return nadj == 1;
}

Function FmuFunction::get_reverse(casadi_int nadj, const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const {
  // Only single directional derivative implemented
  casadi_assert(nadj == 1, "Not implemented");
  // Hack: Inherit parallelization option
  Dict opts1 = opts;
  opts1["parallelization"] = to_string(parallelization_);
  opts1["verbose"] = verbose_;
  opts1["new_jacobian"] = new_hessian_;
  opts1["print_progress"] = print_progress_;
  // Return new instance of class
  Function ret;
  ret.own(new FmuFunction(name, fmu_, inames, onames));
  ret->construct(opts1);
  return ret;
}

bool FmuFunction::has_jac_sparsity(casadi_int oind, casadi_int iind) const {
  // Available in the FMU meta information
  if (out_.at(oind).type == OutputType::REG) {
    if (in_.at(iind).type == InputType::REG) {
      return true;
    } else if (in_.at(iind).type == InputType::ADJ) {
      return true;
    }
  } else if (out_.at(oind).type == OutputType::ADJ) {
    if (in_.at(iind).type == InputType::REG) {
      return true;
    } else if (in_.at(iind).type == InputType::ADJ) {
      return true;
    }
  }
  // Not available
  return false;
}

Sparsity FmuFunction::get_jac_sparsity(casadi_int oind, casadi_int iind,
    bool symmetric) const {
  // Available in the FMU meta information
  if (out_.at(oind).type == OutputType::REG) {
    if (in_.at(iind).type == InputType::REG) {
      return get_fmu(fmu_)->jac_sparsity(out_.at(oind).ind, in_.at(iind).ind);
    } else if (in_.at(iind).type == InputType::ADJ) {
      return Sparsity(nnz_out(oind), nnz_in(iind));
    }
  } else if (out_.at(oind).type == OutputType::ADJ) {
    if (in_.at(iind).type == InputType::REG) {
      return get_fmu(fmu_)->hess_sparsity(out_.at(oind).wrt, in_.at(iind).ind);
    } else if (in_.at(iind).type == InputType::ADJ) {
      return get_fmu(fmu_)->jac_sparsity(in_.at(iind).ind, out_.at(oind).wrt).T();
    }
  }
  // Not available
  casadi_error("Implementation error");
  return Sparsity();
}

Dict FmuFunction::get_stats(void *mem) const {
  // Get the stats from the base classes
  Dict stats = FunctionInternal::get_stats(mem);
  // Get memory object
  FmuMemory* m = static_cast<FmuMemory*>(mem);
  // Get auxilliary variables from Fmu
  get_fmu(fmu_)->get_stats(m, &stats, name_in_, get_ptr(in_));
  // Return stats
  return stats;
}

Fmu2* FmuFunction::get_fmu(const Fmu& fmu) {
  return static_cast<Fmu2*>(fmu.get());
}

#endif  // WITH_FMU

} // namespace casadi
