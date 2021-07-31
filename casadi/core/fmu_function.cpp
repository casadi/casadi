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


#include "fmu_function.hpp"
#include "casadi_misc.hpp"
#include "serializing_stream.hpp"
#include "dae_builder_internal.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

namespace casadi {

#ifdef WITH_FMU

FmuFunction::FmuFunction(const std::string& name, const DaeBuilder& dae,
    const std::vector<std::vector<size_t>>& id_in,
    const std::vector<std::vector<size_t>>& id_out,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out)
  : FunctionInternal(name), dae_(dae), id_in_(id_in), id_out_(id_out) {
  // Names of inputs
  if (!name_in.empty()) {
    casadi_assert(id_in.size() == name_in.size(), "Mismatching number of input names");
    name_in_ = name_in;
  }
  // Names of outputs
  if (!name_out.empty()) {
    casadi_assert(id_out.size() == name_out.size(), "Mismatching number of output names");
    name_out_ = name_out;
  }
}

FmuFunction::~FmuFunction() {
  // Free memory
  clear_mem();
}

void FmuFunction::init(const Dict& opts) {
  // Call the initialization method of the base class
  FunctionInternal::init(opts);

  // Get a pointer to the DaeBuilder class
  casadi_assert(dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(dae_->raw_);

  // Load on first encounter
  if (dae->fmu_ == 0) dae->init_fmu();
}

Sparsity FmuFunction::get_sparsity_in(casadi_int i) {
  return Sparsity::dense(id_in_.at(i).size(), 1);
}

Sparsity FmuFunction::get_sparsity_out(casadi_int i) {
  return Sparsity::dense(id_out_.at(i).size(), 1);
}

int FmuFunction::eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // DaeBuilder instance
  casadi_assert(dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(dae_->raw_);
  // Create instance
  int m = dae->fmu_->checkout();
  // Evaluate fmu
  int flag = eval_fmu(dae, m, arg, res);
  // Release memory object
  dae->fmu_->release(m);
  // Return error flag
  return flag;
}

int FmuFunction::eval_fmu(const DaeBuilderInternal* dae, int mem,
    const double** arg, double** res) const {
  // Set inputs
  for (size_t k = 0; k < id_in_.size(); ++k) {
    for (size_t i = 0; i < id_in_[k].size(); ++i) {
      dae->fmu_->set(mem, id_in_[k][i], arg[k] ? arg[k][i] : 0);
    }
  }
  // Request outputs to be evaluated
  for (size_t k = 0; k < id_out_.size(); ++k) {
    if (res[k]) {
      for (size_t i = 0; i < id_out_[k].size(); ++i) {
        dae->fmu_->request(mem, id_out_[k][i]);
      }
    }
  }
  // Reset solver
  if (dae->fmu_->setup_experiment(mem)) return 1;
  // Initialization mode begins
  if (dae->fmu_->enter_initialization_mode(mem)) return 1;
  // Initialization mode ends
  if (dae->fmu_->exit_initialization_mode(mem)) return 1;
  // Evaluate
  if (dae->fmu_->eval(mem)) return 1;
  // Reset solver
  if (dae->fmu_->reset(mem)) return 1;
  // Get outputs
  for (size_t k = 0; k < id_out_.size(); ++k) {
    if (res[k]) {
      for (size_t i = 0; i < id_out_[k].size(); ++i) {
        dae->fmu_->get(mem, id_out_[k][i], &res[k][i]);
      }
    }
  }
  // Successful return
  return 0;
}

int FmuFunction::eval_jac(const DaeBuilderInternal* dae, int mem,
    const double** arg, double** res) const {
  // Set inputs
  for (size_t k = 0; k < id_in_.size(); ++k) {
    for (size_t i = 0; i < id_in_[k].size(); ++i) {
      dae->fmu_->set(mem, id_in_[k][i], arg[k] ? arg[k][i] : 0);
    }
  }
  // Reset solver
  if (dae->fmu_->setup_experiment(mem)) return 1;
  // Initialization mode begins
  if (dae->fmu_->enter_initialization_mode(mem)) return 1;
  // Initialization mode ends
  if (dae->fmu_->exit_initialization_mode(mem)) return 1;
  // Evaluate
  if (dae->fmu_->eval(mem)) return 1;
  // Loop over function inputs
  for (size_t i1 = 0; i1 < id_in_.size(); ++i1) {
    // Calculate Jacobian, one column at a time
    for (casadi_int i2 = 0; i2 < id_in_[i1].size(); ++i2) {
      // Set seed for column
      dae->fmu_->set_seed(mem, id_in_[i1][i2], 1.);
      // Request all elements of the column
      for (size_t j1 = 0; j1 < id_out_.size(); ++j1) {
        if (res[j1 * id_in_.size() + i1]) {
          for (size_t j2 = 0; j2 < id_out_[j1].size(); ++j2) {
            dae->fmu_->request(mem, id_out_[j1][j2]);
          }
        }
      }
      // Calculate derivatives
      if (dae->fmu_->eval_derivative(mem)) return 1;
      // Loop over function outputs
      for (size_t j1 = 0; j1 < id_out_.size(); ++j1) {
        // Corresponding Jacobian block, if any
        double* J = res[j1 * id_in_.size() + i1];
        if (J) {
          // Shift to right column
          J += id_out_[j1].size() * i2;
          // Get column
          for (size_t j2 = 0; j2 < id_out_[j1].size(); ++j2) {
            dae->fmu_->get_sens(mem, id_out_[j1][j2], J++);
          }
        }
      }
    }
  }
  // Reset solver
  if (dae->fmu_->reset(mem)) return 1;
  // Successful return
  return 0;
}

int FmuFunction::eval_adj(const DaeBuilderInternal* dae, int mem,
    const double** arg, double** res) const {
  // Set inputs
  for (size_t k = 0; k < id_in_.size(); ++k) {
    for (size_t i = 0; i < id_in_[k].size(); ++i) {
      dae->fmu_->set(mem, id_in_[k][i], arg[k] ? arg[k][i] : 0);
    }
  }
  // Setup experiment
  if (dae->fmu_->setup_experiment(mem)) return 1;
  // Initialization mode begins
  if (dae->fmu_->enter_initialization_mode(mem)) return 1;
  // Initialization mode ends
  if (dae->fmu_->exit_initialization_mode(mem)) return 1;
  // Evaluate
  if (dae->fmu_->eval(mem)) return 1;
  // Loop over function inputs
  for (size_t i1 = 0; i1 < id_in_.size(); ++i1) {
    // Sensitivities to be calculated
    double* sens = res[i1];
    // Skip if not requested
    if (sens == 0) continue;
    // Calculate Jacobian, one column at a time
    for (casadi_int i2 = 0; i2 < id_in_[i1].size(); ++i2) {
      // Initialize return to zero
      sens[i2] = 0;
      // Set seed for column i
      dae->fmu_->set_seed(mem, id_in_[i1][i2], 1.);
      // Request all elements of the column, unless corresponding seed is zero
      for (size_t j1 = 0; j1 < id_out_.size(); ++j1) {
        if (arg[id_in_.size() + id_out_.size() + j1]) {
          for (size_t j2 = 0; j2 < id_out_[j1].size(); ++j2) {
            dae->fmu_->request(mem, id_out_[j1][j2]);
          }
        }
      }
      // Calculate derivatives
      if (dae->fmu_->eval_derivative(mem)) return 1;
      // Get sensitivities
      for (size_t j1 = 0; j1 < id_out_.size(); ++j1) {
        const double* seed = arg[id_in_.size() + id_out_.size() + j1];
        if (seed) {
          for (size_t j2 = 0; j2 < id_out_[j1].size(); ++j2) {
            double J_ij;
            dae->fmu_->get_sens(mem, id_out_[j1][j2], &J_ij);
            sens[i2] += seed[j2] * J_ij;
          }
        }
      }
    }
  }
  // Reset solver
  if (dae->fmu_->reset(mem)) return 1;
  // Successful return
  return 0;
}

bool FmuFunction::has_jacobian() const {
  casadi_assert(dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(dae_->raw_);
  return dae->provides_directional_derivative_;
}

Function FmuFunction::get_jacobian(const std::string& name, const std::vector<std::string>& inames,
    const std::vector<std::string>& onames, const Dict& opts) const {
  Function ret;
  ret.own(new FmuFunctionJac(name));
  // Hack: Manually enable finite differenting (pending implementation in class)
  Dict opts2 = opts;
  opts2["enable_fd"] = true;
  ret->construct(opts2);
  return ret;
}

bool FmuFunction::has_reverse(casadi_int nadj) const {
  casadi_assert(dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(dae_->raw_);
  return dae->provides_directional_derivative_ && nadj == 1;
}

Function FmuFunction::get_reverse(casadi_int nadj, const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const {
  casadi_assert(nadj == 1, "Not supported");
  Function ret;
  ret.own(new FmuFunctionAdj(name));
  // Hack: Manually enable finite differenting (pending implementation in class)
  Dict opts2 = opts;
  opts2["enable_fd"] = true;
  ret->construct(opts2);
  return ret;
}

FmuFunctionJac::~FmuFunctionJac() {
  // Free memory
  clear_mem();
}

int FmuFunctionJac::eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // Non-differentiated class
  auto self = derivative_of_.get<FmuFunction>();
  // DaeBuilder instance
  casadi_assert(self->dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(self->dae_->raw_);
  // Create instance
  int m = dae->fmu_->checkout();
  // Evaluate fmu
  int flag = self->eval_jac(dae, m, arg, res);
  // Release memory object
  dae->fmu_->release(m);
  // Return error flag
  return flag;
}

FmuFunctionAdj::~FmuFunctionAdj() {
  // Free memory
  clear_mem();
}

int FmuFunctionAdj::eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // Non-differentiated class
  auto self = derivative_of_.get<FmuFunction>();
  // DaeBuilder instance
  casadi_assert(self->dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(self->dae_->raw_);
  // Create instance
  int m = dae->fmu_->checkout();
  // Evaluate fmu
  int flag = self->eval_adj(dae, m, arg, res);
  // Release memory object
  dae->fmu_->release(m);
  // Return error flag
  return flag;
}

#endif  // WITH_FMU

} // namespace casadi
