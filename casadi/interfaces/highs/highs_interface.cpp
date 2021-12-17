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

#include "highs_interface.hpp"

#include "casadi/core/nlp_tools.hpp"

namespace casadi {

  using namespace std;

  extern "C"
  int CASADI_CONIC_HIGHS_EXPORT
  casadi_register_conic_highs(Conic::Plugin* plugin) {
    plugin->creator = HighsInterface::creator;
    plugin->name = "highs";
    plugin->doc = HighsInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &HighsInterface::options_;
    plugin->deserialize = &HighsInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_HIGHS_EXPORT casadi_load_conic_highs() {
    Conic::registerPlugin(casadi_register_conic_highs);
  }


  HighsInterface::HighsInterface(const std::string& name,
                             const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
  }

  const Options HighsInterface::options_
  = {{&Conic::options_},
     {{"highs",
       {OT_DICT,
        "Options to be passed to HiGHS."
        }},
     }
   };

  // Options with description https://www.maths.ed.ac.uk/hall/HiGHS/HighsOptions.html
  std::list<std::string> HighsInterface::param_bool = {
    "output_flag",
    "log_to_console",
    "write_solution_to_file",
    "mip_detect_symmetry"
  };

  void HighsInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Conic::init(opts);

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="highs") {
        opts_ = op.second;
      }
    }

    // Allocate work vectors
    alloc_w(nx_, true); // g
    alloc_w(nx_, true); // lbx
    alloc_w(nx_, true); // ubx
    alloc_w(na_, true); // lba
    alloc_w(na_, true); // uba
    alloc_w(nnz_in(CONIC_H), true); // H
    alloc_w(nnz_in(CONIC_A), true); // A
  }

  int HighsInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    if (!mem) return 1;
    auto m = static_cast<HighsMemory*>(mem);

    m->add_stat("preprocessing");
    m->add_stat("solver");
    m->add_stat("postprocessing");

    m->colinda.resize(A_.size2()+1);
    m->rowa.resize(A_.nnz());
    m->colindh.resize(H_.size2()+1);
    m->rowh.resize(H_.nnz());
    m->integrality.resize(nx_);

    return 0;
  }

  int HighsInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<HighsMemory*>(mem);

    // Problem has not been solved at this point
    m->return_status = static_cast<int>(HighsStatus::kError);

    m->fstats.at("preprocessing").tic();

    // Get inputs
    double* g=w; w += nx_;
    casadi_copy(arg[CONIC_G], nx_, g);
    double* lbx=w; w += nx_;
    casadi_copy(arg[CONIC_LBX], nx_, lbx);
    double* ubx=w; w += nx_;
    casadi_copy(arg[CONIC_UBX], nx_, ubx);
    double* lba=w; w += na_;
    casadi_copy(arg[CONIC_LBA], na_, lba);
    double* uba=w; w += na_;
    casadi_copy(arg[CONIC_UBA], na_, uba);
    double* H=w; w += nnz_in(CONIC_H);
    casadi_copy(arg[CONIC_H], nnz_in(CONIC_H), H);
    double* A=w; w += nnz_in(CONIC_A);
    casadi_copy(arg[CONIC_A], nnz_in(CONIC_A), A);

    copy_vector(A_.colind(), m->colinda);
    copy_vector(A_.row(), m->rowa);
    copy_vector(H_.colind(), m->colindh);
    copy_vector(H_.row(), m->rowh);


    // Create HiGHS instance and pass problem
    Highs highs;
    HighsStatus status;
    const int matrix_format = 1;
    const int sense = 1;
    const double offset = 0.0;

    // set HiGHS options, option verification is done by HiGHS
    for (auto&& op : opts_) {
      auto it = std::find(param_bool.begin(), param_bool.end(), op.first);
      if(it != param_bool.end() || op.second.getType() == OT_BOOL) {
        casadi_assert(highs.setOptionValue(op.first, op.second.to_bool()) == HighsStatus::kOk,
          "Error setting option '" + op.first + "'.");
      } 
      else if(op.second.getType() == OT_INT){
        casadi_assert(highs.setOptionValue(op.first, static_cast<int>(op.second.to_int())) == HighsStatus::kOk,
          "Error setting option '" + op.first + "'.");
      }
      else if(op.second.getType() == OT_DOUBLE) {
        casadi_assert(highs.setOptionValue(op.first, op.second.to_double()) == HighsStatus::kOk,
          "Error setting option '" + op.first + "'.");
      }
      else if(op.second.getType() == OT_STRING) {
        casadi_assert(highs.setOptionValue(op.first, op.second.to_string()) == HighsStatus::kOk,
          "Error setting option '" + op.first + "'.");
      }
      else {
        casadi_assert(false, "Option type for '" + op.first + "'not supported!");
      }
    }

    // if variables are declared as discrete, set integrality pointer and flag discrete variables
    int* integrality_ptr = nullptr;
    if (!discrete_.empty()) {
      integrality_ptr = get_ptr(m->integrality);
      for (casadi_int i=0; i<nx_; ++i) {
        m->integrality[i] = discrete_.at(i) ? 1 : 0;
      }
    }
    status = highs.passModel(nx_, na_, nnz_in(CONIC_A), nnz_in(CONIC_H),
      matrix_format, matrix_format, sense, offset,
      g, lbx, ubx, lba, uba,
      get_ptr(m->colinda), get_ptr(m->rowa), A,
      get_ptr(m->colindh), get_ptr(m->rowh), H,      
      integrality_ptr);
    
    // check that passing model is successful
    casadi_assert(status == HighsStatus::kOk, "invalid data to build HiGHS model");

    m->fstats.at("preprocessing").toc();
    m->fstats.at("solver").tic();

    // solve incumbent model
    status = highs.run();
    casadi_assert(status == HighsStatus::kOk, "running HiGHS failed");
    
    m->fstats.at("solver").toc();
    m->fstats.at("postprocessing").tic();

    // get primal and dual solution
    HighsSolution solution = highs.getSolution();
    casadi_copy(solution.col_value.data(), nx_, res[CONIC_X]);
    
    if (res[CONIC_LAM_X]) {
      casadi_copy(solution.col_dual.data(), nx_, res[CONIC_LAM_X]);
      casadi_scal(nx_, -1., res[CONIC_LAM_X]);
    }
    if (res[CONIC_LAM_A]) {
      casadi_copy(solution.row_dual.data(), na_, res[CONIC_LAM_A]);
      casadi_scal(na_, -1., res[CONIC_LAM_A]);
    }
    // get information
    const HighsInfo& info = highs.getInfo();

    if (res[CONIC_COST]) {
      *res[CONIC_COST] = info.objective_function_value;
    }

    m->fstats.at("postprocessing").toc();

    HighsModelStatus model_status = highs.getModelStatus();
    m->return_status = static_cast<int>(model_status);
    m->success = model_status==HighsModelStatus::kOptimal;

    if (model_status==HighsModelStatus::kTimeLimit
          || model_status==HighsModelStatus::kIterationLimit)
      m->unified_return_status = SOLVER_RET_LIMITED;

    m->simplex_iteration_count = info.simplex_iteration_count;
    m->simplex_iteration_count = info.simplex_iteration_count;
    m->ipm_iteration_count = info.ipm_iteration_count;
    m->qp_iteration_count = info.qp_iteration_count;
    m->crossover_iteration_count = info.crossover_iteration_count;
    m->primal_solution_status  = info.primal_solution_status;
    m->dual_solution_status = info.dual_solution_status;
    m->basis_validity = info.basis_validity;
    m->mip_dual_bound = info.mip_dual_bound;
    m->mip_gap = info.mip_gap;
    m->num_primal_infeasibilities = info.num_primal_infeasibilities;
    m->max_primal_infeasibility = info.max_primal_infeasibility;
    m->sum_primal_infeasibilities = info.sum_primal_infeasibilities;
    m->num_dual_infeasibilities = info.num_dual_infeasibilities;
    m->max_dual_infeasibility = info.max_dual_infeasibility;
    m->sum_dual_infeasibilities = info.sum_dual_infeasibilities;

    return 0;
  }

  HighsInterface::~HighsInterface() {
    clear_mem();
  }

  HighsMemory::HighsMemory() {
  }

  HighsMemory::~HighsMemory() {
  }


  Dict HighsInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    Highs highs;
    auto m = static_cast<HighsMemory*>(mem);
    stats["return_status"] = highs.modelStatusToString(static_cast<HighsModelStatus>(m->return_status));
    stats["simplex_iteration_count"] = m->simplex_iteration_count;
    stats["simplex_iteration_count"] = m->simplex_iteration_count;
    stats["ipm_iteration_count"] = m->ipm_iteration_count;
    stats["qp_iteration_count"] = m->qp_iteration_count;
    stats["crossover_iteration_count"] = m->crossover_iteration_count;
    stats["primal_solution_status"]  = highs.solutionStatusToString(m->primal_solution_status);
    stats["dual_solution_status"] = highs.solutionStatusToString(m->dual_solution_status);
    stats["basis_validity"] = highs.basisValidityToString(m->basis_validity);
    stats["mip_dual_bound"] = m->mip_dual_bound;
    stats["mip_gap"] = m->mip_gap;
    stats["num_primal_infeasibilities"] = m->num_primal_infeasibilities;
    stats["max_primal_infeasibility"] = m->max_primal_infeasibility;
    stats["sum_primal_infeasibilities"] = m->sum_primal_infeasibilities;
    stats["num_dual_infeasibilities"] = m->num_dual_infeasibilities;
    stats["max_dual_infeasibility"] = m->max_dual_infeasibility;
    stats["sum_dual_infeasibilities"] = m->sum_dual_infeasibilities;
    return stats;
  }

  HighsInterface::HighsInterface(DeserializingStream& s) : Conic(s) {
    s.version("HighsInterface", 1);
    s.unpack("HighsInterface::opts", opts_);
  }

  void HighsInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);

    s.version("HighsInterface", 1);
    s.pack("HighsInterface::opts", opts_);
  }

} // end namespace casadi
