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


#include "qpchasm.hpp"
#include "casadi/core/nlpsol.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_QPCHASM_EXPORT
  casadi_register_conic_qpchasm(Conic::Plugin* plugin) {
    plugin->creator = Qpchasm::creator;
    plugin->name = "qpchasm";
    plugin->doc = Qpchasm::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Qpchasm::options_;
    plugin->deserialize = &Qpchasm::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_QPCHASM_EXPORT casadi_load_conic_qpchasm() {
    Conic::registerPlugin(casadi_register_conic_qpchasm);
  }

  void print_vec(const std::string& str, const double* v, casadi_int n) {
    uout() << str << " =";
    for (casadi_int k = 0; k < n; ++k) uout() << " " << v[k];
    uout() << "\n";
  }

  Qpchasm::Qpchasm(const std::string& name, const std::map<std::string, Sparsity> &st)
    : Conic(name, st) {
  }

  Qpchasm::~Qpchasm() {
    clear_mem();
  }

  const Options Qpchasm::options_
  = {{&Conic::options_},
     {{"max_iter",
       {OT_INT,
        "Maximum number of iterations [1000]."}},
      {"constr_viol_tol",
       {OT_DOUBLE,
        "Constraint violation tolerance [1e-8]."}},
      {"dual_inf_tol",
       {OT_DOUBLE,
        "Dual feasibility violation tolerance [1e-8]"}},
      {"print_header",
       {OT_BOOL,
        "Print header [true]."}},
      {"print_iter",
       {OT_BOOL,
        "Print iterations [true]."}},
      {"print_info",
       {OT_BOOL,
        "Print info [true]."}},
      {"linear_solver",
       {OT_STRING,
        "A custom linear solver creator function [default: qr]"}},
      {"linear_solver_options",
       {OT_DICT,
        "Options to be passed to the linear solver"}},
      {"min_lam",
       {OT_DOUBLE,
        "Smallest multiplier treated as inactive for the initial active set [0]."}}
     }
  };

  void Qpchasm::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    // Assemble KKT system sparsity
    kkt_ = Sparsity::kkt(H_, A_, true, true);

    // Setup memory structure
    set_qp_prob();

    // Default options
    print_iter_ = true;
    print_header_ = true;
    print_info_ = true;
    linear_solver_ = "qr";

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="max_iter") {
        p_.max_iter = op.second;
      } else if (op.first=="constr_viol_tol") {
        p_.constr_viol_tol = op.second;
      } else if (op.first=="dual_inf_tol") {
        p_.dual_inf_tol = op.second;
      } else if (op.first=="min_lam") {
        p_.min_lam = op.second;
      } else if (op.first=="print_iter") {
        print_iter_ = op.second;
      } else if (op.first=="print_header") {
        print_header_ = op.second;
      } else if (op.first=="print_info") {
        print_info_ = op.second;
      } else if (op.first=="linear_solver") {
        linear_solver_ = op.second.to_string();
      } else if (op.first=="linear_solver_options") {
        linear_solver_options_ = op.second;
      }
    }

    // Allocate memory for IP solver
    alloc_w(casadi_ipqp_sz_w(&p_), true);

    // For KKT formation
    alloc_w(kkt_.nnz(), true);
    alloc_iw(A_.size2());
    alloc_w(nx_ + na_);

    // KKT solver
    linsol_ = Linsol("linsol", linear_solver_, kkt_, linear_solver_options_);

    if (print_header_) {
      // Print summary
      print("-------------------------------------------\n");
      print("This is casadi::QPCHASM\n");
      print("Number of variables:                       %9d\n", nx_);
      print("Number of constraints:                     %9d\n", na_);
      print("Number of nonzeros in H:                   %9d\n", H_.nnz());
      print("Number of nonzeros in A:                   %9d\n", A_.nnz());
      print("Number of nonzeros in KKT:                 %9d\n", kkt_.nnz());
    }
  }

  void Qpchasm::set_qp_prob() {
    casadi_ipqp_setup(&p_, nx_, na_);
  }

  int Qpchasm::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<QpchasmMemory*>(mem);
    m->return_status = "";
    return 0;
  }

  int Qpchasm::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<QpchasmMemory*>(mem);
    // Message buffer
    char buf[121];
    // Setup data structure
    casadi_ipqp_data<double> d;
    d.prob = &p_;
    d.g = arg[CONIC_G];
    double* nz_kkt = w; w += kkt_.nnz();
    casadi_ipqp_init(&d, &iw, &w);
    // Pass bounds on z
    casadi_copy(arg[CONIC_LBX], nx_, d.lbz);
    casadi_copy(arg[CONIC_LBA], na_, d.lbz+nx_);
    casadi_copy(arg[CONIC_UBX], nx_, d.ubz);
    casadi_copy(arg[CONIC_UBA], na_, d.ubz+nx_);
    // Pass initial guess
    casadi_copy(arg[CONIC_X0], nx_, d.z);
    casadi_fill(d.z+nx_, na_, nan);
    casadi_copy(arg[CONIC_LAM_X0], nx_, d.lam);
    casadi_copy(arg[CONIC_LAM_A0], na_, d.lam+nx_);
    casadi_fill(d.lam_lbz, p_.nz, 0.);
    casadi_fill(d.lam_ubz, p_.nz, 0.);
    // Checkout a linear solver instance
    int linsol_mem = linsol_.checkout();
    // New QP
    d.next = IPQP_RESET;
    // Reverse communication loop
    while (casadi_ipqp(&d)) {
      switch (d.task) {
      case IPQP_MV:
        // Matrix-vector multiplication
        casadi_clear(d.rz, p_.nz);
        casadi_mv(arg[CONIC_H], H_, d.z, d.rz, 0);
        casadi_mv(arg[CONIC_A], A_, d.lam + p_.nx, d.rz, 1);
        casadi_mv(arg[CONIC_A], A_, d.z, d.rz + p_.nx, 0);
        break;
      case IPQP_PROGRESS:
        // Print progress
        if (print_iter_) {
          if (d.iter % 10 == 0) {
            // Print header
            if (casadi_ipqp_print_header(&d, buf, sizeof(buf))) break;
            uout() << buf << "\n";
          }
          // Print iteration
          if (casadi_ipqp_print_iteration(&d, buf, sizeof(buf))) break;
          uout() << buf << "\n";
          // User interrupt?
          InterruptHandler::check();
        }
        break;
      case IPQP_FACTOR:
        // Form KKT
        casadi_kkt(kkt_, nz_kkt, H_, arg[CONIC_H], A_, arg[CONIC_A],
          d.S, d.D, w, iw);
        // Factorize KKT
        if (linsol_.nfact(nz_kkt, linsol_mem))
          d.status = IPQP_FACTOR_ERROR;
        break;
      case IPQP_SOLVE:
        // Solve KKT
        if (linsol_.solve(nz_kkt, d.linsys, 1, false, linsol_mem))
          d.status = IPQP_SOLVE_ERROR;
        break;
      }
    }
    // Check return flag
    switch (d.status) {
      case IPQP_SUCCESS:
        m->return_status = "success";
        break;
      case IPQP_MAX_ITER:
        m->return_status = "Maximum number of iterations reached";
        m->unified_return_status = SOLVER_RET_LIMITED;
        break;
      case IPQP_NO_SEARCH_DIR:
        m->return_status = "Failed to calculate search direction";
        break;
      case IPQP_MV_ERROR:
        m->return_status = "Matrix-vector evaluation error";
        break;
      case IPQP_FACTOR_ERROR:
        m->return_status = "Linear solver factorization error";
        break;
      case IPQP_SOLVE_ERROR:
        m->return_status = "Linear solver solution error";
        break;
      case IPQP_PROGRESS_ERROR:
        m->return_status = "Printing error";
        break;
    }
    // Release linear solver instance
    linsol_.release(linsol_mem);
    // Get solution
    casadi_copy(d.z, nx_, res[CONIC_X]);
    casadi_copy(d.lam, nx_, res[CONIC_LAM_X]);
    casadi_copy(d.lam+nx_, na_, res[CONIC_LAM_A]);
    // Calculate optimal cost
    if (res[CONIC_COST]) {
      *res[CONIC_COST] = .5 * casadi_bilin(arg[CONIC_H], H_, d.z, d.z)
        + casadi_dot(p_.nx, d.z, d.g);
    }
    // Return
    if (verbose_) casadi_warning(m->return_status);
    m->success = d.status == IPQP_SUCCESS;
    return 0;
  }

  Dict Qpchasm::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<QpchasmMemory*>(mem);
    stats["return_status"] = m->return_status;
    return stats;
  }

  Qpchasm::Qpchasm(DeserializingStream& s) : Conic(s) {
    s.version("Qpchasm", 1);
    s.unpack("Qpchasm::kkt", kkt_);
    s.unpack("Qpchasm::print_iter", print_iter_);
    s.unpack("Qpchasm::print_header", print_header_);
    s.unpack("Qpchasm::print_info", print_info_);
    set_qp_prob();
    s.unpack("Qpchasm::max_iter", p_.max_iter);
    s.unpack("Qpchasm::min_lam", p_.min_lam);
    s.unpack("Qpchasm::constr_viol_tol", p_.constr_viol_tol);
    s.unpack("Qpchasm::dual_inf_tol", p_.dual_inf_tol);
  }

  void Qpchasm::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);

    s.version("Qpchasm", 1);
    s.pack("Qpchasm::kkt", kkt_);
    s.pack("Qpchasm::print_iter", print_iter_);
    s.pack("Qpchasm::print_header", print_header_);
    s.pack("Qpchasm::print_info", print_info_);
    s.pack("Qpchasm::max_iter", p_.max_iter);
    s.pack("Qpchasm::min_lam", p_.min_lam);
    s.pack("Qpchasm::constr_viol_tol", p_.constr_viol_tol);
    s.pack("Qpchasm::dual_inf_tol", p_.dual_inf_tol);
  }

} // namespace casadi
