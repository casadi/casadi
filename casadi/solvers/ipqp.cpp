/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
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


#include "ipqp.hpp"
#include "casadi/core/nlpsol.hpp"

namespace casadi {

  extern "C"
  int CASADI_CONIC_IPQP_EXPORT
  casadi_register_conic_ipqp(Conic::Plugin* plugin) {
    plugin->creator = Ipqp::creator;
    plugin->name = "ipqp";
    plugin->doc = Ipqp::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Ipqp::options_;
    plugin->deserialize = &Ipqp::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_IPQP_EXPORT casadi_load_conic_ipqp() {
    Conic::registerPlugin(casadi_register_conic_ipqp);
  }

  void print_vec(const std::string& str, const double* v, casadi_int n) {
    uout() << str << " =";
    for (casadi_int k = 0; k < n; ++k) uout() << " " << v[k];
    uout() << "\n";
  }

  Ipqp::Ipqp(const std::string& name, const std::map<std::string, Sparsity> &st)
    : Conic(name, st) {
  }

  Ipqp::~Ipqp() {
    clear_mem();
  }

  const Options Ipqp::options_
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
        "A custom linear solver creator function [default: ldl]"}},
      {"linear_solver_options",
       {OT_DICT,
        "Options to be passed to the linear solver"}},
      {"min_lam",
       {OT_DOUBLE,
        "Smallest multiplier treated as inactive for the initial active set [0]."}}
     }
  };

  void Ipqp::init(const Dict& opts) {
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
    linear_solver_ = "ldl";
    // Read user options
    for (auto&& op : opts) {
      if (op.first=="max_iter") {
        p_.max_iter = op.second;
      } else if (op.first=="pr_tol") {
        p_.pr_tol = op.second;
      } else if (op.first=="du_tol") {
        p_.du_tol = op.second;
      } else if (op.first=="co_tol") {
        p_.co_tol = op.second;
      } else if (op.first=="mu_tol") {
        p_.mu_tol = op.second;
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
    // Memory for IP solver
    alloc_w(casadi_ipqp_sz_w(&p_), true);
    // Memory for KKT formation
    alloc_w(kkt_.nnz(), true);
    alloc_iw(A_.size2());
    alloc_w(nx_ + na_);
    // KKT solver
    linsol_ = Linsol("linsol", linear_solver_, kkt_, linear_solver_options_);
    // Print summary
    if (print_header_) {
      print("-------------------------------------------\n");
      print("This is casadi::Ipqp\n");
      print("Linear solver:                   %12s\n", linear_solver_.c_str());
      print("Number of variables:             %12d\n", nx_);
      print("Number of constraints:           %12d\n", na_);
      print("Number of nonzeros in H:         %12d\n", H_.nnz());
      print("Number of nonzeros in A:         %12d\n", A_.nnz());
      print("Number of nonzeros in KKT:       %12d\n", kkt_.nnz());
    }
  }

  void Ipqp::set_qp_prob() {
    casadi_ipqp_setup(&p_, nx_, na_);
  }

  int Ipqp::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<IpqpMemory*>(mem);
    m->return_status = "";
    return 0;
  }

  int Ipqp::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<IpqpMemory*>(mem);
    // Message buffer
    char buf[121];
    // Setup KKT system
    double* nz_kkt = w; w += kkt_.nnz();
    // Checkout a linear solver instance
    int linsol_mem = linsol_.checkout();
    // Setup IP solver
    casadi_ipqp_data<double> d;
    d.prob = &p_;
    casadi_ipqp_init(&d, &iw, &w);
    casadi_ipqp_bounds(&d, arg[CONIC_G],
      arg[CONIC_LBX], arg[CONIC_UBX], arg[CONIC_LBA], arg[CONIC_UBA]);
    casadi_ipqp_guess(&d, arg[CONIC_X0], arg[CONIC_LAM_X0], arg[CONIC_LAM_A0]);
    // Reverse communication loop
    while (casadi_ipqp(&d)) {
      switch (d.task) {
      case IPQP_MV:
        // Matrix-vector multiplication
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
    // Release linear solver instance
    linsol_.release(linsol_mem);
    // Read return status
    m->return_status = casadi_ipqp_return_status(d.status);
    if (d.status == IPQP_MAX_ITER)
      m->d_qp.unified_return_status = SOLVER_RET_LIMITED;
    // Get solution
    casadi_ipqp_solution(&d, res[CONIC_X], res[CONIC_LAM_X], res[CONIC_LAM_A]);
    if (res[CONIC_COST]) {
      *res[CONIC_COST] = .5 * casadi_bilin(arg[CONIC_H], H_, d.z, d.z)
        + casadi_dot(p_.nx, d.z, d.g);
    }
    // Return
    if (verbose_) casadi_warning(m->return_status);
    m->d_qp.success = d.status == IPQP_SUCCESS;
    return 0;
  }

  Dict Ipqp::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<IpqpMemory*>(mem);
    stats["return_status"] = m->return_status;
    return stats;
  }

  Ipqp::Ipqp(DeserializingStream& s) : Conic(s) {
    s.version("Ipqp", 1);
    s.unpack("Ipqp::kkt", kkt_);
    s.unpack("Ipqp::print_iter", print_iter_);
    s.unpack("Ipqp::print_header", print_header_);
    s.unpack("Ipqp::print_info", print_info_);
    s.unpack("Ipqp::linear_solver", linear_solver_);
    s.unpack("Ipqp::linear_solver_options", linear_solver_options_);
    set_qp_prob();
    s.unpack("Ipqp::max_iter", p_.max_iter);
    s.unpack("Ipqp::pr_tol", p_.pr_tol);
    s.unpack("Ipqp::du_tol", p_.du_tol);
    s.unpack("Ipqp::co_tol", p_.co_tol);
    s.unpack("Ipqp::mu_tol", p_.mu_tol);
  }

  void Ipqp::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);

    s.version("Ipqp", 1);
    s.pack("Ipqp::kkt", kkt_);
    s.pack("Ipqp::print_iter", print_iter_);
    s.pack("Ipqp::print_header", print_header_);
    s.pack("Ipqp::print_info", print_info_);
    s.pack("Ipqp::linear_solver", linear_solver_);
    s.pack("Ipqp::linear_solver_options", linear_solver_options_);
    s.pack("Ipqp::max_iter", p_.max_iter);
    s.pack("Ipqp::pr_tol", p_.pr_tol);
    s.pack("Ipqp::du_tol", p_.du_tol);
    s.pack("Ipqp::co_tol", p_.co_tol);
    s.pack("Ipqp::mu_tol", p_.mu_tol);
  }

} // namespace casadi
