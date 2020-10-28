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
      {"print_lincomb",
       {OT_BOOL,
        "Print dependant linear combinations of constraints [false]. "
        "Printed numbers are 0-based indices into the vector of [simple bounds;linear bounds]"}},
      {"min_lam",
       {OT_DOUBLE,
        "Smallest multiplier treated as inactive for the initial active set [0]."}}
     }
  };

  void Qpchasm::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    // Transpose of the Jacobian
    AT_ = A_.T();

    // Assemble KKT system sparsity
    kkt_ = Sparsity::kkt(H_, A_, true, true);

    // Symbolic QR factorization
    kkt_.qr_sparse(sp_v_, sp_r_, prinv_, pc_);

    // Setup memory structure
    set_qp_prob();

    // Default options
    print_iter_ = true;
    print_header_ = true;
    print_info_ = true;
    print_lincomb_ = false;

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
      } else if (op.first=="print_lincomb") {
        print_lincomb_ = op.second;
      }
    }

    // Allocate memory
    casadi_int sz_w, sz_iw;
    casadi_qp_work(&p_, &sz_iw, &sz_w);

    // D_x, D_g
    sz_w += p_.nz;
    // d.S
    sz_w += p_.nz;
    // lam_lbz, lam_ubz
    sz_w += p_.nz;
    sz_w += p_.nz;

    alloc_iw(sz_iw, true);
    alloc_w(sz_w, true);

    if (print_header_) {
      // Print summary
      print("-------------------------------------------\n");
      print("This is casadi::QPCHASM\n");
      print("Number of variables:                       %9d\n", nx_);
      print("Number of constraints:                     %9d\n", na_);
      print("Number of nonzeros in H:                   %9d\n", H_.nnz());
      print("Number of nonzeros in A:                   %9d\n", A_.nnz());
      print("Number of nonzeros in KKT:                 %9d\n", kkt_.nnz());
      print("Number of nonzeros in QR(V):               %9d\n", sp_v_.nnz());
      print("Number of nonzeros in QR(R):               %9d\n", sp_r_.nnz());
    }
  }

  void Qpchasm::set_qp_prob() {
    p_.sp_a = A_;
    p_.sp_h = H_;
    p_.sp_at = AT_;
    p_.sp_kkt = kkt_;
    p_.sp_v = sp_v_;
    p_.sp_r = sp_r_;
    p_.prinv = get_ptr(prinv_);
    p_.pc = get_ptr(pc_);
    casadi_qp_setup(&p_);
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
    casadi_qpip_data<double> d;
    d.prob = &p_;
    d.nz_h = arg[CONIC_H];
    d.g = arg[CONIC_G];
    d.nz_a = arg[CONIC_A];

    // D_x, D_a
    d.D = w; w += p_.nz;
    // d.S
    d.S = w; w += p_.nz;

    // lam_lbx, lam_ubz
    d.lam_lbz = w; w += p_.nz;
    d.lam_ubz = w; w += p_.nz;

    casadi_qp_init(&d, &iw, &w);
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

    // Total number of finite constraints
    casadi_int n_con = 0;

    // Find interior point
    for (casadi_int k = 0; k < p_.nz; ++k) {
      if (d.lbz[k] > -p_.inf) {
        if (d.ubz[k] < p_.inf) {
          // Both upper and lower bounds
          d.z[k] = .5 * (d.lbz[k] + d.ubz[k]);
          if (d.ubz[k] > d.lbz[k] + p_.dmin) {
            d.lam_lbz[k] = 1;
            d.lam_ubz[k] = 1;
            n_con += 2;
          }
        } else {
          // Only lower bound
          d.z[k] = d.lbz[k] + 1.;
          d.lam_lbz[k] = 1;
          n_con++;
        }
      } else {
        if (d.ubz[k] < p_.inf) {
          // Only upper bound
          d.z[k] = d.ubz[k] - 1.;
          d.lam_ubz[k] = 1;
          n_con++;
        } else {
          // Unbounded
          d.z[k] = 0;
        }
      }
    }
    uout() << n_con << " finite constraints\n";

    // Calculate complimentarity measure
    double mu = 0;
    for (casadi_int k = 0; k < p_.nz; ++k) {
      if (d.lbz[k] > -p_.inf) {
        mu += d.lam_lbz[k] * (d.z[k] - d.lbz[k]);
      }
      if (d.ubz[k] < p_.inf) {
        mu += d.lam_ubz[k] * (d.ubz[k] - d.z[k]);
      }
    }
    uout() << "mu = " << mu << "\n";

    // Calculate diagonal entries
    for (casadi_int k = 0; k < p_.nx; ++k) {
      if (d.ubz[k] <= d.lbz[k] + p_.dmin) {
        d.D[k] = -1;
        continue;
      }
      d.D[k] = 0;
      if (d.lbz[k] > -p_.inf) {
        d.D[k] += d.lam_lbz[k] / (d.z[k] - d.lbz[k]);
      }
      if (d.ubz[k] < p_.inf) {
        d.D[k] += d.lam_ubz[k] / (d.ubz[k] - d.z[k]);
      }
    }
    for (casadi_int k = p_.nx; k < p_.nz; ++k) {
      if (d.lbz[k] > -p_.inf) {
        if (d.ubz[k] < p_.inf) {
          if (d.ubz[k] - d.lbz[k] < p_.dmin) {
            // Fixed
            d.D[k] = 0;
          } else {
            // Upper and lower
            d.D[k] = 1. / (d.lam_lbz[k] / (d.z[k] - d.lbz[k])
              + d.lam_ubz[k] / (d.ubz[k] - d.z[k]));
          }
        } else {
          // Only lower
          d.D[k] = (d.z[k] - d.lbz[k]) / d.lam_lbz[k];
        }
      } else {
        if (d.ubz[k] < p_.inf) {
          // Only upper
          d.D[k] = (d.ubz[k] - d.z[k]) / d.lam_ubz[k];
        } else {
          // Neither upper or lower
          d.D[k] = -1;
        }
      }
    }

    uout() << "d.D =";
    for (casadi_int k = 0; k < p_.nz; ++k) uout() << " " << d.D[k];
    uout() << "\n";



    // Calculate scaling factors
    for (casadi_int k = 0; k < p_.nz; ++k) {
      if (d.D[k] < 0) {
        // Eliminate
        d.S[k] = 0;
        d.D[k] = 1;
      } else {
        // Scale
        d.S[k] = fmin(1., std::sqrt(1. / d.D[k]));
        d.D[k] = fmax(1., d.D[k]);
      }
    }

    // Calculate residual equations




    uout() << "d.D =";
    for (casadi_int k = 0; k < p_.nz; ++k) uout() << " " << d.D[k];
    uout() << "\n";
    uout() << "d.S =";
    for (casadi_int k = 0; k < p_.nz; ++k) uout() << " " << d.S[k];
    uout() << "\n";





    // Reset solver
    if (casadi_qp_reset(&d)) return 1;
    while (true) {
      // Prepare QP
      int flag = casadi_qp_prepare(&d);
      // Print iteration progress
      if (print_iter_) {
        if (d.iter % 10 == 0) {
          // Print header
          if (casadi_qp_print_header(&d, buf, sizeof(buf))) break;
          uout() << buf << "\n";
        }
        // Print iteration
        if (casadi_qp_print_iteration(&d, buf, sizeof(buf))) break;
        uout() << buf << "\n";
      }
      // Make an iteration
      flag = flag || casadi_qp_iterate(&d);
      // Print debug info
      if (print_lincomb_) {
        for (casadi_int k=0;k<d.sing;++k) {
          uout() << "lincomb: ";
          casadi_qp_print_colcomb(&d, buf, sizeof(buf), k);
          uout() << buf << "\n";
        }
      }
      if (flag) break;

      // User interrupt
      InterruptHandler::check();
    }
    // Check return flag
    switch (d.status) {
      case QP_SUCCESS:
        m->return_status = "success";
        break;
      case QP_MAX_ITER:
        m->return_status = "Maximum number of iterations reached";
        m->unified_return_status = SOLVER_RET_LIMITED;
        break;
      case QP_NO_SEARCH_DIR:
        m->return_status = "Failed to calculate search direction";
        break;
      case QP_PRINTING_ERROR:
        m->return_status = "Printing error";
        break;
    }
    // Get solution
    casadi_copy(&d.f, 1, res[CONIC_COST]);
    casadi_copy(d.z, nx_, res[CONIC_X]);
    casadi_copy(d.lam, nx_, res[CONIC_LAM_X]);
    casadi_copy(d.lam+nx_, na_, res[CONIC_LAM_A]);
    // Return
    if (verbose_) casadi_warning(m->return_status);
    m->success = d.status == QP_SUCCESS;
    return 0;
  }

  void Qpchasm::codegen_body(CodeGenerator& g) const {
    g.add_auxiliary(CodeGenerator::AUX_QP);
    if (print_iter_) g.add_auxiliary(CodeGenerator::AUX_PRINTF);
    g.local("d", "struct casadi_qp_data");
    g.local("p", "struct casadi_qp_prob");
    g.local("flag", "int");
    if (print_iter_ || print_lincomb_) g.local("buf[121]", "char");
    if (print_lincomb_) g.local("k", "casadi_int");

    // Setup memory structure
    g << "p.sp_a = " << g.sparsity(A_) << ";\n";
    g << "p.sp_h = " << g.sparsity(H_) << ";\n";
    g << "p.sp_at = " << g.sparsity(AT_) << ";\n";
    g << "p.sp_kkt = " << g.sparsity(kkt_) << ";\n";
    g << "p.sp_v = " << g.sparsity(sp_v_) << ";\n";
    g << "p.sp_r = " << g.sparsity(sp_r_) << ";\n";
    g << "p.prinv = " << g.constant(prinv_) << ";\n";
    g << "p.pc =  " << g.constant(pc_) << ";\n";
    g << "casadi_qp_setup(&p);\n";

    // Copy options
    g << "p.max_iter = " << p_.max_iter << ";\n";
    g << "p.min_lam = " << p_.min_lam << ";\n";
    g << "p.constr_viol_tol = " << p_.constr_viol_tol << ";\n";
    g << "p.dual_inf_tol = " << p_.dual_inf_tol << ";\n";

    // Setup data structure
    g << "d.prob = &p;\n";
    g << "d.nz_h = arg[" << CONIC_H << "];\n";
    g << "d.g = arg[" << CONIC_G << "];\n";
    g << "d.nz_a = arg[" << CONIC_A << "];\n";
    g << "casadi_qp_init(&d, &iw, &w);\n";

    g.comment("Pass bounds on z");
    g.copy_default(g.arg(CONIC_LBX), nx_, "d.lbz", "-casadi_inf", false);
    g.copy_default(g.arg(CONIC_LBA), na_, "d.lbz+" + str(nx_), "-casadi_inf", false);
    g.copy_default(g.arg(CONIC_UBX), nx_, "d.ubz", "casadi_inf", false);
    g.copy_default(g.arg(CONIC_UBA), na_, "d.ubz+" + str(nx_), "casadi_inf", false);

    g.comment("Pass initial guess");
    g.copy_default(g.arg(CONIC_X0), nx_, "d.z", "0", false);
    g << g.fill("d.z+"+str(nx_), na_, g.constant(nan)) << "\n";
    g.copy_default(g.arg(CONIC_LAM_X0), nx_, "d.lam", "0", false);
    g.copy_default(g.arg(CONIC_LAM_A0), na_, "d.lam+" + str(nx_), "0", false);

    g.comment("Solve QP");
    g << "if (casadi_qp_reset(&d)) return 1;\n";
    g << "while (1) {\n";
    g << "flag = casadi_qp_prepare(&d);\n";
    if (print_iter_) {
      // Print header
      g << "if (d.iter % 10 == 0) {\n";
      g << "if (casadi_qp_print_header(&d, buf, sizeof(buf))) break;\n";
      g << g.printf("%s\\n", "buf") << "\n";
      g << "}\n";
      // Print iteration
      g << "if (casadi_qp_print_iteration(&d, buf, sizeof(buf))) break;\n";
      g << g.printf("%s\\n", "buf") << "\n";
    }
    g << "if (flag || casadi_qp_iterate(&d)) break;\n";
    if (print_lincomb_) {
      g << "for (k=0;k<d.sing;++k) {\n";
      g << "casadi_qp_print_colcomb(&d, buf, sizeof(buf), k);\n";
      g << g.printf("lincomb: %s\\n", "buf") << "\n";
      g << "}\n";
    }
    g << "}\n";

    g.comment("Get solution");
    g.copy_check("&d.f", 1, g.res(CONIC_COST), false, true);
    g.copy_check("d.z", nx_, g.res(CONIC_X), false, true);
    g.copy_check("d.lam", nx_, g.res(CONIC_LAM_X), false, true);
    g.copy_check("d.lam+"+str(nx_), na_, g.res(CONIC_LAM_A), false, true);

    g << "return d.status != QP_SUCCESS;\n";
  }

  Dict Qpchasm::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<QpchasmMemory*>(mem);
    stats["return_status"] = m->return_status;
    return stats;
  }

  Qpchasm::Qpchasm(DeserializingStream& s) : Conic(s) {
    s.version("Qpchasm", 1);
    s.unpack("Qpchasm::AT", AT_);
    s.unpack("Qpchasm::kkt", kkt_);
    s.unpack("Qpchasm::sp_v", sp_v_);
    s.unpack("Qpchasm::sp_r", sp_r_);
    s.unpack("Qpchasm::prinv", prinv_);
    s.unpack("Qpchasm::pc", pc_);
    s.unpack("Qpchasm::print_iter", print_iter_);
    s.unpack("Qpchasm::print_header", print_header_);
    s.unpack("Qpchasm::print_info", print_info_);
    s.unpack("Qpchasm::print_lincomb_", print_lincomb_);
    set_qp_prob();
    s.unpack("Qpchasm::max_iter", p_.max_iter);
    s.unpack("Qpchasm::min_lam", p_.min_lam);
    s.unpack("Qpchasm::constr_viol_tol", p_.constr_viol_tol);
    s.unpack("Qpchasm::dual_inf_tol", p_.dual_inf_tol);
  }

  void Qpchasm::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);

    s.version("Qpchasm", 1);
    s.pack("Qpchasm::AT", AT_);
    s.pack("Qpchasm::kkt", kkt_);
    s.pack("Qpchasm::sp_v", sp_v_);
    s.pack("Qpchasm::sp_r", sp_r_);
    s.pack("Qpchasm::prinv", prinv_);
    s.pack("Qpchasm::pc", pc_);
    s.pack("Qpchasm::print_iter", print_iter_);
    s.pack("Qpchasm::print_header", print_header_);
    s.pack("Qpchasm::print_info", print_info_);
    s.pack("Qpchasm::print_lincomb_", print_lincomb_);
    s.pack("Qpchasm::max_iter", p_.max_iter);
    s.pack("Qpchasm::min_lam", p_.min_lam);
    s.pack("Qpchasm::constr_viol_tol", p_.constr_viol_tol);
    s.pack("Qpchasm::dual_inf_tol", p_.dual_inf_tol);
  }

} // namespace casadi
