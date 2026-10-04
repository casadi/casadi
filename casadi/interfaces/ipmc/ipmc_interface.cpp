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

#include "ipmc_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/global_options.hpp"
#include "casadi/core/stats_recorder_internal.hpp"
#include "casadi/casadi_c.h"

#include <cmath>
#include <cstring>
#include <algorithm>
#include <sstream>

#include <ipmc_runtime_str.h>

namespace casadi {
  #include "ipmc_base_runtime.hpp"

  extern "C"
  int CASADI_NLPSOL_IPMC_EXPORT
  casadi_register_nlpsol_ipmc(Nlpsol::Plugin* plugin) {
    plugin->creator = IpmcInterface::creator;
    plugin->name = "ipmc";
    plugin->doc = IpmcInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &IpmcInterface::options_;
    plugin->deserialize = &IpmcInterface::deserialize;
    // Slacks stay in the oracle, (x,p,s)->(f,g,f_s): structure detection sees the caller's OCP
    plugin->exposed.handles_slacks = true;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_IPMC_EXPORT casadi_load_nlpsol_ipmc() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_ipmc);
  }

  IpmcInterface::IpmcInterface(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }

  IpmcInterface::~IpmcInterface() {
    clear_mem();
  }

  // The gap-closing identity blocks: transition k's x_{k+1} columns of its dynamics rows
  Sparsity IpmcInterface::identity_sparsity() const {
    std::vector<casadi_ocp_block> I_blocks;
    for (const casadi_ocp_block& b : AB_blocks_) {
      I_blocks.push_back({b.offset_r, b.offset_c+b.cols, b.rows, b.rows});
    }
    return blocksparsity(nat_, nxt_, I_blocks, true);
  }

  Sparsity IpmcInterface::blocksparsity(casadi_int rows, casadi_int cols,
      const std::vector<casadi_ocp_block>& blocks, bool eye) {
    DM r(rows, cols);
    for (auto && b : blocks) {
      if (eye) {
        r(range(b.offset_r, b.offset_r+b.rows),
          range(b.offset_c, b.offset_c+b.cols)) = DM::eye(b.rows);
        casadi_assert_dev(b.rows==b.cols);
      } else {
        r(range(b.offset_r, b.offset_r+b.rows),
        range(b.offset_c, b.offset_c+b.cols)) = DM::zeros(b.rows, b.cols);
      }
    }
    return r.sparsity();
  }

  // IpmcError, and 99 from the mockup libipmc
  static std::string ipmc_soft_error_message(int code) {
    switch (code) {
      case IPMC_ERR_NS: return "a stage reported a negative number of slacks";
      case IPMC_ERR_SOFT_IDX:
        return "a slack index was out of range; this is an internal inconsistency of the "
               "interface, please report it";
      case IPMC_ERR_PENALTY:
        return "a penalty weight (z or Z) of the slack objective f_s is negative; "
               "ipmc needs a convex, non-decreasing penalty";
      case IPMC_ERR_SLACK_BOUND: return "a slack upper bound 'ubs' is not positive";
      case IPMC_ERR_NXC: return "the constant states 'nxc' do not fit every stage";
      case IPMC_ERR_SLACK_HELPER:
        return "the slack-helper description of a lifted problem was rejected; "
               "this is an internal inconsistency of the lift, please report it";
      case IPMC_ERR_BOUND_IDX:
        return "the simple-bound description was rejected; this is an internal "
               "inconsistency of the interface, please report it";
      case IPMC_ERR_MEMORY: return "the memory handed to ipmc is smaller than it asks for";
      case 99: return "the libipmc that was loaded is the mockup, ipmc's symbols "
                      "with no solver behind them. Put a licensed libipmc on the "
                      "library search path (LD_LIBRARY_PATH, PATH or "
                      "DYLD_LIBRARY_PATH) ahead of it";
      default: return "unknown reason";
    }
  }

  static void report_issue(casadi_int i, const std::string& msg) {
    casadi_int idx = i+GlobalOptions::start_index;
    casadi_warning("Structure detection error on row " + str(idx) + ". " + msg);
  }

  const Options IpmcInterface::options_
  = {{&Nlpsol::options_},
     {{"N",
       {OT_INT,
        "OCP horizon"}},
      {"nx",
       {OT_INTVECTOR,
        "Number of states, length N+1"}},
      {"nu",
       {OT_INTVECTOR,
        "Number of controls, length N+1"}},
      {"ng",
       {OT_INTVECTOR,
        "Number of non-dynamic constraints, length N+1"}},
      {"nxc",
       {OT_INT,
        "Number of trailing states per stage that are constant, x_{k+1}=x_k; ipmc "
        "exploits them in the Riccati recursion. Needs structure_detection 'manual' or "
        "'auto' [0]."}},
      {"ipmc",
       {OT_DICT,
        "Options to be passed to ipmc."}},
      {"structure_detection",
       {OT_STRING,
        "Structure detection: none, auto or manual [none]."}},
      {"debug",
       {OT_BOOL,
        "Write the expected and actual Jacobian structure to debug_ipmc_*.mtx [false]."}}
     }
  };

  void IpmcInterface::create_ipmc_functions() {
    const Function& orc = oracle_;
    create_function(orc, "nlp_f", {"x", "p"}, {"f"});
    create_function(orc, "nlp_g", {"x", "p"}, {"g"});
    if (!has_function("nlp_grad_f")) {
      create_function(orc, "nlp_grad_f", {"x", "p"}, {"grad:f:x"});
    }
    if (!has_function("nlp_jac_g")) {
      create_function(orc, "nlp_jac_g", {"x", "p"}, {"g", "jac:g:x"});
    }
    jacg_sp_ = get_function("nlp_jac_g").sparsity_out(1);
    if (exact_hessian_) {
      if (!has_function("nlp_hess_l")) {
        create_function(orc, "nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
                        {"grad:gamma:x", "hess:gamma:x:x"}, {{"gamma", {"f", "g"}}});
      }
      hesslag_sp_ = get_function("nlp_hess_l").sparsity_out(1);
      casadi_assert(hesslag_sp_.is_symmetric(), "Hessian must be symmetric");
    }
  }

  // A submatrix, by index lists
  static std::vector<casadi_int> ones(const std::vector<casadi_int>& v) {
    return std::vector<casadi_int>(v.size(), 1);
  }

  // The caller's partition: offsets of the variables, dynamics rows and path rows of stage k
  static void caller_offsets(casadi_int N, const std::vector<casadi_int>& nx,
      const std::vector<casadi_int>& nu, const std::vector<casadi_int>& ng,
      std::vector<casadi_int>& col_off, std::vector<casadi_int>& dyn_off,
      std::vector<casadi_int>& path_off) {
    col_off.assign(N+1, 0);
    dyn_off.assign(N+1, 0);
    path_off.assign(N+1, 0);
    casadi_int oc = 0, orow = 0;
    for (casadi_int k=0;k<N;++k) {
      col_off[k] = oc; oc += nx[k]+nu[k];
      dyn_off[k] = orow; path_off[k] = orow + nx[k+1];
      orow += nx[k+1]+ng[k];
    }
    col_off[N] = oc;
    dyn_off[N] = orow;
    path_off[N] = orow;
  }

  // Which column relaxes which side of which row; stage-local or lifted
  void IpmcInterface::slack_maps() {
    casadi_int K = N_+1;
    slack_g_lo_.assign(ng_, -1);
    slack_g_up_.assign(ng_, -1);
    slack_x_lo_.assign(nx_, -1);
    slack_x_up_.assign(nx_, -1);
    slack_perm_.clear();
    slack_idx_.clear();
    slack_ns_.assign(K, 0);
    lift_col_.clear();
    lift_ent_.clear();
    n_lift_ = 0;
    if (!slacks_) return;
    casadi_int ns = Nlpsol::ns_, si = GlobalOptions::start_index;
    // The stage of every path row and variable; -1 for a dynamics row
    std::vector<casadi_int> col_off, dyn_off, path_off;
    caller_offsets(N_, nxs_, nus_, ngs_, col_off, dyn_off, path_off);
    std::vector<casadi_int> g_stage(ng_, -1), x_stage(nx_, -1);
    for (casadi_int k=0;k<K;++k) {
      for (casadi_int i=0;i<ngs_[k];++i) g_stage[path_off[k]+i] = k;
      for (casadi_int i=0;i<nxs_[k]+nus_[k];++i) x_stage[col_off[k]+i] = k;
    }
    // The stage of every column: -2 none yet, -1 several
    std::vector<casadi_int> col_stage(ns, -2);
    for (casadi_int up=0; up<2; ++up) {
      // Rows of slack_lo_/slack_up_ are [g; x], unlike z = [x; g]
      const Sparsity& B = up ? slack_up_ : slack_lo_;
      std::string sd = up ? "upper" : "lower";
      const casadi_int* colind = B.colind();
      const casadi_int* row = B.row();
      for (casadi_int j=0;j<ns;++j) {
        for (casadi_int el=colind[j]; el<colind[j+1]; ++el) {
          casadi_int r = row[el];
          bool is_g = r<ng_;
          casadi_int i = is_g ? r : r-ng_;
          std::vector<casadi_int>& map = is_g ? (up ? slack_g_up_ : slack_g_lo_)
                                              : (up ? slack_x_up_ : slack_x_lo_);
          casadi_int k = is_g ? g_stage[i] : x_stage[i];
          std::string what = is_g ? "the " + sd + " side of constraint row g[" + str(i+si) + "]"
                                  : "the " + sd + " bound on x[" + str(i+si) + "]";
          casadi_assert(map[i]<0,
            "Ipmc: " + what + " is relaxed by more than one slack column (" +
            str(map[i]+si) + " and " + str(j+si) + "). ipmc relaxes each side of a row by "
            "at most one slack, so every row of S_lo and of S_up may hold at most one entry.");
          map[i] = j;
          casadi_assert(k>=0,
            "Ipmc: slack column " + str(j+si) + " relaxes " + what + ", which is a "
            "gap-closing (dynamics) row of the detected OCP structure. Only path "
            "constraints and simple bounds can be softened. Structure is: N " + str(N_) +
            ", nx " + str(nxs_) + ", nu " + str(nus_) + ", ng " + str(ngs_) + ".");
          if (col_stage[j]==-2) {
            col_stage[j] = k;
          } else if (col_stage[j]!=k) {
            col_stage[j] = -1;
          }
        }
      }
    }
    // ipmc's slacks: the stage-local columns, stage after stage. A column relaxing nothing
    // is a slack of stage 0; a column spanning several stages is lifted.
    for (casadi_int& k : col_stage) if (k==-2) k = 0;
    std::vector<casadi_int> offs(K+1, 0);
    for (casadi_int j=0;j<ns;++j) if (col_stage[j]>=0) slack_ns_[col_stage[j]]++;
    for (casadi_int k=0;k<K;++k) offs[k+1] = offs[k] + slack_ns_[k];
    slack_perm_.assign(offs[K], -1);
    lift_ent_.assign(ns, -1);
    for (casadi_int j=0;j<ns;++j) {
      if (col_stage[j]>=0) {
        slack_perm_[offs[col_stage[j]]++] = j;
      } else {
        lift_ent_[j] = n_lift_++;
        lift_col_.push_back(j);
      }
    }
    slack_idx_.assign(ns, -1);
    for (casadi_int i=0;i<slack_perm_.size();++i) slack_idx_[slack_perm_[i]] = i;
    if (verbose_) {
      casadi_message("Native slacks: ns " + str(ns) + ", per stage " + str(slack_ns_) +
        ", lifted into helper states " + str(n_lift_) + ".");
    }
  }

  // The handed-over problem: the caller's, with every lifted column a helper state
  void IpmcInterface::lift() {
    casadi_int K = N_+1, ne = n_lift_;
    const std::vector<casadi_int>& nx = nxs_;
    const std::vector<casadi_int>& nu = nus_;
    const std::vector<casadi_int>& ng = ngs_;
    std::vector<casadi_int> col_off, dyn_off, path_off;
    caller_offsets(N_, nx, nu, ng, col_off, dyn_off, path_off);
    // The helper lifting the lower / upper side of caller z entry c, -1 if none
    auto helper = [&](casadi_int c, bool up) -> casadi_int {
      if (!slacks_) return -1;
      casadi_int col = c<nx_ ? (up ? slack_x_up_ : slack_x_lo_)[c]
                             : (up ? slack_g_up_ : slack_g_lo_)[c-nx_];
      return col>=0 ? lift_ent_[col] : -1;
    };
    nxh_ = nx;
    ngh_ = ng;
    lift_src_.clear();
    lift_side_.clear();
    lift_hlp_.clear();
    if (ne==0) {
      nxt_ = nx_;
      nat_ = ng_;
      lift_src_ = range(nx_+ng_);
      lift_side_.assign(nx_+ng_, LIFT_BOTH);
      lift_hlp_.assign(nx_+ng_, -1);
      lift_Px_ = lift_H_ = lift_Pg_ = lift_C_ = lift_B_ = lift_Blo_ = lift_Bup_ = lift_Gz_ = DM();
      lift_lbz0_.clear();
      lift_ubz0_.clear();
      lift_m0_.clear();
      return;
    }
    // Is a caller variable's bound lifted, so that the variable is freed
    auto freed = [&](casadi_int xi) { return helper(xi, false)>=0 || helper(xi, true)>=0; };

    // Variables [x_0, m, u_0, ..., x_N, m]: trailing m makes the helpers constant states
    std::vector<casadi_int> vsrc, vside, px_r, px_c, h_r, h_c, mcopy(K*ne);
    casadi_int a = 0;
    for (casadi_int k=0;k<K;++k) {
      for (casadi_int j=0;j<nx[k]+nu[k];++j) {
        if (j==nx[k]) {
          for (casadi_int e=0;e<ne;++e) {
            mcopy[k*ne+e] = a;
            h_r.push_back(a); h_c.push_back(e);
            vsrc.push_back(-1); vside.push_back(LIFT_NONE);
            a++;
          }
        }
        casadi_int xi = col_off[k]+j;
        px_r.push_back(a); px_c.push_back(xi);
        vsrc.push_back(xi); vside.push_back(freed(xi) ? LIFT_NONE : LIFT_BOTH);
        a++;
      }
      if (nu[k]==0) {
        for (casadi_int e=0;e<ne;++e) {
          mcopy[k*ne+e] = a;
          h_r.push_back(a); h_c.push_back(e);
          vsrc.push_back(-1); vside.push_back(LIFT_NONE);
          a++;
        }
      }
    }
    nxt_ = a;
    lift_Px_ = IM::triplet(px_r, px_c, ones(px_r), nxt_, nx_);
    lift_H_ = IM::triplet(h_r, h_c, ones(h_r), nxt_, ne);
    lift_m0_.assign(mcopy.begin(), mcopy.begin()+ne);

    // Rows, counted by r; m_* the helper terms
    std::vector<casadi_int> rsrc, rside, rhlp, m_r, m_c, m_v;
    // A row carrying caller z entry c, all of it or one half of a split
    auto add = [&](casadi_int k, casadi_int c, casadi_int side, casadi_int e, bool neg) {
      if (e>=0) {
        m_r.push_back(rsrc.size()); m_c.push_back(mcopy[k*ne+e]); m_v.push_back(neg ? -1 : 1);
      }
      rsrc.push_back(c); rside.push_back(side); rhlp.push_back(e>=0 ? 4*e+(neg ? 1 : 0) : -1);
    };
    // A row with a lifted side splits: the lower half c + m_lo, the upper half c - m_up
    auto add_split = [&](casadi_int k, casadi_int c) {
      add(k, c, LIFT_LO, helper(c, false), false);
      add(k, c, LIFT_UP, helper(c, true), true);
      ngh_[k]++;
    };
    for (casadi_int k=0;k<K;++k) {
      if (k<N_) {
        // Dynamics rows of transition k, then m_{k+1} - m_k = 0
        for (casadi_int j=0;j<nx[k+1];++j) add(k, nx_+dyn_off[k]+j, LIFT_BOTH, -1, false);
        for (casadi_int e=0;e<ne;++e) {
          m_r.push_back(rsrc.size()); m_c.push_back(mcopy[k*ne+e]); m_v.push_back(-1);
          m_r.push_back(rsrc.size()); m_c.push_back(mcopy[(k+1)*ne+e]); m_v.push_back(1);
          rsrc.push_back(-1); rside.push_back(LIFT_NONE); rhlp.push_back(-1);
        }
      }
      for (casadi_int i=0;i<ng[k];++i) {
        casadi_int c = nx_+path_off[k]+i;
        if (helper(c, false)>=0 || helper(c, true)>=0) {
          add_split(k, c);
        } else {
          add(k, c, LIFT_BOTH, -1, false);
        }
      }
      // A lifted simple bound becomes a split path row; its own bound is freed
      for (casadi_int i=0;i<nx[k]+nu[k];++i) {
        casadi_int xi = col_off[k]+i;
        if (freed(xi)) {
          add_split(k, xi);
          ngh_[k]++;
        }
      }
    }
    nat_ = rsrc.size();
    for (casadi_int k=0;k<K;++k) nxh_[k] += ne;

    // Per handed-over z entry
    lift_src_ = vsrc;
    lift_src_.insert(lift_src_.end(), rsrc.begin(), rsrc.end());
    lift_side_ = vside;
    lift_side_.insert(lift_side_.end(), rside.begin(), rside.end());
    lift_hlp_.assign(nxt_, -1);
    lift_hlp_.insert(lift_hlp_.end(), rhlp.begin(), rhlp.end());

    // g = Pg g_u + C x, C: the variables a row carries (lifted bounds) and the helper terms
    std::vector<casadi_int> pg_r, pg_c, c_r, c_c, c_v;
    std::vector<casadi_int> pos(nx_);
    for (casadi_int i=0;i<px_r.size();++i) pos[px_c[i]] = px_r[i];
    for (casadi_int r=0;r<nat_;++r) {
      casadi_int c = rsrc[r];
      if (c>=nx_) {
        pg_r.push_back(r); pg_c.push_back(c-nx_);
      } else if (c>=0) {
        c_r.push_back(r); c_c.push_back(pos[c]); c_v.push_back(1);
      }
    }
    IM C = IM::triplet(c_r, c_c, c_v, nat_, nxt_) + IM::triplet(m_r, m_c, m_v, nat_, nxt_);
    lift_Pg_ = IM::triplet(pg_r, pg_c, ones(pg_r), nat_, ng_);
    lift_C_ = C;

    // B, Blo, Bup and the bounds lbz0/ubz0 of whatever carries no caller bound
    std::vector<casadi_int> b_r, b_c, lo_r, lo_c, up_r, up_c, last(ng_, -1);
    lift_lbz0_.assign(nxt_+nat_, 0);
    lift_ubz0_.assign(nxt_+nat_, 0);
    for (casadi_int t=0;t<nxt_+nat_;++t) {
      casadi_int c = lift_src_[t], side = lift_side_[t];
      if (t<nxt_ && side==LIFT_NONE) {
        lift_lbz0_[t] = -inf;
        lift_ubz0_[t] = inf;
      }
      if (c<0 || side==LIFT_NONE) continue;
      if (t>=nxt_ && c>=nx_) last[c-nx_] = t-nxt_;
      b_r.push_back(t); b_c.push_back(c);
      if (side!=LIFT_UP) {
        lo_r.push_back(t); lo_c.push_back(c);
      } else {
        lift_lbz0_[t] = -inf;
      }
      if (side!=LIFT_LO) {
        up_r.push_back(t); up_c.push_back(c);
      } else {
        lift_ubz0_[t] = inf;
      }
    }
    lift_B_ = IM::triplet(b_r, b_c, ones(b_r), nxt_+nat_, nx_+ng_);
    lift_Blo_ = IM::triplet(lo_r, lo_c, ones(lo_r), nxt_+nat_, nx_+ng_);
    lift_Bup_ = IM::triplet(up_r, up_c, ones(up_r), nxt_+nat_, nx_+ng_);
    // Constraint values back: G picks one rewritten row per caller row
    IM G = IM::triplet(range(ng_), last, ones(last), ng_, nat_);
    lift_Gz_ = mtimes(G, horzcat(-C, IM::eye(nat_)));
    if (verbose_) {
      casadi_message("Lifted " + str(ne) + " cross-stage slack column(s) into helper "
        "states; rewritten problem has nx " + str(nxt_) + ", ng " + str(nat_) + ".");
    }
  }

  // Codes: caller nonzero el is el+1, a +-1 of the lift +-one, helper e's penalty zero+1+e
  IM IpmcInterface::jac_codes() const {
    casadi_int one = jacg_sp_.nnz()+1;
    IM J(jacg_sp_, IM(range(1, one)));
    if (n_lift_==0) return J;
    IM B = mtimes(IM(lift_Pg_), mtimes(J, IM(lift_Px_).T()));
    IM A = B + one*IM(lift_C_);
    // A row carries either a caller row or caller variables and helper terms
    casadi_assert_dev(A.nnz()==B.nnz()+lift_C_.nnz());
    return A;
  }

  IM IpmcInterface::hess_codes() const {
    casadi_int zero = hesslag_sp_.nnz()+1;
    IM H(hesslag_sp_, IM(range(1, zero)));
    if (n_lift_==0) return H;
    IM Px = lift_Px_;
    const std::vector<casadi_int>& m0 = lift_m0_;
    IM B = mtimes(Px, mtimes(H, Px.T()));
    IM P = IM::triplet(m0, m0, range(zero+1, zero+1+n_lift_), nxt_, nxt_);
    casadi_assert_dev((B+P).nnz()==B.nnz()+P.nnz());
    return B+P;
  }

  static IM sub(const IM& M, const std::vector<casadi_int>& rr,
      const std::vector<casadi_int>& cc) {
    IM r;
    M.get(r, false, IM(rr), IM(cc));
    return r;
  }

  // Blocks side by side, rows padded to nr; column j of block i to column j of blk[i]
  static void stack(const std::vector<IM>& B, const std::vector<casadi_int>& blk,
      casadi_int nr, Sparsity& sp, std::vector<casadi_int>& code,
      std::vector<casadi_int>& tblk, std::vector<casadi_int>& tcol) {
    std::vector<IM> padded;
    tblk.clear();
    tcol.clear();
    for (size_t i=0;i<B.size();++i) {
      padded.push_back(vertcat(B[i], IM(nr-B[i].size1(), B[i].size2())));
      for (casadi_int j=0;j<B[i].size2();++j) {
        tblk.push_back(blk[i]);
        tcol.push_back(j);
      }
    }
    IM T = padded.empty() ? IM(nr, 0) : horzcat(padded);
    sp = T.sparsity();
    code = T.nonzeros();
  }

  // ipmc's stage blocks, sliced from the code matrices in ipmc's order [u_k; x_k sans helpers]
  void IpmcInterface::build_pack_tables() {
    casadi_int K = N_+1, ne = n_lift_, nc = nxc_+n_lift_;
    // ux_k, as columns of the handed-over problem
    std::vector< std::vector<casadi_int> > ux(K);
    casadi_int nr_max = 0;
    for (casadi_int k=0;k<K;++k) {
      casadi_int c0 = CD_blocks_[k].offset_c;
      ux[k] = range(c0+nxh_[k], c0+nxh_[k]+nus_[k]);
      for (casadi_int j=0;j<nxh_[k]-ne;++j) ux[k].push_back(c0+j);
      nr_max = std::max(nr_max, static_cast<casadi_int>(ux[k].size()));
    }

    // Jacobian; sources au, then 1.0 (code one)
    casadi_int one = jacg_sp_.nnz()+1;
    IM A = jac_codes();

    // BAt[k] = -(non-constant dynamics rows of transition k)': ipmc's rows are x_{k+1} - F
    std::vector<IM> blocks;
    std::vector<casadi_int> blk;
    cchk_.clear();
    for (casadi_int k=0;k<N_;++k) {
      casadi_int r0 = AB_blocks_[k].offset_r, nxd1 = nxh_[k+1]-nc, nxd0 = nxh_[k]-nc;
      casadi_int c0 = CD_blocks_[k].offset_c;
      blocks.push_back(-sub(A, range(r0, r0+nxd1), ux[k]).T());
      blk.push_back(k);
      // Constant-state rows are not handed over: -1 on the own state, 0 elsewhere
      IM C = sub(A, range(r0+nxd1, r0+nxh_[k+1]), range(c0, c0+nxh_[k]+nus_[k]));
      for (casadi_int i=0;i<nc;++i) {
        casadi_assert(C.sparsity().has_nz(i, nxd0+i), "nxc: stage " + str(k+1) + " state "
          + str(nxd1+i) + " was declared constant, but its dynamics row is not x_{k+1}=x_k.");
      }
      const casadi_int* colind = C.colind();
      const casadi_int* row = C.row();
      for (casadi_int cl=0;cl<C.size2();++cl) {
        for (casadi_int el=colind[cl];el<colind[cl+1];++el) {
          casadi_int i = row[el], c = C.nonzeros()[el], ev = cl==nxd0+i ? -1 : 0;
          if (c>0 && c!=one) {
            cchk_.insert(cchk_.end(), {c-1, ev, k, nxd1+i});
          } else {
            // The lift's own entries are right by construction
            casadi_assert_dev((c==one ? 1 : -1)==ev);
          }
        }
      }
    }
    stack(blocks, blk, nr_max, bat_sp_, bat_code_, bat_blk_, bat_col_);

    // The gap-closing identity: every entry of the I blocks must be 1
    ichk_.clear();
    IM I = project(A, identity_sparsity());
    for (casadi_int c : I.nonzeros()) {
      if (c>0 && c!=one) {
        ichk_.push_back(c-1);
      } else {
        casadi_assert_dev(c==0 || c==one);  // 0: absent
      }
    }

    // Gt_ineq: the stored rows are the path rows; a simple bound has no column
    blocks.clear();
    blk.clear();
    for (casadi_int k=0;k<K;++k) {
      casadi_int r0 = CD_blocks_[k].offset_r;
      blocks.push_back(sub(A, range(r0, r0+CD_blocks_[k].rows), ux[k]).T());
      blk.push_back(k);
    }
    stack(blocks, blk, nr_max, gi_sp_, gi_code_, gi_blk_, gi_col_);

    // Hessian; sources hu, then 0.0 (code zero), then helper e's penalty (code zero+1+e)
    blocks.clear();
    blk.clear();
    std::vector<IM> slack;
    if (exact_hessian_) {
      casadi_int zero = hesslag_sp_.nnz()+1;
      IM H = hess_codes();
      for (casadi_int k=0;k<K;++k) {
        casadi_int c0 = CD_blocks_[k].offset_c, nr = ux[k].size();
        // RSQ[k] with a full diagonal: restoration adds to it in place
        IM R = sub(H, ux[k], ux[k]);
        R = project(R, R.sparsity() + Sparsity::diag(nr));
        for (casadi_int& c : R.nonzeros()) if (c==0) c = zero;
        blocks.push_back(R);
        blk.push_back(k);
        // RSQ_slack[k]: the helpers' diagonal, dense for the same reason
        std::vector<casadi_int> hk = range(c0+nxh_[k]-ne, c0+nxh_[k]);
        IM D = project(diag(sub(H, hk, hk)), Sparsity::dense(ne, 1));
        for (casadi_int& c : D.nonzeros()) if (c==0) c = zero;
        slack.push_back(D);
      }
    }
    stack(blocks, blk, nr_max, rsq_sp_, rsq_code_, rsq_blk_, rsq_col_);
    IM S = slack.empty() ? IM(ne, 0) : horzcat(slack);
    rsqs_sp_ = S.sparsity();
    rsqs_code_ = S.nonzeros();
  }

  // Per stage the path rows, then the variables as simple bounds; all of them inequality rows
  void IpmcInterface::build_rows() {
    casadi_int K = N_+1;
    // Stage-local slacks of stage k start at soft_offs[k]
    std::vector<casadi_int> soft_offs(K+1, 0);
    for (casadi_int k=0;slacks_ && k<K;++k) soft_offs[k+1] = soft_offs[k] + slack_ns_[k];
    ineq_z_.clear();
    ineq_lo_.clear();
    ineq_up_.clear();
    d_ng_.assign(K, 0);
    d_ng_ineq_.assign(K, 0);
    d_nb_.assign(K, 0);
    d_ns_.assign(K, 0);
    d_idxb_.clear();
    d_soft_lo_.clear();
    d_soft_up_.clear();
    d_slack_helper_.clear();
    // The stage-local slack of one side of z entry t of stage k, if t carries that side
    auto side = [&](casadi_int k, casadi_int t, bool up, std::vector<casadi_int>& ineq,
        std::vector<ipmc_int>& code) {
      casadi_int c = lift_src_[t], s = lift_side_[t], i = -1;
      if (slacks_ && c>=0 && s!=LIFT_NONE && s!=(up ? LIFT_LO : LIFT_UP)) {
        casadi_int col = c<nx_ ? (up ? slack_x_up_ : slack_x_lo_)[c]
                               : (up ? slack_g_up_ : slack_g_lo_)[c-nx_];
        if (col>=0) i = slack_idx_[col];
      }
      ineq.push_back(i);
      code.push_back(i<0 ? -1 : i-soft_offs[k]);
    };
    // Row t (z-space) of stage k; comp>=0: a simple bound on component comp of [u_k; x_k]
    auto add_row = [&](casadi_int k, casadi_int t, casadi_int comp) {
      casadi_int code = lift_hlp_[t];
      // The upper half of a split row follows its lower half: twins if both carry a helper
      if (code>=0 && code%2==1 && !ineq_z_.empty() && ineq_z_.back()==t-1
          && lift_src_[t-1]==lift_src_[t] && lift_hlp_[t-1]>=0) code += 2;
      ineq_z_.push_back(t);
      side(k, t, false, ineq_lo_, d_soft_lo_);
      side(k, t, true, ineq_up_, d_soft_up_);
      d_slack_helper_.push_back(code);
      d_ng_ineq_[k]++;
      if (comp>=0) {
        d_idxb_.push_back(comp);
        d_nb_[k]++;
      }
    };
    for (casadi_int k=0;k<K;++k) {
      for (casadi_int i=0;i<CD_blocks_[k].rows;++i) {
        add_row(k, nxt_+CD_blocks_[k].offset_r+i, -1);
      }
      // The simple bounds
      for (casadi_int j=0;j<CD_blocks_[k].cols;++j) {
        casadi_int a = CD_blocks_[k].offset_c+j;
        add_row(k, a, j<nxh_[k] ? nus_[k]+j : j-nxh_[k]);
      }
      if (slacks_) d_ns_[k] = slack_ns_[k];
    }
    d_nu_.assign(nus_.begin(), nus_.end());
    d_nx_.assign(nxh_.begin(), nxh_.end());
  }

  // The handed-over partition and its blocks
  void IpmcInterface::build() {
    slack_maps();
    lift();
    const std::vector<casadi_int>& nx = nxh_;
    const std::vector<casadi_int>& nu = nus_;
    const std::vector<casadi_int>& ng = ngh_;
    AB_blocks_.clear();
    CD_blocks_.clear();
    casadi_int offset_r = 0, offset_c = 0;
    for (casadi_int k=0;k<N_;++k) {
      AB_blocks_.push_back({offset_r,        offset_c,            nx[k+1], nx[k]+nu[k]});
      CD_blocks_.push_back({offset_r+nx[k+1], offset_c,           ng[k], nx[k]+nu[k]});
      offset_c+= nx[k]+nu[k];
      offset_r+= nx[k+1]+ng[k];
    }
    CD_blocks_.push_back({offset_r, offset_c,           ng[N_], nx[N_]+nu[N_]});
  }

  void IpmcInterface::detect_structure(std::set<casadi_int>& errors) {
    const Sparsity& A_ = jacg_sp_;
    casadi_int na_ = ng_;
    const std::vector<casadi_int>& nx = nxs_;
    const std::vector<casadi_int>& ng = ngs_;
    const std::vector<casadi_int>& nu = nus_;
    casadi_assert(!equality_.empty(),
      "Structure detection auto requires the 'equality' option to be set");
    // General strategy: look for the x_{k+1} diagonal part in A

    // Find the right-most column for each row in A -> A_skyline
    // Find the second-to-right-most column -> A_skyline2
    // Find the left-most column -> A_bottomline
    Sparsity AT = A_.T();
    std::vector<casadi_int> A_skyline;
    std::vector<casadi_int> A_skyline2;
    std::vector<casadi_int> A_bottomline;

    std::vector<casadi_int> AT_colind = AT.get_colind();
    std::vector<casadi_int> AT_row = AT.get_row();
    for (casadi_int i=0;i<AT.size2();++i) {
      casadi_int pivot = AT_colind.at(i+1);
      if (pivot>AT_colind.at(i)) {
        A_bottomline.push_back(AT_row.at(AT_colind.at(i)));
      } else {
        A_bottomline.push_back(-1);
      }
      if (pivot>AT_colind.at(i)) {
        A_skyline.push_back(AT_row.at(pivot-1));
        if (pivot>AT_colind.at(i)+1) {
          A_skyline2.push_back(AT_row.at(pivot-2));
        } else {
          A_skyline2.push_back(-1);
        }
      } else {
        A_skyline.push_back(-1);
        A_skyline2.push_back(-1);
      }
    }

    casadi_assert(equality_[0],
     "Constraint Jacobian must start with gap-closing constraint "
     "(tagged 'true' in equality vector).");

    casadi_int pivot = A_skyline[0]; // Current right-most element
    casadi_int start_pivot = pivot; // First right-most element that started the stage
    casadi_int prev_start_pivot = 0;

    bool walking = true;

    nxs_.push_back(1);
    nus_.push_back(0);
    ngs_.push_back(0);
    for (casadi_int i=1;i<na_;++i) { // Loop over all rows
      bool is_gap_closing = true;
      if (A_bottomline[i]!=-1 && A_bottomline[i]<prev_start_pivot) {
        errors.insert(i);
        report_issue(i, "Constraint found depending on a state of the previous interval.");
      }
      if (equality_[i]) {
        // A candidate for a gap-closing constraint must tagged as equality
        if (A_skyline[i]>pivot+1) { // Jump to a diagonal in the future
          if (A_bottomline[i]!=-1 && A_bottomline[i]<start_pivot) {
            errors.insert(i);
            report_issue(i, "Constraint found depending on a state of the previous interval.");
          }
          nxs_.push_back(1);
          nus_.push_back(A_skyline[i]-pivot-1); // Size of jump equals number of states
          ngs_.push_back(0);
          prev_start_pivot = start_pivot;
          start_pivot = A_skyline[i];
          pivot = A_skyline[i];
          walking = true;
        } else if (A_skyline[i]==pivot+1) { // Walking the diagonal
          if (A_skyline2[i]<start_pivot) { // Free of below-diagonal entries?
            if (A_bottomline[i]>=prev_start_pivot) { // We must depend on at least one state
              pivot++;
              nxs_.back()++;
              walking = true;
            } else {
              if (A_bottomline[i]!=-1 && A_bottomline[i]<start_pivot) {
                errors.insert(i);
                report_issue(i, "Constraint found depending "
                  "on a state of the previous interval.");
              }
              is_gap_closing = false;
            }
          } else {
            nxs_.push_back(1);
            nus_.push_back(0);
            ngs_.push_back(0);
            if (A_bottomline[i]!=-1 && A_bottomline[i]<start_pivot) {
              errors.insert(i);
              report_issue(i, "Gap-closing constraint found depending "
                "on a state of the previous interval.");
            }
            prev_start_pivot = start_pivot;
            start_pivot = A_skyline[i];
            pivot = A_skyline[i];
            walking = true;
          }
        } else {
          is_gap_closing = false;
        }
      } else {
        is_gap_closing = false;
      }

      if (!is_gap_closing) {
        if (walking) {
          if (A_skyline[i]>=start_pivot) {
            nxs_.push_back(0);
            nus_.push_back(0);
            ngs_.push_back(0);
            walking = false;
          }
        }
        ngs_.back()++; // non-gap-closing constraint detected
      }

    }

    if (nxs_.back()!=0) {
      nxs_.push_back(0);
      nus_.push_back(0);
      ngs_.push_back(0);
    }

    // Set nx0==nx1 unless not allowed
    nxs_.insert(nxs_.begin(), std::min(A_skyline[0], nxs_.front()));

    // Patch loose ends
    nus_.front() += std::max(A_skyline[0]-nxs_.front(), static_cast<casadi_int>(0));
    nus_.back() += nx_-sum(nu)-sum(nx);

    casadi_assert_dev(nxs_.back()==0);
    nxs_.pop_back();

    casadi_assert_dev(nx.size()==nu.size());
    casadi_assert_dev(nx.size()==ng.size());

    casadi_assert_dev(sum(ng)+sum(nx)==na_+nx.front());
    casadi_assert_dev(sum(nx)+sum(nu)==nx_);

    N_ = nxs_.size()-1;
  }

  void IpmcInterface::settle_slack_penalty() {
    slacks_ = slack_native_ && Nlpsol::ns_ > 0;
    fs_z_.clear();
    fs_Z_.clear();
    if (slacks_) {
      // Not create_function: used once here, so no codegen dependency
      Function fs = oracle_.factory(name_ + "_fs",
        {"s", "p"}, {"f_s", "grad:f_s:s", "hess:f_s:s:s"});
      casadi_assert(!fs.has_free(),
        "Ipmc: the slack penalty f_s depends on " + str(fs.get_free()) + ", which is "
        "neither s nor p. Pass expand_slacks=True to fall back to the reference "
        "expansion.");
      Sparsity hsp = fs.sparsity_out(2);
      const casadi_int* hcolind = hsp.colind();
      const casadi_int* hrow = hsp.row();
      for (casadi_int c=0; c<hsp.size2(); ++c) {
        for (casadi_int el=hcolind[c]; el<hcolind[c+1]; ++el) {
          casadi_assert(hrow[el]==c,
            "Ipmc: the slack penalty f_s is not separable: hess(f_s,s,s) has a "
            "structural entry at (" + str(hrow[el]) + ", " + str(c) + "). ipmc can only "
            "represent a penalty of the form sum_j (z_j s_j + 1/2 Z_j s_j^2), so every "
            "slack must be penalised on its own. Use a separable quadratic (L1/L2/Huber) "
            "penalty, or pass expand_slacks=True to fall back to the reference expansion.");
        }
      }
      casadi_assert(fs.jac_sparsity(2, 0).nnz()==0,
        "Ipmc: the slack penalty f_s is not quadratic in s: hess(f_s,s,s) still depends "
        "on s. ipmc can only represent sum_j (z_j s_j + 1/2 Z_j s_j^2). Use a quadratic "
        "penalty, or pass expand_slacks=True to fall back to the reference expansion.");
      casadi_assert(fs.jac_sparsity(1, 1).nnz()==0 && fs.jac_sparsity(2, 1).nnz()==0,
        "Ipmc: the slack penalty f_s depends on the parameter p. z and Z are evaluated "
        "once, when the solver is built, and handed to ipmc as numbers, so a penalty "
        "retuned through p cannot be honoured: the solve would keep using the weights p "
        "held at construction. Make the weights literal constants, rebuild the solver "
        "when they change, or pass expand_slacks=True to fall back to the reference "
        "expansion, which carries f_s into the NLP itself and therefore tracks p.");
      // z = grad(f_s,s) at s=0, Z = diag hess(f_s,s,s)
      std::vector<DM> fs_arg = {DM::zeros(fs.size1_in(0), fs.size2_in(0)),
                                DM::zeros(fs.size1_in(1), fs.size2_in(1))};
      std::vector<DM> fs_res = fs(fs_arg);
      fs_z_ = densify(fs_res[1]).nonzeros();
      fs_Z_ = densify(diag(fs_res[2])).nonzeros();
      casadi_assert_dev(fs_z_.size()==static_cast<size_t>(Nlpsol::ns_));
      casadi_assert_dev(fs_Z_.size()==static_cast<size_t>(Nlpsol::ns_));
      for (casadi_int i=0; i<Nlpsol::ns_; ++i) {
        casadi_assert(std::isfinite(fs_z_[i]) && std::isfinite(fs_Z_[i]),
          "Ipmc: the slack penalty f_s has a non-finite gradient or Hessian at s=0 "
          "(slack column " + str(i+GlobalOptions::start_index) + "). Use a penalty that "
          "is finite there, or pass expand_slacks=True to fall back to the reference "
          "expansion.");
      }
    }
  }

  // Which column relaxes which side of which row; ipmc's numbering of the columns

  void IpmcInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    casadi_int struct_cnt=0;

    // Default options
    StructureDetection structure_detection = STRUCTURE_NONE;
    bool debug = false;
    nxc_ = 0;
    bool nxc_given = false;

    calc_g_ = true;
    calc_f_ = true;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="N") {
        N_ = op.second;
        struct_cnt++;
      } else if (op.first=="nx") {
        nxs_ = op.second;
        struct_cnt++;
      } else if (op.first=="nu") {
        nus_ = op.second;
        struct_cnt++;
      } else if (op.first=="ng") {
        ngs_ = op.second;
        struct_cnt++;
      } else if (op.first=="nxc") {
        // Not in struct_cnt: also valid with 'auto'
        nxc_ = op.second;
        nxc_given = true;
      } else if (op.first=="ipmc") {
        opts_ = op.second;
      } else if (op.first=="structure_detection") {
        std::string v = op.second;
        if (v=="auto") {
          structure_detection = STRUCTURE_AUTO;
        } else if (v=="manual") {
          structure_detection = STRUCTURE_MANUAL;
        } else if (v=="none") {
          structure_detection = STRUCTURE_NONE;
        } else {
          casadi_error("Unknown option for structure_detection: '" + v + "'.");
        }
      } else if (op.first=="debug") {
        debug = op.second;
      }
    }

    // Do we need second order derivatives?
    exact_hessian_ = true;
    auto hessian_approximation = opts_.find("hessian_approximation");
    if (hessian_approximation!=opts_.end()) {
      exact_hessian_ = hessian_approximation->second == "exact";
    }

    create_ipmc_functions();

    settle_slack_penalty();

    // Keep list of erroring rows
    std::set<casadi_int> errors;

    if (struct_cnt>0) {
      casadi_assert(structure_detection == STRUCTURE_MANUAL,
        "You must set structure_detection to 'manual' if you set N, nx, nu, ng.");
    }

    if (structure_detection==STRUCTURE_MANUAL) {
      casadi_assert(struct_cnt==4,
        "You must set all of N, nx, nu, ng.");
    } else if (structure_detection==STRUCTURE_NONE) {
      N_ = 0;
      nxs_ = {0};
      nus_ = {nx_};
      ngs_ = {ng_};
    } else if (structure_detection==STRUCTURE_AUTO) {
      detect_structure(errors);
    }

    casadi_assert(nxs_.size()==N_+1, "nx must have length N+1.");
    casadi_assert(nus_.size()==N_+1, "nu must have length N+1.");
    casadi_assert(ngs_.size()==N_+1, "ng must have length N+1.");
    // Checked here rather than by ipmc, to name the offending stage
    if (nxc_given) {
      casadi_assert(structure_detection!=STRUCTURE_NONE,
        "nxc needs structure_detection 'manual' or 'auto': with 'none' there is "
        "no horizon and no state to be constant.");
      for (casadi_int k=0;k<N_+1;++k) {
        casadi_assert(nxc_>=0 && nxc_<=nxs_[k],
          "nxc = " + str(nxc_) + " is not in [0, nx[" + str(k) +
          "]] = [0, " + str(nxs_[k]) + "].");
      }
    }

    // From here on, the partition of the handed-over problem
    build();
    const std::vector<casadi_int>& nx = nxh_;
    const std::vector<casadi_int>& nu = nus_;
    const std::vector<casadi_int>& ng = ngh_;

    Sparsity A_ = jac_codes().sparsity();
    if (debug) {
      A_.to_file("debug_ipmc_actual.mtx");
    }

    if (verbose_) {
      casadi_message("Using structure: N " + str(N_) + ", nx " + str(nx) + ", "
            "nu " + str(nu) + ", ng " + str(ng) + ".");
    }

    Sparsity ABsp = blocksparsity(nat_, nxt_, AB_blocks_);
    Sparsity CDsp = blocksparsity(nat_, nxt_, CD_blocks_);
    Sparsity Isp = identity_sparsity();

    Sparsity total = ABsp + CDsp + Isp;

    if (debug) {
      // A, B, C, D for debugging only
      std::vector< casadi_ocp_block > A_blocks, B_blocks, C_blocks, D_blocks;
      for (casadi_int k=0;k<=N_;++k) {
        const casadi_ocp_block& cd = CD_blocks_[k];
        if (k<N_) {
          const casadi_ocp_block& ab = AB_blocks_[k];
          A_blocks.push_back({ab.offset_r, ab.offset_c, ab.rows, nx[k]});
          B_blocks.push_back({ab.offset_r, ab.offset_c+nx[k], ab.rows, nu[k]});
        }
        C_blocks.push_back({cd.offset_r, cd.offset_c, cd.rows, nx[k]});
        D_blocks.push_back({cd.offset_r, cd.offset_c+nx[k], cd.rows, nu[k]});
      }
      total.to_file("debug_ipmc_expected.mtx");
      blocksparsity(nat_, nxt_, A_blocks).to_file("debug_ipmc_A.mtx");
      blocksparsity(nat_, nxt_, B_blocks).to_file("debug_ipmc_B.mtx");
      blocksparsity(nat_, nxt_, C_blocks).to_file("debug_ipmc_C.mtx");
      blocksparsity(nat_, nxt_, D_blocks).to_file("debug_ipmc_D.mtx");
      Isp.to_file("debug_ipmc_I.mtx");
      std::vector<casadi_int> errors_vec(errors.begin(), errors.end());
      std::vector<casadi_int> colind = {0, static_cast<casadi_int>(errors_vec.size())};
      Sparsity(nat_, 1, colind, errors_vec).to_file("debug_ipmc_errors.mtx");
    }

    casadi_assert(errors.empty() && (A_ + total).nnz() == total.nnz(),
      "Ipmc: specified structure of A does not correspond to what the interface can handle. "
      "Structure is: N " + str(N_) + ", nx " + str(nx) + ", nu " + str(nu) + ", "
      "ng " + str(ng) + ".\n"
      "Note that debug_ipmc_expected.mtx and debug_ipmc_actual.mtx are written "
      "to the current directory when 'debug' option is true.\n"
      "These can be read with Sparsity.from_file(...)."
      "For a ready-to-use script, "
      "see https://gist.github.com/jgillis/dec56fa16c90a8e4a69465e8422c5459");
    casadi_assert_dev(total.nnz() == ABsp.nnz() + CDsp.nnz() + Isp.nnz());

    build_rows();
    build_pack_tables();
    set_ipmc_prob();

    // Allocate memory
    casadi_int sz_arg, sz_res, sz_w, sz_iw;
    casadi_ipmc_work(&p_, &sz_arg, &sz_res, &sz_iw, &sz_w);
    if (n_lift_>0) {
      // the lifted problem's own z/lbz/ubz/lam
      casadi_int a2, r2, i2, w2;
      casadi_nlpsol_work(&p_nlp_lift_, &a2, &r2, &i2, &w2);
      sz_arg += a2; sz_res += r2; sz_iw += i2; sz_w += w2;
    }

    alloc_arg(sz_arg, true);
    alloc_res(sz_res, true);
    alloc_iw(sz_iw, true);
    alloc_w(sz_w, true);
  }

  void IpmcInterface::push_options(IpmcSolver* solver) const {
    for (const auto& kv : opts_) {
      switch (ipmc_option_type(kv.first.c_str())) {
        case 0:
          ipmc_set_option_double(solver, kv.first.c_str(), kv.second);
          break;
        case 1:
          ipmc_set_option_int(solver, kv.first.c_str(), kv.second.to_int());
          break;
        case 2:
          ipmc_set_option_bool(solver, kv.first.c_str(), kv.second.to_bool());
          break;
        case 3:
          {
            std::string s = kv.second.to_string();
            ipmc_set_option_string(solver, kv.first.c_str(), s.c_str());
          }
          break;
        case -1:
          casadi_error("Unknown ipmc option '" + kv.first + "'.");
        default:
          casadi_error("Unknown type of ipmc option '" + kv.first + "'.");
      }
    }
  }

  int IpmcInterface::init_mem(void* mem) const {
    if (Nlpsol::init_mem(mem)) return 1;
    auto m = static_cast<IpmcMemory*>(mem);
    casadi_ipmc_init_mem(&m->d);
    // The solver lives in m->ipmc_block for as long as the memory object does
    m->ipmc_block.assign((ipmc_memsize_ + sizeof(double) - 1)/sizeof(double), 0.0);
    IpmcError err = IPMC_OK;
    m->d.solver = ipmc_create_in(&desc_, get_ptr(m->ipmc_block),
      m->ipmc_block.size()*sizeof(double), &err);
    casadi_assert(m->d.solver, "Ipmc: the solver refused the problem description: " +
      ipmc_soft_error_message(err) + " (IpmcError == " + str(static_cast<int>(err)) + ").");
    m->d.layout = ipmc_get_layout(m->d.solver);
    ipmc_set_output(m->d.solver, &casadi_c_logger_write, &casadi_c_logger_flush);
    push_options(m->d.solver);
    return 0;
  }

  void IpmcInterface::set_work(void* mem, const double**& arg, double**& res,
      casadi_int*& iw, double*& w) const {
    auto m = static_cast<IpmcMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    m->d.prob = &p_;
    m->d.nlp = &m->d_nlp;
    m->d.rewrite.user = &m->d_nlp;

    // The lifted problem's own z/lbz/ubz/lam; set_ipmc_prob(CodeGenerator&) mirrors it
    if (n_lift_>0) {
      m->d_nlp_lift = m->d_nlp;
      m->d_nlp_lift.prob = &p_nlp_lift_;
      casadi_nlpsol_set_work(&m->d_nlp_lift, &arg, &res, &iw, &w);
      m->d.nlp = &m->d_nlp_lift;
    }

    casadi_ipmc_set_work(&m->d, &arg, &res, &iw, &w);

    // Slacks live in NlpsolMemory's slack_* scratch (never null), not in z
    if (slacks_) {
      m->d.slack_s = m->slack_s;
      m->d.slack_lam_s = m->slack_lam_s;
      m->d.slack_ubs = m->slack_ubs;
    } else {
      m->d.slack_s = nullptr;
      m->d.slack_lam_s = nullptr;
      m->d.slack_ubs = nullptr;
    }
  }

  static std::string ipmc_bounds(double lo, double up) {
    std::stringstream ss;
    ss << "lower bound " << lo << ", upper bound " << up;
    return ss.str();
  }

  std::string IpmcInterface::caller_entry(casadi_int z) const {
    return "constraint row g[" + str(lift_src_[z]-nx_+GlobalOptions::start_index) + "]";
  }

  int IpmcInterface::solve(void* mem) const {
    auto m = static_cast<IpmcMemory*>(mem);
    casadi_ipmc_data<double>* d = &m->d;

    // Caller's bounds and x0 onto the lifted problem; no-op without lift
    casadi_ipmc_rewrite_expand(&p_.rewrite, &d->rewrite, d->nlp, d->slack_ubs, d->slack_s);
    if (casadi_ipmc_hand_over(d)) {
      casadi_int i = d->error_index, si = GlobalOptions::start_index;
      const double* lbz = d->nlp->lbz;
      const double* ubz = d->nlp->ubz;
      switch (d->error) {
        case CASADI_IPMC_GAP_BOUNDS:
          casadi_error("Ipmc: " + caller_entry(i) + " closes a gap of the dynamics and needs "
            "equal finite bounds (" + ipmc_bounds(lbz[i], ubz[i]) + ").");
        case CASADI_IPMC_UBS:
          casadi_error("Ipmc: slack column " + str(i+si) + " has the upper bound ubs = " +
            str(d->slack_ubs[i]) + ", but relaxes a finite bound. Its rows sit in one stage, "
            "where ipmc carries the slack itself and needs a positive bound (inf for none); "
            "to make the rows hard, leave the column out of 'S'. A column spanning several "
            "stages may take ubs = 0.");
        default:  // CASADI_IPMC_PENALTY
          casadi_error("Ipmc: the solver refused the slack penalty: " +
            ipmc_soft_error_message(static_cast<int>(i)) + " (IpmcError == " + str(i) + ").");
      }
    }

    // casadi_ipmc_set_x fills the oracle inputs through d->arg
    casadi_assert_dev(d->arg==m->arg);

    // The solve loop.  codegen_body() emits a copy of it; keep them in step.
    {
      const IpmcLayout* str = d->layout;
      IpmcEval e;
      casadi_int i, n_ux;
      int stop = 0;

      d->unified_return_status = SOLVER_RET_UNKNOWN;
      d->success = 0;

      // stop: set by the iteration callback
      ipmc_start(d->solver);
      while (!stop && ipmc_step(d->solver, &e)!=IPMC_DONE) {
        switch (e.request) {
        // Point to the caller's, oracle, output back; without lift, lift_* are no-ops
        case IPMC_EVAL_OBJ:
          casadi_ipmc_set_x(d, e.x);
          m->res[0] = e.obj;
          calc_function(m, "nlp_f");
          casadi_ipmc_lift_obj(&p_, d, e.obj);
          *e.obj *= e.obj_scale;
          break;
        case IPMC_EVAL_OBJ_GRAD:
          casadi_ipmc_set_x(d, e.x);
          m->res[0] = d->gu;
          calc_function(m, "nlp_grad_f");
          casadi_ipmc_lift_grad_f(&p_, d);
          casadi_ipmc_write_primal_data(&p_, d->g, e.grad, str);
          casadi_scal(p_.nlp->nx, e.obj_scale, e.grad);
          break;
        case IPMC_EVAL_CONSTR_VIOL:
          casadi_ipmc_set_x(d, e.x);
          m->res[0] = d->gu;
          calc_function(m, "nlp_g");
          casadi_ipmc_lift_g(&p_, d);
          if (!fcallback_.is_null()) casadi_ipmc_snapshot_g(d, str, e.x);
          casadi_ipmc_pack_constr_viol(d, str, e.cv);
          break;
        case IPMC_EVAL_CONSTR_JAC:
          casadi_ipmc_set_x(d, e.x);
          m->res[0] = d->gu;
          m->res[1] = d->au;
          calc_function(m, "nlp_jac_g");
          casadi_ipmc_lift_g(&p_, d);
          if (casadi_ipmc_pack_constr_jac(d, str, e.BAt, e.b, e.Gt_ineq, e.g_ineq)) {
            casadi_int i = d->error_index;
            if (d->error==CASADI_IPMC_IDENTITY) casadi_error("Structure mismatch: gap-closing "
              "constraints must be like this: x_{k+1}-F(xk,uk).");
            casadi_error("nxc: stage " + casadi::str(cchk_[4*i+2]+1) + " state "
              + casadi::str(cchk_[4*i+3])
              + " was declared constant, but its dynamics row is not x_{k+1}=x_k.");
          }
          break;
        case IPMC_EVAL_LAG_HESS:
          // lam is an input to nlp_hess_l
          casadi_ipmc_set_x(d, e.x);
          casadi_ipmc_read_lam(d, str, e.lam);
          casadi_ipmc_lift_lam(&p_, d);
          m->arg[2] = &e.obj_scale;
          m->arg[3] = d->lamu;
          m->res[0] = d->gu;
          m->res[1] = d->hu;
          calc_function(m, "nlp_hess_l");
          casadi_ipmc_lift_hess_l(&p_, d, e.obj_scale);
          casadi_ipmc_pack_lag_hess(d, str, e.lam, e.obj_scale, e.RSQ, e.RSQ_slack, e.rq);
          break;
        // Iteration callback, VM only; g from the snapshot taken at IPMC_EVAL_CONSTR_VIOL
        case IPMC_POST_ITERATION:
          if (fcallback_.is_null()) break;
          // Retake the snapshot if ipmc stepped back to an earlier iterate
          n_ux = str->ux_offs[str->K-1] + str->nu[str->K-1] + str->nx[str->K-1];
          for (i=0;i<n_ux;++i) if (d->cb_x[i]!=e.x[i]) break;
          if (i<n_ux) {
            casadi_ipmc_set_x(d, e.x);
            m->res[0] = d->gu;
            calc_function(m, "nlp_g");
            casadi_ipmc_lift_g(&p_, d);
            casadi_ipmc_snapshot_g(d, str, e.x);
          }
          // A restoration iterate reports the last ordinary iterate's duals
          casadi_ipmc_report_iterate(d, str, e.x, e.lam, e.objective/e.obj_scale);
          // Throwing out of the loop is safe: ipmc_start resets all state
          if (callback(m)) {
            ipmc_abort(d->solver);
            stop = 1;
          }
          break;
        // Unreachable; keeps -Wswitch quiet
        case IPMC_DONE:
          break;
        }
      }
      casadi_ipmc_finish(d);
    }

    casadi_ipmc_rewrite_collect(&p_.rewrite, &d->rewrite, d->nlp, d->slack_s, d->slack_lam_s);

    m->success = d->success;
    m->unified_return_status = static_cast<UnifiedReturnStatus>(d->unified_return_status);

    if (m->sink) {
      casadi_ipmc_stats_set_outcome(m->sink, m->call, d->return_status,
        d->unified_return_status, d->success, d->stats.iterations_count);
    }

    return 0;
  }

  Dict IpmcInterface::get_stats(void* mem) const {
    Dict stats = Nlpsol::get_stats(mem);
    auto m = static_cast<IpmcMemory*>(mem);
    Dict ipmc;
    ipmc["compute_sd_time"] = m->d.stats.compute_sd_time;
    ipmc["duinf_time"] = m->d.stats.duinf_time;
    ipmc["eval_hess_time"] = m->d.stats.eval_hess_time;
    ipmc["eval_jac_time"] = m->d.stats.eval_jac_time;
    ipmc["eval_cv_time"] = m->d.stats.eval_cv_time;
    ipmc["eval_grad_time"] = m->d.stats.eval_grad_time;
    ipmc["eval_obj_time"] = m->d.stats.eval_obj_time;
    ipmc["initialization_time"] = m->d.stats.initialization_time;
    ipmc["time_total"] = m->d.stats.time_total;
    ipmc["eval_hess_count"] = m->d.stats.eval_hess_count;
    ipmc["eval_jac_count"] = m->d.stats.eval_jac_count;
    ipmc["eval_cv_count"] = m->d.stats.eval_cv_count;
    ipmc["eval_grad_count"] = m->d.stats.eval_grad_count;
    ipmc["eval_obj_count"] = m->d.stats.eval_obj_count;
    ipmc["iterations_count"] = m->d.stats.iterations_count;
    ipmc["restoration_iterations_count"] = m->d.stats.restoration_iterations_count;
    ipmc["return_flag"] = m->d.stats.return_flag;
    stats["ipmc"] = ipmc;
    stats["iter_count"] = m->d.stats.iterations_count;
    stats["nx"] = nxs_;
    stats["nu"] = nus_;
    stats["ng"] = ngs_;
    stats["nxc"] = nxc_;
    // Helper states the lift added, one per cross-stage slack column
    stats["n_lift"] = n_lift_;
    stats["N"] = N_;
    stats["return_status"] = casadi_ipmc_return_status_string(m->d.return_status);
    return stats;
  }

  // File-scope static of the generated function, see codegen_declarations
  static std::string ipmc_static(CodeGenerator& g, const std::string& fname,
      const std::string& name) {
    return g.shorthand(fname + "_ipmc_" + name);
  }

  // The file-scope array of a description field, 0 if empty
  static std::string ipmc_array_ref(CodeGenerator& g, const std::string& fname,
      const std::string& name, const std::vector<ipmc_int>& v) {
    return v.empty() ? std::string("0") : ipmc_static(g, fname, name);
  }

  // Define the file-scope array of a description field, nothing if empty
  static void ipmc_array_def(CodeGenerator& g, const std::string& fname,
      const std::string& name, const std::vector<ipmc_int>& v) {
    if (v.empty()) return;
    g << "static const ipmc_int " << ipmc_static(g, fname, name) << "[" << v.size() << "] = {";
    for (size_t i=0;i<v.size();++i) g << (i ? ", " : "") << v[i];
    g << "};\n";
  }

  std::string IpmcInterface::codegen_desc(CodeGenerator& g) const {
    std::string fname = codegen_name(g, false);
    g.local("desc", "struct IpmcProblem");
    // Fields not set below are zero
    g << "desc = " << ipmc_static(g, fname, "desc0") << ";\n";
    g << "desc.K = " << desc_.K << ";\n";
    g << "desc.nu = " << ipmc_array_ref(g, fname, "nu", d_nu_) << ";\n";
    g << "desc.nx = " << ipmc_array_ref(g, fname, "nx", d_nx_) << ";\n";
    g << "desc.ng = " << ipmc_array_ref(g, fname, "ng", d_ng_) << ";\n";
    g << "desc.ng_ineq = " << ipmc_array_ref(g, fname, "ng_ineq", d_ng_ineq_) << ";\n";
    g << "desc.nxc = " << desc_.nxc << ";\n";
    if (desc_.ns) {
      g << "desc.ns = " << ipmc_array_ref(g, fname, "ns", d_ns_) << ";\n";
      g << "desc.soft_lo = " << ipmc_array_ref(g, fname, "soft_lo", d_soft_lo_) << ";\n";
      g << "desc.soft_up = " << ipmc_array_ref(g, fname, "soft_up", d_soft_up_) << ";\n";
    }
    if (desc_.slack_helper) {
      g << "desc.nxc_slack = " << desc_.nxc_slack << ";\n";
      g << "desc.slack_helper = " << ipmc_array_ref(g, fname, "slack_helper", d_slack_helper_)
        << ";\n";
    }
    g << "desc.nb = " << ipmc_array_ref(g, fname, "nb", d_nb_) << ";\n";
    g << "desc.idxb = " << ipmc_array_ref(g, fname, "idxb", d_idxb_) << ";\n";
    return "desc";
  }

  // The solver, once per memory object, in a static block sized at codegen time
  void IpmcInterface::codegen_init_mem(CodeGenerator& g) const {
    std::string fname = codegen_name(g, false);
    std::string block = ipmc_static(g, fname, "block");
    std::string mem = codegen_mem(g);
    std::string desc = codegen_desc(g);
    g.local("err", "IpmcError");
    g << "casadi_ipmc_init_mem(&" << mem << ");\n";
    g << mem << ".solver = ipmc_create_in(&" << desc << ", " << block << "[mem], sizeof("
      << block << "[0]), &err);\n";
    g << "if (!" << mem << ".solver) {\n"
      << "CASADI_PRINTF(\"Ipmc: the solver refused the problem description (IpmcError == %d; "
      << "8: the block sized when this code was generated is too small for the ipmc linked "
      << "in).\\n\", "
      << "(int) err);\n"
      << "return 1;\n"
      << "}\n";
    g << mem << ".layout = ipmc_get_layout(" << mem << ".solver);\n";
    g << "ipmc_set_output(" << mem << ".solver, &" << ipmc_static(g, fname, "cb_write") << ", &"
      << ipmc_static(g, fname, "cb_flush") << ");\n";
    for (const auto& kv : opts_) {
      std::string call = "ipmc_set_option_";
      switch (ipmc_option_type(kv.first.c_str())) {
        case 0:
          g << call << "double(" << mem << ".solver, \"" << kv.first << "\", "
            << g.constant(kv.second.to_double()) << ");\n";
          break;
        case 1:
          g << call << "int(" << mem << ".solver, \"" << kv.first << "\", "
            << kv.second.to_int() << ");\n";
          break;
        case 2:
          g << call << "bool(" << mem << ".solver, \"" << kv.first << "\", "
            << static_cast<int>(kv.second.to_bool()) << ");\n";
          break;
        case 3:
          g << call << "string(" << mem << ".solver, \"" << kv.first << "\", \""
            << kv.second.to_string() << "\");\n";
          break;
        case -1:
          casadi_error("Unknown ipmc option '" + kv.first + "'.");
        default:
          casadi_error("Unknown type of ipmc option '" + kv.first + "'.");
      }
    }
    g << "return 0;\n";
  }

  void IpmcInterface::codegen_declarations(CodeGenerator& g) const {
    Nlpsol::codegen_declarations(g);
    g.add_auxiliary(CodeGenerator::AUX_NLP);
    g.add_auxiliary(CodeGenerator::AUX_INF);
    g.add_auxiliary(CodeGenerator::AUX_MAX);
    g.add_auxiliary(CodeGenerator::AUX_COPY);
    g.add_auxiliary(CodeGenerator::AUX_SCAL);
    g.add_auxiliary(CodeGenerator::AUX_OCP_BLOCK);
    g.add_auxiliary(CodeGenerator::AUX_SCALED_COPY);
    g.add_auxiliary(CodeGenerator::AUX_AXPY);
    g.add_auxiliary(CodeGenerator::AUX_CLEAR);
    g.add_auxiliary(CodeGenerator::AUX_MV);
    g.add_auxiliary(CodeGenerator::AUX_PRINTF);
    g.add_dependency(get_function("nlp_f"));
    g.add_dependency(get_function("nlp_grad_f"));
    g.add_dependency(get_function("nlp_g"));
    g.add_dependency(get_function("nlp_jac_g"));
    if (exact_hessian_) {
      g.add_dependency(get_function("nlp_hess_l"));
    }
    g.add_include("ipmc/ipmc.h");

    // ipmc's output, per generated function
    std::string fname = codegen_name(g, false);
    g << "void " << ipmc_static(g, fname, "cb_write") << "(const char* msg, int num) {\n";
    g.flush(g.body);
    g.scope_enter();
    g << "CASADI_PRINTF(\"%.*s\", num, msg);\n";
    g.scope_exit();
    g << "}\n";

    g << "void " << ipmc_static(g, fname, "cb_flush") << "(void) {\n";
    g.flush(g.body);
    g.scope_enter();
    g.scope_exit();
    g << "}\n";

    // The IpmcProblem arrays, and the block every memory object builds its solver in
    g << "static const struct IpmcProblem " << ipmc_static(g, fname, "desc0") << ";\n";
    ipmc_array_def(g, fname, "nu", d_nu_);
    ipmc_array_def(g, fname, "nx", d_nx_);
    ipmc_array_def(g, fname, "ng", d_ng_);
    ipmc_array_def(g, fname, "ng_ineq", d_ng_ineq_);
    if (desc_.ns) {
      ipmc_array_def(g, fname, "ns", d_ns_);
      ipmc_array_def(g, fname, "soft_lo", d_soft_lo_);
      ipmc_array_def(g, fname, "soft_up", d_soft_up_);
    }
    if (desc_.slack_helper) ipmc_array_def(g, fname, "slack_helper", d_slack_helper_);
    ipmc_array_def(g, fname, "nb", d_nb_);
    ipmc_array_def(g, fname, "idxb", d_idxb_);
    g << "static double " << ipmc_static(g, fname, "block") << "[CASADI_MAX_NUM_THREADS]["
      << (ipmc_memsize_ + sizeof(double) - 1)/sizeof(double) << "];\n";
    g.flush(g.body);
  }

  void IpmcInterface::codegen_body(CodeGenerator& g) const {
    codegen_body_enter(g);
    // The runtime, once per generated file
    if (g.auxiliaries.str().find("struct casadi_ipmc_data {")==std::string::npos) {
      g.auxiliaries << g.sanitize_source(ipmc_runtime_str, {"casadi_real"});
    }
    if (g.stats()
        && g.auxiliaries.str().find("casadi_ipmc_stats_set_outcome(")==std::string::npos) {
      g.auxiliaries << g.sanitize_source(ipmc_base_runtime_str, {"casadi_real"});
    }

    g.local("d", "struct casadi_ipmc_data*");
    g.init_local("d", "&" + codegen_mem(g));
    g.local("p", "struct casadi_ipmc_prob");
    set_ipmc_prob(g);

    g << "casadi_ipmc_set_work(d, &arg, &res, &iw, &w);\n";

    // Mirror of set_work(); locals come from Nlpsol::codegen_slack_native_enter
    if (slacks_) {
      g << "d->slack_s = nlp_slack_s;\n";
      g << "d->slack_lam_s = nlp_slack_lam_s;\n";
      g << "d->slack_ubs = nlp_slack_ubs;\n";
    } else {
      g << "d->slack_s = 0;\n";
      g << "d->slack_lam_s = 0;\n";
      g << "d->slack_ubs = 0;\n";
    }

    // See solve(); the refusals name the caller's entry, as there
    g << "casadi_ipmc_rewrite_expand(&p.rewrite, &d->rewrite, d->nlp, "
         "d->slack_ubs, d->slack_s);\n";
    std::string caller = g.constant(lift_src_);
    std::string nx = str(nx_), si = str(GlobalOptions::start_index);
    g << "if (casadi_ipmc_hand_over(d)) {\n";
    // CASADI_IPMC_UBS needs stage-local slacks
    if (!slack_perm_.empty()) g << "if (d->error==CASADI_IPMC_UBS) {\n"
      << "CASADI_PRINTF(\"Ipmc: slack column %d has an upper bound ubs <= 0, but relaxes a "
         "finite bound. Its rows sit in one stage, where ipmc carries the slack itself and "
         "needs a positive bound (inf for none); to make the rows hard, leave the column out "
         "of 'S'. A column spanning several stages may take ubs = 0.\\n\", "
         "(int) d->error_index + " << si << ");\n"
      << "} else ";
    g << "if (d->error==CASADI_IPMC_PENALTY) {\n"
      << "CASADI_PRINTF(\"Ipmc: the solver refused the slack penalty (IpmcError == %d).\\n\", "
         "(int) d->error_index);\n"
      << "} else {\n"
      << "CASADI_PRINTF(\"Ipmc: constraint row g[%d] closes a gap of the dynamics and needs "
         "equal finite bounds.\\n\", (int) (" << caller << "[d->error_index]-" << nx << ") + "
      << si << ");\n"
      << "}\n"
      << "return 1;\n"
      << "}\n";
    // Copy of the solve loop in solve(); keep them in step
    g.local("str", "const struct IpmcLayout", "*");
    g.local("e", "struct IpmcEval");
    g << "str = d->layout;\n";
    g << "d->unified_return_status = " << static_cast<int>(SOLVER_RET_UNKNOWN) << ";\n";
    g << "d->success = 0;\n";
    g << "ipmc_start(d->solver);\n";
    g << "while (ipmc_step(d->solver, &e)!=IPMC_DONE) {\n";
    g << "switch (e.request) {\n";

    g << "case IPMC_EVAL_OBJ:\n";
    g << "casadi_ipmc_set_x(d, e.x);\n";
    g << "d->res[0] = e.obj;\n";
    g << g(get_function("nlp_f"), "d->arg", "d->res", "d->iw", "d->w") << ";\n";
    g << "casadi_ipmc_lift_obj(&p, d, e.obj);\n";
    g << "*e.obj *= e.obj_scale;\n";
    g << "break;\n";

    g << "case IPMC_EVAL_OBJ_GRAD:\n";
    g << "casadi_ipmc_set_x(d, e.x);\n";
    g << "d->res[0] = d->gu;\n";
    g << g(get_function("nlp_grad_f"), "d->arg", "d->res", "d->iw", "d->w") << ";\n";
    g << "casadi_ipmc_lift_grad_f(&p, d);\n";
    g << "casadi_ipmc_write_primal_data(&p, d->g, e.grad, str);\n";
    g << "casadi_scal(" << nxt_ << ", e.obj_scale, e.grad);\n";
    g << "break;\n";

    g << "case IPMC_EVAL_CONSTR_VIOL:\n";
    g << "casadi_ipmc_set_x(d, e.x);\n";
    g << "d->res[0] = d->gu;\n";
    g << g(get_function("nlp_g"), "d->arg", "d->res", "d->iw", "d->w") << ";\n";
    g << "casadi_ipmc_lift_g(&p, d);\n";
    g << "casadi_ipmc_pack_constr_viol(d, str, e.cv);\n";
    g << "break;\n";

    g << "case IPMC_EVAL_CONSTR_JAC:\n";
    g << "casadi_ipmc_set_x(d, e.x);\n";
    g << "d->res[0] = d->gu;\n";
    g << "d->res[1] = d->au;\n";
    g << g(get_function("nlp_jac_g"), "d->arg", "d->res", "d->iw", "d->w") << ";\n";
    g << "casadi_ipmc_lift_g(&p, d);\n";
    g << "if (casadi_ipmc_pack_constr_jac(d, str, e.BAt, e.b, e.Gt_ineq, e.g_ineq)) {\n";
    g << "if (d->error==CASADI_IPMC_IDENTITY) {\n"
      << "CASADI_PRINTF(\"Structure mismatch: gap-closing constraints must be like this: "
         "x_{k+1}-F(xk,uk).\\n\");\n"
      << "} else {\n"
      << "CASADI_PRINTF(\"nxc: stage %d state %d was declared constant, but its dynamics row is "
         "not x_{k+1}=x_k.\\n\", (int) p.cchk[4*d->error_index+2]+1, "
         "(int) p.cchk[4*d->error_index+3]);\n"
      << "}\n";
    g << "return 1;\n";
    g << "}\n";
    g << "break;\n";

    g << "case IPMC_EVAL_LAG_HESS:\n";
    g << "casadi_ipmc_set_x(d, e.x);\n";
    g << "casadi_ipmc_read_lam(d, str, e.lam);\n";
    g << "casadi_ipmc_lift_lam(&p, d);\n";
    g << "d->arg[2] = &e.obj_scale;\n";
    g << "d->arg[3] = d->lamu;\n";
    g << "d->res[0] = d->gu;\n";
    g << "d->res[1] = d->hu;\n";
    g << g(get_function("nlp_hess_l"), "d->arg", "d->res", "d->iw", "d->w") << ";\n";
    g << "casadi_ipmc_lift_hess_l(&p, d, e.obj_scale);\n";
    g << "casadi_ipmc_pack_lag_hess(d, str, e.lam, e.obj_scale, e.RSQ, e.RSQ_slack, e.rq);\n";
    g << "break;\n";

    // No iteration callback in generated C
    g << "case IPMC_POST_ITERATION:\n";
    g << "case IPMC_DONE:\n";
    g << "break;\n";

    g << "}\n";
    g << "}\n";
    g << "casadi_ipmc_finish(d);\n";
    if (g.stats()) {
      g << "casadi_ipmc_stats_set_outcome(sink, call, d->return_status, "
        << "d->unified_return_status, d->success, d->stats.iterations_count);\n";
    }
    g << "casadi_ipmc_rewrite_collect(&p.rewrite, &d->rewrite, d->nlp, "
         "d->slack_s, d->slack_lam_s);\n";

    codegen_body_exit(g);

    if (error_on_fail_) {
      g << "return d->unified_return_status;\n";
    } else {
      g << "return 0;\n";
    }
  }

  static std::vector<casadi_int> ipmc_blocks_pack(const std::vector<casadi_ocp_block>& blocks) {
    size_t N = blocks.size();
    std::vector<casadi_int> ret(4*N+1);
    casadi_int* r = get_ptr(ret);
    *r++ = N;
    for (casadi_int i=0;i<N;++i) {
      *r++ = blocks[i].offset_r;
      *r++ = blocks[i].offset_c;
      *r++ = blocks[i].rows;
      *r++ = blocks[i].cols;
    }
    return ret;
  }

  // Every constant of casadi_ipmc_prob as (field, C name, value); an empty vector is null
  // Every constant of casadi_ipmc_prob as (field, C name, value); an empty vector is null
  template<typename V>
  void IpmcInterface::prob_fields(casadi_ipmc_prob<double>& p, V& v) const {
    v(p.nx, "nx", nxh_);
    v(p.nu, "nu", nus_);
    v(p.N, "N", N_);

    v(p.n_soft, "n_soft", static_cast<casadi_int>(slack_perm_.size()));
    v(p.slack_perm, "slack_perm", slack_perm_);
    v(p.fs_z, "fs_z", fs_z_);
    v(p.fs_Z, "fs_Z", fs_Z_);

    v(p.rewrite.ne, "rewrite.ne", n_lift_);
    v(p.rewrite.nx, "rewrite.nx", nx_);
    v(p.rewrite.ng, "rewrite.ng", ng_);
    v(p.rewrite.m0, "rewrite.m0", lift_m0_);
    v(p.rewrite.col, "rewrite.col", lift_col_);
    v(p.rewrite.Px_sp, "rewrite.Px_sp", lift_Px_.sparsity());
    v(p.rewrite.Px, "rewrite.Px", lift_Px_.nonzeros());
    v(p.rewrite.Pg_sp, "rewrite.Pg_sp", lift_Pg_.sparsity());
    v(p.rewrite.Pg, "rewrite.Pg", lift_Pg_.nonzeros());
    v(p.rewrite.C_sp, "rewrite.C_sp", lift_C_.sparsity());
    v(p.rewrite.C, "rewrite.C", lift_C_.nonzeros());
    v(p.rewrite.B_sp, "rewrite.B_sp", lift_B_.sparsity());
    v(p.rewrite.B, "rewrite.B", lift_B_.nonzeros());
    v(p.rewrite.Blo_sp, "rewrite.Blo_sp", lift_Blo_.sparsity());
    v(p.rewrite.Blo, "rewrite.Blo", lift_Blo_.nonzeros());
    v(p.rewrite.Bup_sp, "rewrite.Bup_sp", lift_Bup_.sparsity());
    v(p.rewrite.Bup, "rewrite.Bup", lift_Bup_.nonzeros());
    v(p.rewrite.Gz_sp, "rewrite.Gz_sp", lift_Gz_.sparsity());
    v(p.rewrite.Gz, "rewrite.Gz", lift_Gz_.nonzeros());
    v(p.rewrite.H_sp, "rewrite.H_sp", lift_H_.sparsity());
    v(p.rewrite.H, "rewrite.H", lift_H_.nonzeros());
    v(p.rewrite.lbz0, "rewrite.lbz0", lift_lbz0_);
    v(p.rewrite.ubz0, "rewrite.ubz0", lift_ubz0_);
    v(p.rewrite.nnz_au, "rewrite.nnz_au", jacg_sp_.nnz());
    v(p.rewrite.nnz_hu, "rewrite.nnz_hu", exact_hessian_ ? hesslag_sp_.nnz() : 0);

    v(p.n_ineq, "n_ineq", static_cast<casadi_int>(ineq_z_.size()));
    v(p.ineq_z, "ineq_z", ineq_z_);
    v(p.ineq_lo, "ineq_lo", ineq_lo_);
    v(p.ineq_up, "ineq_up", ineq_up_);

    v(p.bat_sp, "bat_sp", bat_sp_);
    v(p.bat_code, "bat_code", bat_code_);
    v(p.bat_blk, "bat_blk", bat_blk_);
    v(p.bat_col, "bat_col", bat_col_);
    v(p.rsq_sp, "rsq_sp", rsq_sp_);
    v(p.rsq_code, "rsq_code", rsq_code_);
    v(p.rsq_blk, "rsq_blk", rsq_blk_);
    v(p.rsq_col, "rsq_col", rsq_col_);
    v(p.rsqs_sp, "rsqs_sp", rsqs_sp_);
    v(p.rsqs_code, "rsqs_code", rsqs_code_);
    v(p.gi_sp, "gi_sp", gi_sp_);
    v(p.gi_code, "gi_code", gi_code_);
    v(p.gi_blk, "gi_blk", gi_blk_);
    v(p.gi_col, "gi_col", gi_col_);
    v(p.n_ichk, "n_ichk", static_cast<casadi_int>(ichk_.size()));
    v(p.ichk, "ichk", ichk_);
    v(p.n_cchk, "n_cchk", static_cast<casadi_int>(cchk_.size()/4));
    v(p.cchk, "cchk", cchk_);
  }

  struct IpmcPointAt {
    void operator()(casadi_int& f, const char*, casadi_int v) { f = v; }
    void operator()(const casadi_int*& f, const char*, const Sparsity& v) { f = v; }
    template<typename T>
    void operator()(const T*& f, const char*, const std::vector<T>& v) { f = get_ptr(v); }
  };

  // prob_fields() visitor of the codegen mirror: emit the value as a constant
  struct IpmcEmit {
    CodeGenerator& g;
    void operator()(casadi_int, const char* n, casadi_int v) {
      g << "p." << n << " = " << v << ";\n";
    }
    void operator()(const casadi_int*, const char* n, const Sparsity& v) {
      g << "p." << n << " = " << g.sparsity(v) << ";\n";
    }
    template<typename T>
    void operator()(const T*, const char* n, const std::vector<T>& v) {
      g << "p." << n << " = " << (v.empty() ? std::string("0") : g.constant(v)) << ";\n";
    }
  };

  void IpmcInterface::set_ipmc_prob() {
    p_nlp_lift_ = p_nlp_;
    p_nlp_lift_.nx = nxt_;
    p_nlp_lift_.ng = nat_;
    p_nlp_lift_.detect_bounds.ng = 0;
    p_.nlp = n_lift_>0 ? &p_nlp_lift_ : &p_nlp_;
    p_.AB = get_ptr(AB_blocks_);
    p_.CD = get_ptr(CD_blocks_);
    IpmcPointAt v;
    prob_fields(p_, v);

    // The IpmcProblem, and the memory a solver for it needs
    memset(&desc_, 0, sizeof(desc_));
    desc_.K = N_+1;
    desc_.nu = get_ptr(d_nu_);
    desc_.nx = get_ptr(d_nx_);
    desc_.ng = get_ptr(d_ng_);
    desc_.ng_ineq = get_ptr(d_ng_ineq_);
    desc_.nxc = nxc_ + n_lift_;
    if (!slack_perm_.empty()) {
      desc_.ns = get_ptr(d_ns_);
      desc_.soft_lo = get_ptr(d_soft_lo_);
      desc_.soft_up = get_ptr(d_soft_up_);
    }
    if (n_lift_>0) {
      desc_.nxc_slack = n_lift_;
      desc_.slack_helper = get_ptr(d_slack_helper_);
    }
    desc_.nb = get_ptr(d_nb_);
    desc_.idxb = get_ptr(d_idxb_);
    IpmcError err = IPMC_OK;
    ipmc_memsize_ = ipmc_memsize(&desc_, &err);
    casadi_assert(ipmc_memsize_>0, "Ipmc: the solver refused the problem description: " +
      ipmc_soft_error_message(err) + " (IpmcError == " + str(static_cast<int>(err)) + ").");
  }

  static void codegen_unpack_block(CodeGenerator& g, const std::string& name,
      const std::vector<casadi_ocp_block>& blocks) {
    casadi_int sz = blocks.size();
    if (sz==0) sz++;
    std::string n = "block_" + name + "[" + str(sz) + "]";
    g.local(n, "static struct casadi_ocp_block");
    g << "p." << name << " = block_" + name + ";\n";
    g << "casadi_unpack_ocp_blocks(" << "p." << name
    << ", " << g.constant(ipmc_blocks_pack(blocks)) << ");\n";
  }

  void IpmcInterface::set_ipmc_prob(CodeGenerator& g) const {
    g << "d->nlp = &d_nlp;\n";
    g << "d->rewrite.user = &d_nlp;\n";
    g << "d->prob = &p;\n";
    if (n_lift_>0) {
      // Mirror of set_work(); w consumption must match
      g.local("p_nlp_lift", "struct casadi_nlpsol_prob");
      g.local("d_nlp_lift", "struct casadi_nlpsol_data");
      g << "p_nlp_lift.nx = " << nxt_ << ";\n";
      g << "p_nlp_lift.ng = " << nat_ << ";\n";
      g << "p_nlp_lift.np = " << np_ << ";\n";
      g << "p_nlp_lift.detect_bounds.ng = 0;\n";
      // All of d_nlp, then the overrides
      g << "d_nlp_lift = d_nlp;\n";
      g << "d_nlp_lift.prob = &p_nlp_lift;\n";
      g << "casadi_nlpsol_set_work(&d_nlp_lift, &arg, &res, &iw, &w);\n";
      g << "d->nlp = &d_nlp_lift;\n";
      g << "p.nlp = &p_nlp_lift;\n";
    } else {
      g << "p.nlp = &p_nlp;\n";
    }

    codegen_unpack_block(g, "AB", AB_blocks_);
    codegen_unpack_block(g, "CD", CD_blocks_);
    IpmcEmit v{g};
    casadi_ipmc_prob<double> unused;
    prob_fields(unused, v);
  }

  IpmcInterface::IpmcInterface(DeserializingStream& s) : Nlpsol(s) {
    s.version("IpmcInterface", 7);
    s.unpack("IpmcInterface::jacg_sp", jacg_sp_);
    s.unpack("IpmcInterface::hesslag_sp", hesslag_sp_);
    s.unpack("IpmcInterface::exact_hessian", exact_hessian_);
    s.unpack("IpmcInterface::opts", opts_);

    // The caller's partition
    s.unpack("IpmcInterface::nxs", nxs_);
    s.unpack("IpmcInterface::nus", nus_);
    s.unpack("IpmcInterface::ngs", ngs_);
    s.unpack("IpmcInterface::nxc", nxc_);
    s.unpack("IpmcInterface::N", N_);

    // Native soft constraints; the maps follow from S, stored by Nlpsol
    s.unpack("IpmcInterface::slacks", slacks_);
    if (slacks_) {
      s.unpack("IpmcInterface::fs_z", fs_z_);
      s.unpack("IpmcInterface::fs_Z", fs_Z_);
    }

    // The slack maps, the lift, ipmc's rows and the pack tables are rebuilt, not stored
    build();
    build_rows();
    build_pack_tables();
    set_ipmc_prob();
  }

  void IpmcInterface::serialize_body(SerializingStream& s) const {
    Nlpsol::serialize_body(s);
    s.version("IpmcInterface", 7);

    s.pack("IpmcInterface::jacg_sp", jacg_sp_);
    s.pack("IpmcInterface::hesslag_sp", hesslag_sp_);
    s.pack("IpmcInterface::exact_hessian", exact_hessian_);
    s.pack("IpmcInterface::opts", opts_);

    s.pack("IpmcInterface::nxs", nxs_);
    s.pack("IpmcInterface::nus", nus_);
    s.pack("IpmcInterface::ngs", ngs_);
    s.pack("IpmcInterface::nxc", nxc_);
    s.pack("IpmcInterface::N", N_);

    s.pack("IpmcInterface::slacks", slacks_);
    if (slacks_) {
      s.pack("IpmcInterface::fs_z", fs_z_);
      s.pack("IpmcInterface::fs_Z", fs_Z_);
    }
  }

} // namespace casadi
