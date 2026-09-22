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
#include "../../core/global_options.hpp"
#include "../../core/casadi_interrupt.hpp"
#include "../../core/convexify.hpp"
#include "casadi/casadi_c.h"

#include <cmath>
#include <ctime>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>

#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.mutex.h>
#else // CASADI_WITH_THREAD_MINGW
#include <mutex>
#endif // CASADI_WITH_THREAD_MINGW
#endif //CASADI_WITH_THREAD

#include <ipmc_runtime_str.h>

namespace casadi {

  extern "C"
  int CASADI_NLPSOL_IPMC_EXPORT
  casadi_register_nlpsol_ipmc(Nlpsol::Plugin* plugin) {
    plugin->creator = IpmcInterface::creator;
    plugin->name = "ipmc";
    plugin->doc = IpmcInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &IpmcInterface::options_;
    plugin->deserialize = &IpmcInterface::deserialize;
    // Native soft constraints: the oracle stays (x,p,s)->(f,g,f_s), nx_==nxu_
    // and ng_==ngu_, so structure detection keeps seeing the user's own OCP.
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

  // Decodes IpmcError; 5 and 7 are ours (see casadi_ipmc_presolve), and 99
  // is the mockup's (casadi/mockups, ipmc/src/ipmc.c)
  std::string ipmc_soft_error_message(int code) {
    switch (code) {
      case 1: return "a stage reported a negative number of slack pairs";
      case 2: return "a slack index was out of range -- this is also what a slack shared "
                     "across stages looks like to ipmc";
      case 3: return "a penalty weight (z or Z) of the slack objective f_s is negative; "
                     "ipmc needs a convex, non-decreasing penalty";
      case 4: return "a slack upper bound 'ubs' is negative";
      case 5: return "'ubs' is zero on one side of a slack column and positive on the "
                     "other. ipmc needs a strictly positive slack upper bound, and a "
                     "one-sided pin cannot be expressed as a slack PAIR. Set both sides "
                     "of that column to zero (which makes the rows hard), or give the "
                     "pinned side a positive bound";
      case 7: return "out of memory while allocating the problem-structure cache";
      /* Not ipmc's: the mockup libipmc's only way of saying what it is.  Only
         reachable when a mockup ends up on the library search path, since a
         casadi built against it does not ship one -- the plugin normally just
         fails to load and 'ipmc' is reported as unavailable. */
      case 99: return "the libipmc that was loaded is the MOCKUP -- ipmc's symbols "
                      "with no solver behind them. Put a licensed libipmc on the "
                      "library search path (LD_LIBRARY_PATH, PATH or "
                      "DYLD_LIBRARY_PATH) ahead of it";
      default: return "unknown reason";
    }
  }

  void report_issue(casadi_int i, const std::string& msg) {
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
        "Number of CONSTANT states -- states whose dynamics are the identity, "
        "x_{k+1}=x_k -- counted from the END of the state vector; one number "
        "for the whole horizon.  Write the problem exactly as you would "
        "without this option; nxc only tells the interface which trailing "
        "block of A and B it may drop on the way to ipmc, which then exploits "
        "the identity in the Riccati recursion instead of factorizing it. "
        "0 <= nxc <= min_k nx[k].  Needs structure_detection 'manual' or "
        "'auto'; absent or 0 is bit-for-bit the behaviour without it."}},
      {"ipmc",
       {OT_DICT,
        "Options to be passed to ipmc"}},
      {"structure_detection",
       {OT_STRING,
        "NONE | auto | manual"}},
      {"convexify_strategy",
       {OT_STRING,
        "NONE|regularize|eigen-reflect|eigen-clip. "
        "Strategy to convexify the Lagrange Hessian before passing it to the solver."}},
      {"convexify_margin",
       {OT_DOUBLE,
        "When using a convexification strategy, make sure that "
        "the smallest eigenvalue is at least this (default: 1e-7)."}},
      {"debug",
       {OT_BOOL,
        "Produce debug information (default: false)"}}
     }
  };

  /* The five functions the runtime evaluates, of whichever problem is being
     handed over.  Factored out because the oracle is not always known at the
     same moment: only a problem that CAN be rewritten (native slacks, and a
     block structure to be cross-stage with respect to) has to wait for the
     structure detection first. */
  void IpmcInterface::create_ipmc_functions(const Function& orc) {
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

  /* ---- the lifted problem's patterns and nonzero maps -------------------
     The rewritten problem's quantities are, per rewritten row,
       g[src]                        a row that was left alone
       g[src] + m_e / g[src] - m_e   the two halves of a split row
       x[src] + m_e / x[src] - m_e   the same for a softened simple bound
       m_{k+1}[e] - m_k[e]           a helper's own dynamics, in casadi's
                                     gap-closing convention so the identity
                                     check on constant states covers it
     and f + sum_e (z_e m_0[e] + 1/2 Z_e m_0[e]^2) for the objective.  All
     of it is the caller's expression plus something LINEAR with constant
     +-1 coefficients, so the jacobian and Hessian of the lifted problem
     are the caller's nonzeros rearranged plus constants -- one map per
     pattern, built here once and applied per solve by the runtime's
     lift_jac_g / lift_hess_l.  (A symbolic lift_oracle was tried instead
     and cost more than it saved; see git history.) */
  void IpmcInterface::lift_patterns() {
    // caller g row -> the rewritten rows carrying it (one, or a split pair)
    std::vector<casadi_int> nrows(ng_, 0), rmap(2*ng_, -1);
    for (casadi_int r=0;r<nat_;++r) {
      if (lift_kind_[r]!=3 && lift_row_[r]>=nx_) {
        casadi_int j = lift_row_[r]-nx_;
        rmap[2*j+nrows[j]] = r;
        nrows[j]++;
      }
    }
    // the caller's nonzeros, into every rewritten row that carries their row
    std::vector<casadi_int> arow, acol, asrc;
    const casadi_int* colind = jacg_sp_.colind();
    const casadi_int* row = jacg_sp_.row();
    for (casadi_int c=0;c<nx_;++c) {
      for (casadi_int el=colind[c]; el<colind[c+1]; ++el) {
        casadi_int j = row[el];
        for (casadi_int t=0;t<nrows[j];++t) {
          arow.push_back(rmap[2*j+t]);
          acol.push_back(lift_x_[c]);
          asrc.push_back(el);
        }
      }
    }
    // the constant entries: -1 encodes +1.0, -2 encodes -1.0
    for (casadi_int r=0;r<nat_;++r) {
      casadi_int kind = lift_kind_[r], fl = lift_mflat_[r];
      if (kind==3) {
        arow.push_back(r); acol.push_back(lift_m_[fl]);          asrc.push_back(-2);
        arow.push_back(r); acol.push_back(lift_m_[fl+n_lift_]);  asrc.push_back(-1);
      } else if (kind==1 || kind==2) {
        if (lift_row_[r]<nx_) {
          arow.push_back(r); acol.push_back(lift_x_[lift_row_[r]]); asrc.push_back(-1);
        }
        arow.push_back(r); acol.push_back(lift_m_[fl]); asrc.push_back(kind==1 ? -1 : -2);
      }
    }
    // triplet's mapping (invert_mapping=false) is NONZERO -> triplet index
    std::vector<casadi_int> mapping;
    jacg_lift_sp_ = Sparsity::triplet(nat_, nxt_, arow, acol, mapping, false);
    // no two triplets may collide: the caller's g does not depend on m, and
    // the added rows are new
    casadi_assert_dev(jacg_lift_sp_.nnz()==static_cast<casadi_int>(arow.size()));
    lift_a_src_.assign(jacg_lift_sp_.nnz(), 0);
    for (size_t nz=0;nz<mapping.size();++nz) lift_a_src_[nz] = asrc[mapping[nz]];

    if (exact_hessian_) {
      // the caller's Hessian in its new columns, plus the stage-0 penalty
      // diagonal of every helper (-2-e encodes component e)
      std::vector<casadi_int> hrow, hcol, hsrc;
      const casadi_int* hcolind = hesslag_sp_.colind();
      const casadi_int* hr = hesslag_sp_.row();
      for (casadi_int c=0;c<nx_;++c) {
        for (casadi_int el=hcolind[c]; el<hcolind[c+1]; ++el) {
          hrow.push_back(lift_x_[hr[el]]);
          hcol.push_back(lift_x_[c]);
          hsrc.push_back(el);
        }
      }
      for (casadi_int e=0;e<n_lift_;++e) {
        hrow.push_back(lift_m_[e]); hcol.push_back(lift_m_[e]); hsrc.push_back(-2-e);
      }
      mapping.clear();
      hess_lift_sp_ = Sparsity::triplet(nxt_, nxt_, hrow, hcol, mapping, false);
      casadi_assert_dev(hess_lift_sp_.nnz()==static_cast<casadi_int>(hrow.size()));
      lift_h_src_.assign(hess_lift_sp_.nnz(), 0);
      for (size_t nz=0;nz<mapping.size();++nz) lift_h_src_[nz] = hsrc[mapping[nz]];
    }
  }

  void IpmcInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    casadi_int struct_cnt=0;

    // Default options
    std::string convexify_strategy = "none";
    double convexify_margin = 1e-7;
    casadi_int max_iter_eig = 200;
    structure_detection_ = STRUCTURE_NONE;
    debug_ = false;
    /* nxc_ is only validated when the user actually declared it: 0 is both
       the default and a legal declaration. */
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
        /* deliberately NOT counted in struct_cnt: nxc is orthogonal to who
           found the partition, and is usable with 'auto' too */
        nxc_ = op.second;
        nxc_given = true;
      } else if (op.first=="convexify_strategy") {
        convexify_strategy = op.second.to_string();
      } else if (op.first=="convexify_margin") {
        convexify_margin = op.second;
      } else if (op.first=="max_iter_eig") {
        max_iter_eig = op.second;
      } else if (op.first=="ipmc") {
        opts_ = op.second;
      } else if (op.first=="structure_detection") {
        std::string v = op.second;
        if (v=="auto") {
          structure_detection_ = STRUCTURE_AUTO;
        } else if (v=="manual") {
          structure_detection_ = STRUCTURE_MANUAL;
        } else if (v=="none") {
          structure_detection_ = STRUCTURE_NONE;
        } else {
          casadi_error("Unknown option for structure_detection: '" + v + "'.");
        }
      } else if (op.first=="debug") {
        debug_ = op.second;
      }
    }

    // Do we need second order derivatives?
    exact_hessian_ = true;
    auto hessian_approximation = opts_.find("hessian_approximation");
    if (hessian_approximation!=opts_.end()) {
      exact_hessian_ = hessian_approximation->second == "exact";
    }

    /* The functions are always the CALLER's: a lifted (cross-stage slack)
       problem reuses them too, its rewritten quantities produced numerically
       in the solve loop by the runtime's lift_* routines. */
    create_ipmc_functions(oracle_);

    /* ---- native soft constraints: settle the penalty, once ---------------
       ipmc's model is  zl*sl + 1/2 Zl sl^2 + zu*su + 1/2 Zu su^2  per slack
       pair, i.e. f_s must be a separable quadratic in s with constant
       coefficients.  All three conditions are structural, so all three are
       decided here rather than hoped for at solve time.  Once they hold, z
       and Z are plain numbers, evaluated right here and handed to the
       runtime as constants -- which is why nothing after this point, the
       symbolic rewrite and generated C included, ever touches f_s. */
    slacks_ = slack_native_ && Nlpsol::ns_ > 0;
    fs_z_.clear();
    fs_Z_.clear();
    if (slacks_) {
      // Deliberately NOT create_function: used here and thrown away, so no
      // codegen dependency, nothing to serialise, no per-solve call.
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
      // z = grad(f_s,s) at s=0 and Z = diag hess(f_s,s,s), both constant by
      // the assertions above, so one evaluation at s=0 settles them for good
      std::vector<DM> fs_arg = {DM::zeros(fs.size1_in(0), fs.size2_in(0)),
                                DM::zeros(fs.size1_in(1), fs.size2_in(1))};
      std::vector<DM> fs_res = fs(fs_arg);
      fs_z_ = densify(fs_res[1]).nonzeros();
      fs_Z_ = densify(diag(fs_res[2])).nonzeros();
      casadi_assert_dev(fs_z_.size()==static_cast<size_t>(2*Nlpsol::ns_));
      casadi_assert_dev(fs_Z_.size()==static_cast<size_t>(2*Nlpsol::ns_));
      for (casadi_int i=0; i<2*Nlpsol::ns_; ++i) {
        casadi_assert(std::isfinite(fs_z_[i]) && std::isfinite(fs_Z_[i]),
          "Ipmc: the slack penalty f_s has a non-finite gradient or Hessian at s=0 "
          "(entry " + str(i) + " of 2*ns). Use a penalty that is finite there, or pass "
          "expand_slacks=True to fall back to the reference expansion.");
      }
    }

    const std::vector<casadi_int>& nx = nxs_;
    const std::vector<casadi_int>& ng = ngs_;
    const std::vector<casadi_int>& nu = nus_;

    // Keep list of erroring rows
    std::set<casadi_int> errors;

    /* The CALLER's constraint Jacobian pattern, which is what the staircase
       is detected in.  Only 'auto' reads it; 'manual' was handed the
       partition and 'none' has none. */
    Sparsity A_;
    if (structure_detection_==STRUCTURE_AUTO) {
      A_ = jacg_sp_;
    } else {
      A_ = Sparsity(ng_, nx_);
    }
    casadi_int na_ = ng_;

    if (struct_cnt>0) {
      casadi_assert(structure_detection_ == STRUCTURE_MANUAL,
        "You must set structure_detection to 'manual' if you set N, nx, nu, ng.");
    }

    if (structure_detection_==STRUCTURE_MANUAL) {
      casadi_assert(struct_cnt==4,
        "You must set all of N, nx, nu, ng.");
    } else if (structure_detection_==STRUCTURE_NONE) {
      N_ = 0;
      nxs_ = {0};
      nus_ = {nx_};
      ngs_ = {ng_};
    } else if (structure_detection_==STRUCTURE_AUTO) {
      casadi_assert(!equality_.empty(),
        "Structure detection auto requires the 'equality' option to be set");
      /* General strategy: look for the x_{k+1} diagonal part in A */

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

    casadi_assert(nx.size()==N_+1, "nx must have length N+1.");
    casadi_assert(nu.size()==N_+1, "nu must have length N+1.");
    casadi_assert(ng.size()==N_+1, "ng must have length N+1.");

    /* ================ cross-stage slack columns, and the rewrite ==========
       slack_S_ rows are [S_g (ngu_ rows) ; S_x (nxu_ rows)], which is NOT the
       interface's own z ordering ([x ; g]).  Row r < ngu_ is caller
       constraint row g[r]; row r >= ngu_ is caller variable x[r-ngu_].  Both
       are mapped to a stage through the partition just found: stage k owns
       path rows [path_off[k], path_off[k]+ng[k]) and variables
       [col_off[k], col_off[k]+nx[k]+nu[k]); a constraint row outside every
       path block is a gap-closing (dynamics) row.                          */
    n_lift_ = 0;
    nxt_ = nx_;
    nat_ = ng_;
    nxs_u_ = nxs_;
    ngs_u_ = ngs_;
    lift_x_.clear(); lift_xfree_.clear(); lift_row_.clear();
    lift_kind_.clear(); lift_ent_.clear(); lift_mflat_.clear();
    lift_m_.clear(); lift_slot_.clear();
    lift_a_src_.clear(); lift_h_src_.clear();
    jacg_lift_sp_ = Sparsity();
    hess_lift_sp_ = Sparsity();

    std::vector<casadi_int> col_off(N_+1, 0), path_off(N_+1, 0), dyn_off(N_+1, 0);
    {
      casadi_int oc = 0, orow = 0;
      for (casadi_int k=0;k<N_;++k) {
        col_off[k] = oc; oc += nx[k]+nu[k];
        dyn_off[k] = orow; path_off[k] = orow + nx[k+1];
        orow += nx[k+1]+ng[k];
      }
      col_off[N_] = oc;
      path_off[N_] = orow;
    }
    // caller g row / caller x -> the slack column softening it, -1 when hard
    std::vector<casadi_int> ug_col(ngu_, -1), ux_col(nxu_, -1);
    // slack column -> its stage, -1 when it spans several and must be lifted;
    // and -> the first of its two helper components, -1 when it is not lifted
    std::vector<casadi_int> col_stage, col_ent;
    if (slacks_) {
      casadi_int ns = Nlpsol::ns_;
      std::vector<casadi_int> g_stage(ngu_, -1), x_stage(nxu_, -1);
      for (casadi_int k=0;k<N_+1;++k) {
        for (casadi_int i=0;i<ng[k];++i) g_stage[path_off[k]+i] = k;
        for (casadi_int i=0;i<nx[k]+nu[k];++i) x_stage[col_off[k]+i] = k;
      }
      col_stage.assign(ns, -1);
      col_ent.assign(ns, -1);
      const casadi_int* colind = slack_S_.colind();
      const casadi_int* row = slack_S_.row();
      casadi_int si = GlobalOptions::start_index;
      for (casadi_int j=0;j<ns;++j) {
        std::vector<casadi_int> stages;
        for (casadi_int el=colind[j]; el<colind[j+1]; ++el) {
          casadi_int r = row[el];
          casadi_int k;
          std::string what;
          if (r < ngu_) {
            k = g_stage[r];
            what = "constraint row g[" + str(r+si) + "]";
            casadi_assert(ug_col[r]<0,
              "Ipmc: " + what + " is softened by more than one slack column (" +
              str(ug_col[r]+si) + " and " + str(j+si) + "). ipmc maps every inequality "
              "row to at most one slack pair, so the columns of S must have disjoint rows.");
            ug_col[r] = j;
          } else {
            casadi_int i = r - ngu_;
            k = x_stage[i];
            what = "the simple bound on x[" + str(i+si) + "]";
            casadi_assert(ux_col[i]<0,
              "Ipmc: " + what + " is softened by more than one slack column (" +
              str(ux_col[i]+si) + " and " + str(j+si) + "). ipmc maps every inequality "
              "row to at most one slack pair, so the columns of S must have disjoint rows.");
            ux_col[i] = j;
          }
          casadi_assert(k>=0,
            "Ipmc: slack column " + str(j+si) + " softens " + what + ", which is a "
            "gap-closing (dynamics) row of the detected OCP structure. Only path "
            "constraints and simple bounds can be softened. Structure is: N " + str(N_) +
            ", nx " + str(nxs_) + ", nu " + str(nus_) + ", ng " + str(ngs_) + ".");
          if (std::find(stages.begin(), stages.end(), k)==stages.end()) stages.push_back(k);
        }
        /* A column local to ONE stage stays an ipmc slack pair: its rank-1
           Schur correction is stage-local, so the recursion can absorb it.
           A column spanning SEVERAL stages -- an L-infinity (peak) penalty --
           couples every stage it touches and cannot be expressed through
           ipmc's stage-local soft_idx, so it is LIFTED: one helper state
           component per softened side, see the runtime header. */
        if (stages.size()>1) {
          col_ent[j] = n_lift_;
          lift_slot_.push_back(j);      // lower side
          lift_slot_.push_back(ns+j);   // upper side
          n_lift_ += 2;
        } else {
          col_stage[j] = stages.empty() ? 0 : stages[0];
        }
      }
    }

    /* ---- the rewrite itself. Nothing below this point knows about slacks:
       what comes out is an ordinary OCP, one stage wider and a few rows
       longer, and it goes down the un-lifted path unchanged. ------------- */
    if (n_lift_>0) {
      // variables: [x_0, m, u_0, x_1, m, u_1, ..., x_N, m]; the helpers go at
      // the END of the state block, which is what makes them eligible as nxc
      lift_x_.assign(nx_, -1);
      lift_m_.assign((N_+1)*n_lift_, -1);
      casadi_int a = 0;
      for (casadi_int k=0;k<N_+1;++k) {
        for (casadi_int j=0;j<nx[k];++j) lift_x_[col_off[k]+j] = a++;
        for (casadi_int e=0;e<n_lift_;++e) lift_m_[k*n_lift_+e] = a++;
        for (casadi_int j=0;j<nu[k];++j) lift_x_[col_off[k]+nx[k]+j] = a++;
      }
      nxt_ = a;
      lift_xfree_.assign(nx_, 0);
      std::vector<casadi_int> ng_new(N_+1, 0);
      for (casadi_int k=0;k<N_+1;++k) {
        if (k<N_) {
          // the gap-closing rows of transition k, then m_{k+1} - m_k = 0
          // AFTER them, so the dynamics row order matches the state order
          for (casadi_int j=0;j<nx[k+1];++j) {
            lift_kind_.push_back(0); lift_row_.push_back(nx_+dyn_off[k]+j);
            lift_ent_.push_back(-1); lift_mflat_.push_back(-1);
          }
          for (casadi_int e=0;e<n_lift_;++e) {
            lift_kind_.push_back(3); lift_row_.push_back(-1);
            lift_ent_.push_back(-1); lift_mflat_.push_back(k*n_lift_+e);
          }
        }
        // the path rows of stage k; a row softened by a lifted column SPLITS
        // into c + m_l >= lb and c - m_u <= ub, because the two sides need
        // different helper states
        for (casadi_int i=0;i<ng[k];++i) {
          casadi_int r = path_off[k]+i, j = ug_col[r];
          if (j>=0 && col_ent[j]>=0) {
            for (casadi_int h=0;h<2;++h) {
              lift_kind_.push_back(1+h); lift_row_.push_back(nx_+r);
              lift_ent_.push_back(col_ent[j]+h);
              lift_mflat_.push_back(k*n_lift_+col_ent[j]+h);
            }
            ng_new[k] += 2;
          } else {
            lift_kind_.push_back(0); lift_row_.push_back(nx_+r);
            lift_ent_.push_back(-1); lift_mflat_.push_back(-1);
            ng_new[k] += 1;
          }
        }
        // a softened SIMPLE bound couples its variable to m once relaxed, so
        // it cannot stay a simple bound: it becomes a path row of the stage
        // that owns the variable, and lbx/ubx are opened to +-inf
        for (casadi_int i=0;i<nx[k]+nu[k];++i) {
          casadi_int xi = col_off[k]+i, j = ux_col[xi];
          if (j>=0 && col_ent[j]>=0) {
            lift_xfree_[xi] = 1;
            for (casadi_int h=0;h<2;++h) {
              lift_kind_.push_back(1+h); lift_row_.push_back(xi);
              lift_ent_.push_back(col_ent[j]+h);
              lift_mflat_.push_back(k*n_lift_+col_ent[j]+h);
            }
            ng_new[k] += 2;
          }
        }
      }
      nat_ = lift_kind_.size();
      for (casadi_int k=0;k<N_+1;++k) nxs_[k] += n_lift_;
      ngs_ = ng_new;
      if (verbose_) {
        casadi_message("Lifted " + str(n_lift_/2) + " cross-stage slack column(s) into "
          + str(n_lift_) + " helper state components; rewritten problem has nx "
          + str(nxt_) + ", ng " + str(nat_) + ".");
      }
    }

    convexify_ = false;
    if (exact_hessian_ && convexify_strategy!="none") {
      convexify_ = true;
      Dict convexify_opts;
      convexify_opts["strategy"] = convexify_strategy;
      convexify_opts["margin"] = convexify_margin;
      convexify_opts["max_iter_eig"] = max_iter_eig;
      convexify_opts["verbose"] = verbose_;
      hesslag_sp_ = Convexify::setup(convexify_data_, hesslag_sp_, convexify_opts);
    }
    /* the lifted patterns and nonzero maps, from the (convexified) caller
       patterns; what is validated and handed over below is the lifted
       problem's jacobian pattern */
    if (n_lift_>0) lift_patterns();
    A_ = n_lift_>0 ? jacg_lift_sp_ : jacg_sp_;
    na_ = nat_;
    if (debug_) {
      A_.to_file("debug_ipmc_actual.mtx");
    }

    /* The constant-state declaration.  It is checked here rather than left to
       ipmc_create so that the message can name the stage that has no room for
       it; ipmc only answers with an error code. */
    if (nxc_given) {
      casadi_assert(structure_detection_!=STRUCTURE_NONE,
        "nxc needs structure_detection 'manual' or 'auto': with 'none' there is "
        "no horizon and no state to be constant.");
      for (casadi_int k=0;k<N_+1;++k) {
        casadi_assert(nxc_>=0 && nxc_<=nxs_u_[k],
          "nxc = " + str(nxc_) + " is not in [0, nx[" + str(k) +
          "]] = [0, " + str(nxs_u_[k]) + "].");
      }
    }

    if (verbose_) {
      casadi_message("Using structure: N " + str(N_) + ", nx " + str(nx) + ", "
            "nu " + str(nu) + ", ng " + str(ng) + ".");
    }

    // For debugging purposes
    std::vector< casadi_ocp_block > A_blocks, B_blocks, C_blocks, D_blocks;

    /* Disassemble A input into:
       A B I
       C D
           A B I
           C D
               C D
    */
    casadi_int offset_r = 0, offset_c = 0;
    for (casadi_int k=0;k<N_;++k) { // Loop over blocks
      AB_blocks_.push_back({offset_r,        offset_c,            nx[k+1], nx[k]+nu[k]});
      CD_blocks_.push_back({offset_r+nx[k+1], offset_c,           ng[k], nx[k]+nu[k]});
      A_blocks.push_back({offset_r, offset_c, nx[k+1], nx[k]});
      B_blocks.push_back({offset_r, offset_c+nx[k], nx[k+1], nu[k]});
      C_blocks.push_back({offset_r+nx[k+1], offset_c, ng[k], nx[k]});
      D_blocks.push_back({offset_r+nx[k+1], offset_c+nx[k], ng[k], nu[k]});
      offset_c+= nx[k]+nu[k];
      I_blocks_.push_back({offset_r, offset_c, nx[k+1], nx[k+1]});
      offset_r+= nx[k+1]+ng[k];
    }
    CD_blocks_.push_back({offset_r, offset_c,           ng[N_], nx[N_]+nu[N_]});
    C_blocks.push_back({offset_r, offset_c,           ng[N_], nx[N_]});
    D_blocks.push_back({offset_r, offset_c+nx[N_],           ng[N_], nu[N_]});

    casadi_int offset = 0;
    AB_offsets_.push_back(0);
    for (auto e : AB_blocks_) {
      offset += e.rows*e.cols;
      AB_offsets_.push_back(offset);
    }
    offset = 0;
    CD_offsets_.push_back(0);
    for (auto e : CD_blocks_) {
      offset += e.rows*e.cols;
      CD_offsets_.push_back(offset);
    }

    ABsp_ = blocksparsity(nat_, nxt_, AB_blocks_);
    CDsp_ = blocksparsity(nat_, nxt_, CD_blocks_);
    Isp_ = blocksparsity(nat_, nxt_, I_blocks_, true);

    Sparsity total = ABsp_ + CDsp_ + Isp_;

    if (debug_) {
      total.to_file("debug_ipmc_expected.mtx");
      blocksparsity(nat_, nxt_, A_blocks).to_file("debug_ipmc_A.mtx");
      blocksparsity(nat_, nxt_, B_blocks).to_file("debug_ipmc_B.mtx");
      blocksparsity(nat_, nxt_, C_blocks).to_file("debug_ipmc_C.mtx");
      blocksparsity(nat_, nxt_, D_blocks).to_file("debug_ipmc_D.mtx");
      Isp_.to_file("debug_ipmc_I.mtx");
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
    casadi_assert_dev(total.nnz() == ABsp_.nnz() + CDsp_.nnz() + Isp_.nnz());

    /* Disassemble H input into:
       Q S'
       S R
           Q S'
           S R
    */
    offset = 0;
    for (casadi_int k=0;k<N_+1;++k) { // Loop over blocks
      RSQ_blocks_.push_back({offset, offset,       nx[k]+nu[k], nx[k]+nu[k]});
      offset+= nx[k]+nu[k];
    }
    RSQsp_ = blocksparsity(nxt_, nxt_, RSQ_blocks_);

    offset = 0;
    RSQ_offsets_.push_back(0);
    for (auto e : RSQ_blocks_) {
      offset += e.rows*e.cols;
      RSQ_offsets_.push_back(offset);
    }

    /* ---- the slack maps, in the REWRITTEN row/variable space -------------
       The classification itself was done above, before the rewrite, because
       it is what decided the rewrite; all that is left here is to state the
       surviving stage-local columns in the indices ipmc will actually use.
       A lifted column is not a ipmc slack pair at all: it carries
       slack_stage_ = -1 and appears in no map, while its rows -- now split
       and hard -- carry -1 in slack_g_/slack_x_ like any other hard row. */
    if (slacks_) {
      casadi_int ns = Nlpsol::ns_;
      slack_stage_ = col_stage;
      slack_g_.assign(nat_, -1);
      slack_x_.assign(nxt_, -1);
      if (n_lift_>0) {
        for (casadi_int r=0;r<nat_;++r) {
          if (lift_kind_[r]==0) slack_g_[r] = ug_col[lift_row_[r]-nx_];
        }
        for (casadi_int i=0;i<nx_;++i) {
          if (!lift_xfree_[i]) slack_x_[lift_x_[i]] = ux_col[i];
        }
      } else {
        slack_g_ = ug_col;
        slack_x_ = ux_col;
      }
      // Group the surviving columns by stage: slack_perm_ is ipmc's own order
      slack_ns_.assign(N_+1, 0);
      for (casadi_int j=0;j<ns;++j) {
        if (slack_stage_[j]>=0) slack_ns_[slack_stage_[j]]++;
      }
      slack_offs_.assign(N_+2, 0);
      for (casadi_int k=0;k<N_+1;++k) slack_offs_[k+1] = slack_offs_[k] + slack_ns_[k];
      std::vector<casadi_int> cnt(slack_offs_.begin(), slack_offs_.end()-1);
      slack_perm_.assign(slack_offs_[N_+1], -1);
      for (casadi_int j=0;j<ns;++j) {
        if (slack_stage_[j]>=0) slack_perm_[cnt[slack_stage_[j]]++] = j;
      }
      if (verbose_) {
        casadi_message("Native slacks: ns " + str(ns) + ", per stage " + str(slack_ns_) +
          ", lifted helper states " + str(n_lift_) + ".");
      }
    }

    set_ipmc_prob();

    // Allocate memory
    casadi_int sz_arg, sz_res, sz_w, sz_iw;
    casadi_ipmc_work(&p_, &sz_arg, &sz_res, &sz_iw, &sz_w);
    if (n_lift_>0) {
      /* the REWRITTEN problem's own z/lbz/ubz/lam; casadi_ipmc_work sizes
         everything else off p_.nlp, which already IS the rewritten one */
      casadi_int a2, r2, i2, w2;
      casadi_nlpsol_work(&p_nlp_lift_, &a2, &r2, &i2, &w2);
      sz_arg += a2; sz_res += r2; sz_iw += i2; sz_w += w2;
    }

    alloc_arg(sz_arg, true);
    alloc_res(sz_res, true);
    alloc_iw(sz_iw, true);
    alloc_w(sz_w, true);
  }

  int IpmcInterface::init_mem(void* mem) const {
    if (Nlpsol::init_mem(mem)) return 1;
    if (!mem) return 1;
    auto m = static_cast<IpmcMemory*>(mem);
    m->self = this;
    casadi_ipmc_init_mem(&m->d);

    return 0;
  }

  void IpmcInterface::free_mem(void* mem) const {
    auto m = static_cast<IpmcMemory*>(mem);
    casadi_ipmc_free_mem(&m->d);
    delete static_cast<IpmcMemory*>(mem);
  }

  /** \brief Set the (persistent) work vectors */
  void IpmcInterface::set_work(void* mem, const double**& arg, double**& res,
                                casadi_int*& iw, double*& w) const {
    auto m = static_cast<IpmcMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    m->d.prob = &p_;
    m->d.nlp = &m->d_nlp;
    m->d.rewrite.user = &m->d_nlp;

    /* When a cross-stage slack column was lifted, the problem handed to ipmc
       is not the caller's: it gets its own z/lbz/ubz/lam, filled by
       casadi_ipmc_lift_expand at the top of every solve.  Mirrored by
       set_ipmc_prob(CodeGenerator&); the w consumption must match. */
    if (n_lift_>0) {
      m->d_nlp_lift = m->d_nlp;
      m->d_nlp_lift.prob = &p_nlp_lift_;
      casadi_nlpsol_set_work(&m->d_nlp_lift, &arg, &res, &iw, &w);
      m->d.nlp = &m->d_nlp_lift;
    }

    casadi_ipmc_set_work(&m->d, &arg, &res, &iw, &w);

    // Native slacks live on NlpsolMemory, not in casadi_nlpsol_data::z.
    // m->s / m->lam_s point straight at res[...] and may be null; the
    // slack_* scratch buffers never are.
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

  int IpmcInterface::solve(void* mem) const {
    auto m = static_cast<IpmcMemory*>(mem);

    /* Cache the solver: presolve (re)creates it only when needed.  Must be
       read AFTER presolve: a re-created solver starts with ipmc's defaults,
       so deciding "is this a new solver" from d.solver==0 on entry would
       silently drop every plugin option from the second solve onwards. */
    bool new_solver;
    /* The caller's problem onto the one that was handed over.  A no-op
       unless a cross-stage slack column was lifted, and then a pure vector
       scatter: bounds and ubs are runtime inputs, so the rewritten bounds
       cannot be baked in at construction the way the expressions were. */
    casadi_ipmc_rewrite_expand(&p_.rewrite, &m->d.rewrite, m->d.nlp,
                                 m->d.slack_ubs, m->d.slack_s, m->d.slack_lam_s);
    /* ipmc's diagnostic output goes through a file-static write callback, so
       two threads inside ipmc at once can send each other's messages to the
       wrong stream.  Serialize the create, the call that can print before
       the solver even exists. */
    {
#ifdef CASADI_WITH_THREAD
      static std::mutex mutex_ipmc_create;
      std::lock_guard<std::mutex> lock(mutex_ipmc_create);
#endif //CASADI_WITH_THREAD
      casadi_ipmc_presolve(&m->d);
    }
    new_solver = m->d.solver_created;

    if (m->d.solver == 0) {
      casadi_error("Ipmc: the solver refused the problem description: " +
        ipmc_soft_error_message(m->d.soft_error) +
        " (IpmcError == " + str(m->d.soft_error) + ").");
    }

    // Set options only when a new solver was created (options persist across solves)
    if (new_solver) {
      for (const auto& kv : opts_) {
        switch (ipmc_option_type(kv.first.c_str())) {
          case 0:
            ipmc_set_option_double(m->d.solver, kv.first.c_str(), kv.second);
            break;
          case 1:
            ipmc_set_option_int(m->d.solver, kv.first.c_str(), kv.second.to_int());
            break;
          case 2:
            ipmc_set_option_bool(m->d.solver, kv.first.c_str(), kv.second.to_bool());
            break;
          case 3:
            {
              std::string s = kv.second.to_string();
              ipmc_set_option_string(m->d.solver, kv.first.c_str(), s.c_str());
            }
            break;
          case -1:
            casadi_error("Ipmc option not supported: " + kv.first);
          default:
            casadi_error("Unknown option type.");
        }
      }
    }

    /* The solve loop: ipmc asks, we answer.  IpmcInterface::codegen_body
       emits the SAME loop as C; the two are deliberately separate copies, so
       each can make its oracle calls the way its own world does --
       calc_function() here, a direct call to the generated function there --
       and no function pointer has to live in casadi_ipmc_prob.  Keep the two
       copies in step; casadi_ipmc_finish, which needs no oracle, is
       shared. */
    {
      casadi_ipmc_data<double>* d = &m->d;
      const IpmcLayout* str = d->layout;
      IpmcEval e_;
      casadi_int i, n_ux;
      int stop = 0;

      d->unified_return_status = SOLVER_RET_UNKNOWN;
      d->success = 0;

      /* stop is the iteration callback's answer; it sits in the loop
         CONDITION because a break inside the switch would just step ipmc
         again. */
      ipmc_start(d->solver);
      while (!stop && ipmc_step(d->solver, &e_)!=IPMC_DONE) {
        switch (e_.request) {
        /* EVERY oracle call the solve makes is below, written out where the
           request arrives: read the point into d->x, gather the CALLER's
           point/multipliers into the u-buffers, call, lift the output into
           the rewritten problem's quantity, then hand it to a pack_ routine
           that only rearranges it into ipmc's layout.  The u-buffers alias
           the lifted ones and every lift_ is a no-op when nothing was
           rewritten. */
        case IPMC_EVAL_OBJ:
          casadi_ipmc_read_primal_data(&p_, e_.x, d->x, str);
          casadi_ipmc_lift_x(&p_, d);
          m->arg[0] = d->xu;
          m->arg[1] = d->nlp->p;
          m->res[0] = e_.obj;
          calc_function(m, "nlp_f");
          casadi_ipmc_lift_obj(&p_, d, e_.obj);
          *e_.obj *= e_.obj_scale;
          break;
        case IPMC_EVAL_OBJ_GRAD:
          casadi_ipmc_read_primal_data(&p_, e_.x, d->x, str);
          casadi_ipmc_lift_x(&p_, d);
          m->arg[0] = d->xu;
          m->arg[1] = d->nlp->p;
          m->res[0] = d->gu;
          calc_function(m, "nlp_grad_f");
          casadi_ipmc_lift_grad_f(&p_, d);
          /* nlp_grad_f writes in the caller's column order; ipmc reads stage
             by stage, and the scaling is ipmc's */
          casadi_ipmc_write_primal_data(&p_, d->g, e_.grad, str);
          casadi_scal(p_.nlp->nx, e_.obj_scale, e_.grad);
          break;
        case IPMC_EVAL_CONSTR_VIOL:
          casadi_ipmc_read_primal_data(&p_, e_.x, d->x, str);
          casadi_ipmc_lift_x(&p_, d);
          m->arg[0] = d->xu;
          m->arg[1] = d->nlp->p;
          m->res[0] = d->gu;
          calc_function(m, "nlp_g");
          casadi_ipmc_lift_g(&p_, d);
          if (!fcallback_.is_null()) casadi_ipmc_snapshot_g(d, str, e_.x);
          casadi_ipmc_pack_constr_viol(d, str, e_.cv);
          break;
        case IPMC_EVAL_CONSTR_JAC:
          casadi_ipmc_read_primal_data(&p_, e_.x, d->x, str);
          casadi_ipmc_lift_x(&p_, d);
          m->arg[0] = d->xu;
          m->arg[1] = d->nlp->p;
          m->res[0] = d->gu;
          m->res[1] = d->au;
          calc_function(m, "nlp_jac_g");
          casadi_ipmc_lift_g(&p_, d);
          casadi_ipmc_lift_jac_g(&p_, d);
          casadi_ipmc_pack_constr_jac(d, str, e_.BAbt, e_.Ggt, e_.Ggt_ineq);
          break;
        case IPMC_EVAL_LAG_HESS:
          /* lam is an INPUT to nlp_hess_l, so the multipliers are scattered
             into d->lam -- and folded onto the caller's rows -- before it
             runs */
          casadi_ipmc_read_primal_data(&p_, e_.x, d->x, str);
          casadi_ipmc_lift_x(&p_, d);
          casadi_ipmc_read_lam(d, str, e_.lam);
          casadi_ipmc_lift_lam(&p_, d);
          m->arg[0] = d->xu;
          m->arg[1] = d->nlp->p;
          m->arg[2] = &e_.obj_scale;
          m->arg[3] = d->lamu;
          m->res[0] = d->gu;
          m->res[1] = d->hu;
          calc_function(m, "nlp_hess_l");
          casadi_ipmc_lift_hess_l(&p_, d, e_.obj_scale);
          casadi_ipmc_pack_lag_hess(d, str, e_.lam, e_.RSQrqt);
          break;
        /* ipmc hands an accepted iterate over whole: x, lam and the
           objective are guaranteed to be that one point (see ipmc.h).  Two
           amendments are ours, not ipmc's: the objective carries obj_scale
           and, when a column was lifted, the penalty on the helper columns;
           and g is the snapshot, because ipmc holds the constraint violation
           rather than g itself.  Both are undone in
           casadi_ipmc_report_iterate, along with the move out of the
           rewritten z-space.  The emitted twin has no case for this at all:
           generated C has no casadi Function to call back into. */
        case IPMC_POST_ITERATION:
          if (fcallback_.is_null()) break;
          /* ipmc can step back to an earlier iterate (a rejected watchdog
             step) and redo the iteration from there.  Our snapshot cannot
             come back with it, so it is checked against the iterate being
             reported and retaken when it belongs to some other point --
             one nlp_g on that path only, never once per iteration. */
          n_ux = str->ux_offs[str->K-1] + str->nu[str->K-1] + str->nx[str->K-1];
          for (i=0;i<n_ux;++i) if (d->cb_x[i]!=e_.x[i]) break;
          if (i<n_ux) {
            casadi_ipmc_read_primal_data(&p_, e_.x, d->x, str);
            casadi_ipmc_lift_x(&p_, d);
            m->arg[0] = d->xu;
            m->arg[1] = d->nlp->p;
            m->res[0] = d->gu;
            calc_function(m, "nlp_g");
            casadi_ipmc_lift_g(&p_, d);
            casadi_ipmc_snapshot_g(d, str, e_.x);
          }
          /* No multipliers on a restoration iterate (they would be the
             restoration problem's); the duals stay at the last accepted
             ordinary iterate's, as in ipopt's callback. */
          casadi_ipmc_report_iterate(d, str, e_.x, e_.lam, e_.objective/e_.obj_scale);
          /* Nlpsol::callback implements the whole policy: it rethrows
             KeyboardInterruptException, honours
             iteration_callback_ignore_errors, and returns non-zero to stop.
             A throw unwinding out of this loop is safe: ipmc keeps no state
             a later ipmc_start does not reset. */
          if (callback(m)) {
            ipmc_abort(d->solver);
            stop = 1;
          }
          break;
        /* unreachable (the loop condition tested for it), but naming the
           case keeps -Wswitch quiet */
        case IPMC_DONE:
          break;
        }
      }
      casadi_ipmc_finish(d);
    }

    casadi_ipmc_rewrite_collect(&p_.rewrite, &m->d.rewrite, m->d.nlp,
                                  m->d.slack_s, m->d.slack_lam_s);

    m->success = m->d.success;
    m->unified_return_status = static_cast<UnifiedReturnStatus>(m->d.unified_return_status);

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
    stats["nx"] = nxs_u_;
    stats["nu"] = nus_;
    stats["ng"] = ngs_u_;
    stats["nxc"] = nxc_;
    /* How many helper state components the cross-stage (L-infinity) columns
       were lifted into, so a test can assert the lifting actually ran.  The
       structure above is the CALLER's, not the rewritten one. */
    stats["n_lift"] = n_lift_;
    stats["N"] = N_;
    stats["return_status"] = m->d.return_status;
    return stats;
  }

  void IpmcInterface::codegen_init_mem(CodeGenerator& g) const {
    g << "casadi_ipmc_init_mem(&" + codegen_mem(g) + ");\n";
    g << "return 0;\n";
  }

  void IpmcInterface::codegen_free_mem(CodeGenerator& g) const {
    g << "casadi_ipmc_free_mem(&" + codegen_mem(g) + ");\n";
  }

void IpmcInterface::codegen_declarations(CodeGenerator& g) const {
  Nlpsol::codegen_declarations(g);
  g.add_auxiliary(CodeGenerator::AUX_NLP);
  g.add_auxiliary(CodeGenerator::AUX_INF);
  g.add_auxiliary(CodeGenerator::AUX_MAX);
  g.add_auxiliary(CodeGenerator::AUX_COPY);
  g.add_auxiliary(CodeGenerator::AUX_PROJECT);
  g.add_auxiliary(CodeGenerator::AUX_SCAL);
  g.add_auxiliary(CodeGenerator::AUX_SPARSITY);
  g.add_auxiliary(CodeGenerator::AUX_OCP_BLOCK);
  g.add_auxiliary(CodeGenerator::AUX_DENSIFY);
  g.add_auxiliary(CodeGenerator::AUX_SPARSIFY);
  g.add_auxiliary(CodeGenerator::AUX_SCALED_COPY);
  g.add_auxiliary(CodeGenerator::AUX_AXPY);
  g.add_auxiliary(CodeGenerator::AUX_PRINTF);
  g.add_dependency(get_function("nlp_f"));
  g.add_dependency(get_function("nlp_grad_f"));
  g.add_dependency(get_function("nlp_g"));
  g.add_dependency(get_function("nlp_jac_g"));
  if (exact_hessian_) {
    g.add_dependency(get_function("nlp_hess_l"));
  }
  g.add_include("ipmc/ipmc.h");
  g.add_include("stdlib.h", false);

  std::string name = "ipmc_cb_write";
  std::string f = g.shorthand(name);

  g << "void " << f
    << "(const char* msg, int num) {\n";
  g.flush(g.body);
  g.scope_enter();
  g << "CASADI_PRINTF(\"%.*s\", num, msg);\n";
  g.scope_exit();
  g << "}\n";

  name = "ipmc_cb_flush";
  f = g.shorthand(name);

  g << "void " << f
    << "(void) {\n";
  g.flush(g.body);
  g.scope_enter();
  g.scope_exit();
  g << "}\n";
}

void IpmcInterface::codegen_body(CodeGenerator& g) const {
  codegen_body_enter(g);
  g.auxiliaries << g.sanitize_source(ipmc_runtime_str, {"casadi_real"});

  g.local("d", "struct casadi_ipmc_data*");
  g.init_local("d", "&" + codegen_mem(g));
  g.local("p", "struct casadi_ipmc_prob");
  set_ipmc_prob(g);

  g << "casadi_ipmc_set_work(d, &arg, &res, &iw, &w);\n";

  // Mirror of IpmcInterface::set_work's native-slack branch; the locals are
  // emitted by Nlpsol::codegen_slack_native_enter from codegen_body_enter.
  if (slacks_) {
    g << "d->slack_s = nlp_slack_s;\n";
    g << "d->slack_lam_s = nlp_slack_lam_s;\n";
    g << "d->slack_ubs = nlp_slack_ubs;\n";
  } else {
    g << "d->slack_s = 0;\n";
    g << "d->slack_lam_s = 0;\n";
    g << "d->slack_ubs = 0;\n";
  }

  g << "casadi_ipmc_rewrite_expand(&p.rewrite, &d->rewrite, d->nlp, "
       "d->slack_ubs, d->slack_s, d->slack_lam_s);\n";

  // Cache the solver: presolve (re)creates it only when needed
  g << "{\n";
  g << "int new_solver;\n";
  // see IpmcInterface::solve: ipmc's diagnostic output goes through a
  // file-static write callback, so the create is serialized across threads.
  if (g.thread_safe()) {
    Function F = shared_from_this<Function>();
    std::string mutex_name = codegen_name(g, false) + "_ipmc_create_mutex";
    g.define_local_mutex(F, mutex_name);
    std::string mtx = g.local_mutex(F, mutex_name);
    g << "CASADI_MUTEX_LOCK(&" << mtx << ");\n";
    g << "casadi_ipmc_presolve(d);\n";
    g << "CASADI_MUTEX_UNLOCK(&" << mtx << ");\n";
  } else {
    g << "casadi_ipmc_presolve(d);\n";
  }
  // see IpmcInterface::solve: only presolve knows whether it re-created
  g << "new_solver = d->solver_created;\n";
  g << "if (new_solver) {\n";

  for (const auto& kv : opts_) {
    switch (ipmc_option_type(kv.first.c_str())) {
      case 0:
        g << "ipmc_set_option_double(d->solver, \"" + kv.first + "\", "
              + g.constant(kv.second.to_double()) + ");\n";
        break;
      case 1:
        g << "ipmc_set_option_int(d->solver, \"" + kv.first + "\", "
              + str(kv.second.to_int()) + ");\n";
        break;
      case 2:
        g << "ipmc_set_option_bool(d->solver, \"" + kv.first + "\", "
              + str(static_cast<int>(kv.second.to_bool())) + ");\n";
        break;
      case 3:
        {
          std::string s = kv.second.to_string();
          g << "ipmc_set_option_bool(d->solver, \"" + kv.first + "\", \""
              + s + "\");\n";
        }
        break;
      case -1:
        casadi_error("Ipmc option not supported: " + kv.first);
      default:
        casadi_error("Unknown option type.");
    }
  }

  g << "}\n";
  g << "}\n";
  // ipmc_create can refuse the description (an inconsistent soft structure) or
  // simply run out of memory, so this is not conditional on slacks_.
  g << "if (d->solver == 0) {\n";
  g << "  CASADI_PRINTF(\"Ipmc: the solver refused the problem description "
       "(IpmcError == %d).\\n\", d->soft_error);\n";
  g << "  return 1;\n";
  g << "}\n";
  /* The solve loop: the emitted twin of the C++ loop in
     IpmcInterface::solve -- deliberately a separate copy, see the comment
     there.  Keep the two in step. */
  g.local("str", "const struct IpmcLayout", "*");
  g.local("e_", "struct IpmcEval");
  g << "str = d->layout;\n";
  g << "d->unified_return_status = " << static_cast<int>(SOLVER_RET_UNKNOWN) << ";\n";
  g << "d->success = 0;\n";
  g << "ipmc_start(d->solver);\n";
  g << "while (ipmc_step(d->solver, &e_)!=IPMC_DONE) {\n";
  g << "switch (e_.request) {\n";

  g << "case IPMC_EVAL_OBJ:\n";
  g << "casadi_ipmc_read_primal_data(&p, e_.x, d->x, str);\n";
  g << "casadi_ipmc_lift_x(&p, d);\n";
  g << "d->arg[0] = d->xu;\n";
  g << "d->arg[1] = d->nlp->p;\n";
  g << "d->res[0] = e_.obj;\n";
  g << g(get_function("nlp_f"), "d->arg", "d->res", "d->iw", "d->w") << ";\n";
  g << "casadi_ipmc_lift_obj(&p, d, e_.obj);\n";
  g << "*e_.obj *= e_.obj_scale;\n";
  g << "break;\n";

  g << "case IPMC_EVAL_OBJ_GRAD:\n";
  g << "casadi_ipmc_read_primal_data(&p, e_.x, d->x, str);\n";
  g << "casadi_ipmc_lift_x(&p, d);\n";
  g << "d->arg[0] = d->xu;\n";
  g << "d->arg[1] = d->nlp->p;\n";
  g << "d->res[0] = d->gu;\n";
  g << g(get_function("nlp_grad_f"), "d->arg", "d->res", "d->iw", "d->w") << ";\n";
  g << "casadi_ipmc_lift_grad_f(&p, d);\n";
  g << "casadi_ipmc_write_primal_data(&p, d->g, e_.grad, str);\n";
  g << "casadi_scal(" << nxt_ << ", e_.obj_scale, e_.grad);\n";
  g << "break;\n";

  g << "case IPMC_EVAL_CONSTR_VIOL:\n";
  g << "casadi_ipmc_read_primal_data(&p, e_.x, d->x, str);\n";
  g << "casadi_ipmc_lift_x(&p, d);\n";
  g << "d->arg[0] = d->xu;\n";
  g << "d->arg[1] = d->nlp->p;\n";
  g << "d->res[0] = d->gu;\n";
  g << g(get_function("nlp_g"), "d->arg", "d->res", "d->iw", "d->w") << ";\n";
  g << "casadi_ipmc_lift_g(&p, d);\n";
  g << "casadi_ipmc_pack_constr_viol(d, str, e_.cv);\n";
  g << "break;\n";

  g << "case IPMC_EVAL_CONSTR_JAC:\n";
  g << "casadi_ipmc_read_primal_data(&p, e_.x, d->x, str);\n";
  g << "casadi_ipmc_lift_x(&p, d);\n";
  g << "d->arg[0] = d->xu;\n";
  g << "d->arg[1] = d->nlp->p;\n";
  g << "d->res[0] = d->gu;\n";
  g << "d->res[1] = d->au;\n";
  g << g(get_function("nlp_jac_g"), "d->arg", "d->res", "d->iw", "d->w") << ";\n";
  g << "casadi_ipmc_lift_g(&p, d);\n";
  g << "casadi_ipmc_lift_jac_g(&p, d);\n";
  g << "casadi_ipmc_pack_constr_jac(d, str, e_.BAbt, e_.Ggt, e_.Ggt_ineq);\n";
  g << "break;\n";

  g << "case IPMC_EVAL_LAG_HESS:\n";
  // lam is an INPUT to nlp_hess_l, so the scatter and fold run BEFORE the call
  g << "casadi_ipmc_read_primal_data(&p, e_.x, d->x, str);\n";
  g << "casadi_ipmc_lift_x(&p, d);\n";
  g << "casadi_ipmc_read_lam(d, str, e_.lam);\n";
  g << "casadi_ipmc_lift_lam(&p, d);\n";
  g << "d->arg[0] = d->xu;\n";
  g << "d->arg[1] = d->nlp->p;\n";
  g << "d->arg[2] = &e_.obj_scale;\n";
  g << "d->arg[3] = d->lamu;\n";
  g << "d->res[0] = d->gu;\n";
  g << "d->res[1] = d->hu;\n";
  g << g(get_function("nlp_hess_l"), "d->arg", "d->res", "d->iw", "d->w") << ";\n";
  g << "casadi_ipmc_lift_hess_l(&p, d, e_.obj_scale);\n";
  g << "casadi_ipmc_pack_lag_hess(d, str, e_.lam, e_.RSQrqt);\n";
  g << "break;\n";

  /* No IPMC_POST_ITERATION branch, and no `stop': generated C has no casadi
     Function to report an iterate to, so nothing can ask it to stop early.
     Naming the two cases keeps -Wswitch quiet. */
  g << "case IPMC_POST_ITERATION:\n";
  g << "case IPMC_DONE:\n";
  g << "break;\n";

  g << "}\n";
  g << "}\n";
  g << "casadi_ipmc_finish(d);\n";
  g << "casadi_ipmc_rewrite_collect(&p.rewrite, &d->rewrite, d->nlp, "
       "d->slack_s, d->slack_lam_s);\n";

  codegen_body_exit(g);

  if (error_on_fail_) {
    g << "return d->unified_return_status;\n";
  } else {
    g << "return 0;\n";
  }
}

std::vector<casadi_int> ipmc_blocks_pack(const std::vector<casadi_ocp_block>& blocks) {
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

void IpmcInterface::set_ipmc_prob() {
  /* Everything below describes the problem HANDED OVER, which is the
     rewritten one when a cross-stage slack column was lifted; the lift_*
     maps at the bottom are what carries the caller's z-space onto it. */
  p_nlp_lift_ = p_nlp_;
  p_nlp_lift_.nx = nxt_;
  p_nlp_lift_.ng = nat_;
  p_nlp_lift_.detect_bounds.ng = 0;
  p_.nlp = n_lift_>0 ? &p_nlp_lift_ : &p_nlp_;
  p_.nx  = get_ptr(nxs_);
  p_.nu  = get_ptr(nus_);
  /* The caller's declared constants occupy [nx-nxc, nx) of every stage and
     the helper components [nx, nx+n_lift): adjacent, hence ONE contiguous
     trailing block of nxc + n_lift -- exactly ipmc's one-number contract.
     The helper rows are real rows of the rewritten Jacobian, so
     casadi_ipmc_pack_BAbt's identity check covers them with no special
     case. */
  p_.nxc = nxc_ + n_lift_;
  p_.ABsp = ABsp_;
  p_.AB_offsets = get_ptr(AB_offsets_);
  p_.CDsp = CDsp_;
  p_.CD_offsets = get_ptr(CD_offsets_);
  p_.RSQsp = RSQsp_;
  p_.RSQ_offsets = get_ptr(RSQ_offsets_);
  p_.Isp = Isp_;
  p_.I_offsets = get_ptr(I_offsets_);

  p_.AB = get_ptr(AB_blocks_);
  p_.CD = get_ptr(CD_blocks_);
  p_.RSQ = get_ptr(RSQ_blocks_);
  p_.I = get_ptr(I_blocks_);
  p_.N = N_;

  /* the LIFTED patterns when a column was lifted: the oracle's outputs are
     rearranged onto them by the runtime's lift_jac_g / lift_hess_l */
  p_.sp_a = n_lift_>0 ? jacg_lift_sp_ : jacg_sp_;
  p_.sp_h = n_lift_>0 && exact_hessian_ ? hess_lift_sp_ : hesslag_sp_;

  p_.ns = slacks_ ? Nlpsol::ns_ : 0;
  if (slacks_) {
    p_.slack_g = get_ptr(slack_g_);
    p_.slack_x = get_ptr(slack_x_);
    p_.slack_perm = get_ptr(slack_perm_);
    p_.slack_offs = get_ptr(slack_offs_);
    p_.fs_z = get_ptr(fs_z_);
    p_.fs_Z = get_ptr(fs_Z_);
  } else {
    p_.slack_g = nullptr;
    p_.slack_x = nullptr;
    p_.slack_perm = nullptr;
    p_.slack_offs = nullptr;
    p_.fs_z = nullptr;
    p_.fs_Z = nullptr;
  }

  p_.rewrite.ne = n_lift_;
  p_.rewrite.nx = nx_;
  p_.rewrite.ng = ng_;
  p_.rewrite.K = N_+1;
  p_.rewrite.x = lift_x_.empty() ? nullptr : get_ptr(lift_x_);
  p_.rewrite.xfree = lift_xfree_.empty() ? nullptr : get_ptr(lift_xfree_);
  p_.rewrite.row = lift_row_.empty() ? nullptr : get_ptr(lift_row_);
  p_.rewrite.kind = lift_kind_.empty() ? nullptr : get_ptr(lift_kind_);
  p_.rewrite.ent = lift_ent_.empty() ? nullptr : get_ptr(lift_ent_);
  p_.rewrite.mflat = lift_mflat_.empty() ? nullptr : get_ptr(lift_mflat_);
  p_.rewrite.m = lift_m_.empty() ? nullptr : get_ptr(lift_m_);
  p_.rewrite.slot = lift_slot_.empty() ? nullptr : get_ptr(lift_slot_);
  p_.rewrite.a_src = lift_a_src_.empty() ? nullptr : get_ptr(lift_a_src_);
  p_.rewrite.h_src = lift_h_src_.empty() ? nullptr : get_ptr(lift_h_src_);
  p_.rewrite.nnz_au = jacg_sp_.nnz();
  p_.rewrite.nnz_hu = exact_hessian_ ? hesslag_sp_.nnz() : 0;

  p_.write = &casadi_c_logger_write;
  p_.flush = &casadi_c_logger_flush;

  casadi_ipmc_setup(&p_);
}

  void codegen_unpack_block(CodeGenerator& g, const std::string& name,
      const std::vector<casadi_ocp_block>& blocks) {
    casadi_int sz = blocks.size();
    if (sz==0) sz++;
    std::string n = "block_" + name + "[" + str(sz) + "]";
    g.local(n, "static struct casadi_ocp_block");
    g << "p." << name << " = block_" + name + ";\n";
    g << "casadi_unpack_ocp_blocks(" << "p." << name
    << ", " << g.constant(ipmc_blocks_pack(blocks)) << ");\n";
  }

  void unpack_block(const std::vector<casadi_int>& p, std::vector<casadi_ocp_block>& blocks) {
    const casadi_int* packed = get_ptr(p);
    casadi_int N = *packed++;
    blocks.resize(N);
    for (casadi_int i=0;i<N;++i) {
        blocks[i].offset_r = *packed++;
        blocks[i].offset_c = *packed++;
        blocks[i].rows = *packed++;
        blocks[i].cols = *packed++;
    }
  }

void IpmcInterface::set_ipmc_prob(CodeGenerator& g) const {
  if (jacg_sp_.size1()>0 && jacg_sp_.nnz()==0) {
    casadi_error("Empty sparsity pattern not supported in IPMC C interface");
  }
  g << "d->nlp = &d_nlp;\n";
  g << "d->rewrite.user = &d_nlp;\n";
  g << "d->prob = &p;\n";
  if (n_lift_>0) {
    // the REWRITTEN problem, and its own z/lbz/ubz/lam (mirror of
    // IpmcInterface::set_work; the w consumption must match exactly)
    g.local("p_nlp_lift", "struct casadi_nlpsol_prob");
    g.local("d_nlp_lift", "struct casadi_nlpsol_data");
    g << "p_nlp_lift.nx = " << nxt_ << ";\n";
    g << "p_nlp_lift.ng = " << nat_ << ";\n";
    g << "p_nlp_lift.np = " << np_ << ";\n";
    g << "p_nlp_lift.detect_bounds.ng = 0;\n";
    /* Exact mirror of IpmcInterface::set_work: copy the whole struct, then
       override what differs.  A bare local with only .prob set would leave
       every other field -- p among them -- as stack garbage the solve loop
       reads through d->nlp on every oracle call. */
    g << "d_nlp_lift = d_nlp;\n";
    g << "d_nlp_lift.prob = &p_nlp_lift;\n";
    g << "casadi_nlpsol_set_work(&d_nlp_lift, &arg, &res, &iw, &w);\n";
    g << "d->nlp = &d_nlp_lift;\n";
    g << "p.nlp = &p_nlp_lift;\n";
  } else {
    g << "p.nlp = &p_nlp;\n";
  }

  g << "p.nx = " << g.constant(nxs_) << ";\n";
  g << "p.nu = " << g.constant(nus_) << ";\n";
  g << "p.nxc = " << nxc_ + n_lift_ << ";\n";
  g << "p.ABsp = " << g.sparsity(ABsp_) << ";\n";
  g << "p.AB_offsets = " << g.constant(AB_offsets_) << ";\n";
  g << "p.CDsp = " << g.sparsity(CDsp_) << ";\n";
  g << "p.CD_offsets = " << g.constant(CD_offsets_) << ";\n";
  g << "p.RSQsp = " << g.sparsity(RSQsp_) << ";\n";
  g << "p.RSQ_offsets = " << g.constant(RSQ_offsets_) << ";\n";
  g << "p.Isp = " << g.sparsity(Isp_) << ";\n";
  g << "p.I_offsets = " << g.constant(I_offsets_) << ";\n";

  codegen_unpack_block(g, "AB", AB_blocks_);
  codegen_unpack_block(g, "CD", CD_blocks_);
  codegen_unpack_block(g, "RSQ", RSQ_blocks_);
  codegen_unpack_block(g, "I", I_blocks_);
  g << "p.N = " << N_ << ";\n";

  g << "p.ns = " << (slacks_ ? Nlpsol::ns_ : 0) << ";\n";
  if (slacks_) {
    g << "p.slack_g = " << g.constant(slack_g_) << ";\n";
    g << "p.slack_x = " << g.constant(slack_x_) << ";\n";
    g << "p.slack_perm = " << g.constant(slack_perm_) << ";\n";
    g << "p.slack_offs = " << g.constant(slack_offs_) << ";\n";
    g << "p.fs_z = " << g.constant(fs_z_) << ";\n";
    g << "p.fs_Z = " << g.constant(fs_Z_) << ";\n";
  } else {
    g << "p.slack_g = 0;\n";
    g << "p.slack_x = 0;\n";
    g << "p.slack_perm = 0;\n";
    g << "p.slack_offs = 0;\n";
    g << "p.fs_z = 0;\n";
    g << "p.fs_Z = 0;\n";
  }

  g << "p.rewrite.ne = " << n_lift_ << ";\n";
  g << "p.rewrite.nx = " << nx_ << ";\n";
  g << "p.rewrite.ng = " << ng_ << ";\n";
  g << "p.rewrite.K = " << N_+1 << ";\n";
  g << "p.rewrite.x = " << (lift_x_.empty() ? "0" : g.constant(lift_x_)) << ";\n";
  g << "p.rewrite.xfree = " << (lift_xfree_.empty() ? "0" : g.constant(lift_xfree_)) << ";\n";
  g << "p.rewrite.row = " << (lift_row_.empty() ? "0" : g.constant(lift_row_)) << ";\n";
  g << "p.rewrite.kind = " << (lift_kind_.empty() ? "0" : g.constant(lift_kind_)) << ";\n";
  g << "p.rewrite.ent = " << (lift_ent_.empty() ? "0" : g.constant(lift_ent_)) << ";\n";
  g << "p.rewrite.mflat = " << (lift_mflat_.empty() ? "0" : g.constant(lift_mflat_)) << ";\n";
  g << "p.rewrite.m = " << (lift_m_.empty() ? "0" : g.constant(lift_m_)) << ";\n";
  g << "p.rewrite.slot = " << (lift_slot_.empty() ? "0" : g.constant(lift_slot_)) << ";\n";
  g << "p.rewrite.a_src = " << (lift_a_src_.empty() ? "0" : g.constant(lift_a_src_)) << ";\n";
  g << "p.rewrite.h_src = " << (lift_h_src_.empty() ? "0" : g.constant(lift_h_src_)) << ";\n";
  g << "p.rewrite.nnz_au = " << jacg_sp_.nnz() << ";\n";
  g << "p.rewrite.nnz_hu = " << (exact_hessian_ ? hesslag_sp_.nnz() : 0) << ";\n";

  g << "p.sp_a = " << g.sparsity(n_lift_>0 ? jacg_lift_sp_ : jacg_sp_) << ";\n";
  if (exact_hessian_) {
    g << "p.sp_h = " << g.sparsity(n_lift_>0 ? hess_lift_sp_ : hesslag_sp_) << ";\n";
  } else {
    g << "p.sp_h = 0;\n";
  }

  g << "p.write = &" << g.shorthand("ipmc_cb_write") << ";\n";
  g << "p.flush = &" << g.shorthand("ipmc_cb_flush") << ";\n";

  g << "casadi_ipmc_setup(&p);\n";

}

IpmcInterface::IpmcInterface(DeserializingStream& s) : Nlpsol(s) {
  s.version("IpmcInterface", 1);
  s.unpack("IpmcInterface::jacg_sp", jacg_sp_);
  s.unpack("IpmcInterface::hesslag_sp", hesslag_sp_);
  s.unpack("IpmcInterface::exact_hessian", exact_hessian_);
  s.unpack("IpmcInterface::opts", opts_);
  s.unpack("IpmcInterface::convexify", convexify_);

  s.unpack("IpmcInterface::Isp", Isp_);
  s.unpack("IpmcInterface::ABsp", ABsp_);
  s.unpack("IpmcInterface::CDsp", CDsp_);
  s.unpack("IpmcInterface::RSQsp", RSQsp_);

  std::vector<casadi_int> AB_blocks;
  s.unpack("IpmcInterface::AB_blocks", AB_blocks);
  unpack_block(AB_blocks, AB_blocks_);
  std::vector<casadi_int> CD_blocks;
  s.unpack("IpmcInterface::CD_blocks", CD_blocks);
  unpack_block(CD_blocks, CD_blocks_);
  std::vector<casadi_int> RSQ_blocks;
  s.unpack("IpmcInterface::RSQ_blocks", RSQ_blocks);
  unpack_block(RSQ_blocks, RSQ_blocks_);
  std::vector<casadi_int> I_blocks;
  s.unpack("IpmcInterface::I_blocks", I_blocks);
  unpack_block(I_blocks, I_blocks_);

  s.unpack("IpmcInterface::nxs", nxs_);
  s.unpack("IpmcInterface::nus", nus_);
  s.unpack("IpmcInterface::ngs", ngs_);
  s.unpack("IpmcInterface::nxc", nxc_);
  s.unpack("IpmcInterface::N", N_);

  casadi_int structure_detection;
  s.unpack("IpmcInterface::structure_detection", structure_detection);
  structure_detection_ = static_cast<StructureDetection>(structure_detection);

  s.unpack("IpmcInterface::AB_offsets", AB_offsets_);
  s.unpack("IpmcInterface::CD_offsets", CD_offsets_);
  s.unpack("IpmcInterface::RSQ_offsets", RSQ_offsets_);
  s.unpack("IpmcInterface::I_offsets", I_offsets_);
  s.unpack("IpmcInterface::debug", debug_);

  // Native soft constraints
  s.unpack("IpmcInterface::slacks", slacks_);
  if (slacks_) {
    s.unpack("IpmcInterface::slack_g", slack_g_);
    s.unpack("IpmcInterface::slack_x", slack_x_);
    s.unpack("IpmcInterface::slack_stage", slack_stage_);
    s.unpack("IpmcInterface::slack_perm", slack_perm_);
    s.unpack("IpmcInterface::slack_ns", slack_ns_);
    s.unpack("IpmcInterface::slack_offs", slack_offs_);
    s.unpack("IpmcInterface::fs_z", fs_z_);
    s.unpack("IpmcInterface::fs_Z", fs_Z_);
  }

  // The cross-stage lifting.  Without it the caller's partition is the one
  // that was handed over.
  nxt_ = nx_;
  nat_ = ng_;
  nxs_u_ = nxs_;
  ngs_u_ = ngs_;
  s.unpack("IpmcInterface::n_lift", n_lift_);
  if (n_lift_>0) {
    s.unpack("IpmcInterface::nxt", nxt_);
    s.unpack("IpmcInterface::nat", nat_);
    s.unpack("IpmcInterface::nxs_u", nxs_u_);
    s.unpack("IpmcInterface::ngs_u", ngs_u_);
    s.unpack("IpmcInterface::lift_x", lift_x_);
    s.unpack("IpmcInterface::lift_xfree", lift_xfree_);
    s.unpack("IpmcInterface::lift_row", lift_row_);
    s.unpack("IpmcInterface::lift_kind", lift_kind_);
    s.unpack("IpmcInterface::lift_ent", lift_ent_);
    s.unpack("IpmcInterface::lift_mflat", lift_mflat_);
    s.unpack("IpmcInterface::lift_m", lift_m_);
    s.unpack("IpmcInterface::lift_slot", lift_slot_);
    /* the patterns and nonzero maps are deterministic in what the stream
       carries, so they are rebuilt rather than stored */
    lift_patterns();
  }

  set_ipmc_prob();
}

void IpmcInterface::serialize_body(SerializingStream &s) const {
  Nlpsol::serialize_body(s);
  s.version("IpmcInterface", 1);

  s.pack("IpmcInterface::jacg_sp", jacg_sp_);
  s.pack("IpmcInterface::hesslag_sp", hesslag_sp_);
  s.pack("IpmcInterface::exact_hessian", exact_hessian_);
  s.pack("IpmcInterface::opts", opts_);
  s.pack("IpmcInterface::convexify", convexify_);

  s.pack("IpmcInterface::Isp", Isp_);
  s.pack("IpmcInterface::ABsp", ABsp_);
  s.pack("IpmcInterface::CDsp", CDsp_);
  s.pack("IpmcInterface::RSQsp", RSQsp_);

  s.pack("IpmcInterface::AB_blocks", ipmc_blocks_pack(AB_blocks_));
  s.pack("IpmcInterface::CD_blocks", ipmc_blocks_pack(CD_blocks_));
  s.pack("IpmcInterface::RSQ_blocks", ipmc_blocks_pack(RSQ_blocks_));
  s.pack("IpmcInterface::I_blocks", ipmc_blocks_pack(I_blocks_));

  s.pack("IpmcInterface::nxs", nxs_);
  s.pack("IpmcInterface::nus", nus_);
  s.pack("IpmcInterface::ngs", ngs_);
  s.pack("IpmcInterface::nxc", nxc_);
  s.pack("IpmcInterface::N", N_);
  s.pack("IpmcInterface::structure_detection", static_cast<casadi_int>(structure_detection_));
  s.pack("IpmcInterface::AB_offsets", AB_offsets_);
  s.pack("IpmcInterface::CD_offsets", CD_offsets_);
  s.pack("IpmcInterface::RSQ_offsets", RSQ_offsets_);
  s.pack("IpmcInterface::I_offsets", I_offsets_);
  s.pack("IpmcInterface::debug", debug_);

  s.pack("IpmcInterface::slacks", slacks_);
  if (slacks_) {
    s.pack("IpmcInterface::slack_g", slack_g_);
    s.pack("IpmcInterface::slack_x", slack_x_);
    s.pack("IpmcInterface::slack_stage", slack_stage_);
    s.pack("IpmcInterface::slack_perm", slack_perm_);
    s.pack("IpmcInterface::slack_ns", slack_ns_);
    s.pack("IpmcInterface::slack_offs", slack_offs_);
    s.pack("IpmcInterface::fs_z", fs_z_);
    s.pack("IpmcInterface::fs_Z", fs_Z_);
  }

  s.pack("IpmcInterface::n_lift", n_lift_);
  if (n_lift_>0) {
    s.pack("IpmcInterface::nxt", nxt_);
    s.pack("IpmcInterface::nat", nat_);
    s.pack("IpmcInterface::nxs_u", nxs_u_);
    s.pack("IpmcInterface::ngs_u", ngs_u_);
    s.pack("IpmcInterface::lift_x", lift_x_);
    s.pack("IpmcInterface::lift_xfree", lift_xfree_);
    s.pack("IpmcInterface::lift_row", lift_row_);
    s.pack("IpmcInterface::lift_kind", lift_kind_);
    s.pack("IpmcInterface::lift_ent", lift_ent_);
    s.pack("IpmcInterface::lift_mflat", lift_mflat_);
    s.pack("IpmcInterface::lift_m", lift_m_);
    s.pack("IpmcInterface::lift_slot", lift_slot_);
  }
}

} // namespace casadi
