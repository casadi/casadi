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

#include "hpipm_interface.hpp"
#include <numeric>
#include <cstring>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_HPIPM_EXPORT
  casadi_register_conic_hpipm(Conic::Plugin* plugin) {
    plugin->creator = HpipmInterface::creator;
    plugin->name = "hpipm";
    plugin->doc = HpipmInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &HpipmInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_HPIPM_EXPORT casadi_load_conic_hpipm() {
    Conic::registerPlugin(casadi_register_conic_hpipm);
  }

  HpipmInterface::HpipmInterface(const std::string& name,
                                     const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
  }

  HpipmInterface::~HpipmInterface() {
    clear_mem();
  }

  const Options HpipmInterface::options_
  = {{&Conic::options_},
     {{"N",
       {OT_INT,
        "OCP horizon"}},
      {"nx",
       {OT_INTVECTOR,
        "Number of states, length N+1"}},
      {"nu",
       {OT_INTVECTOR,
        "Number of controls, length N"}},
      {"ng",
       {OT_INTVECTOR,
        "Number of non-dynamic constraints, length N+1"}},
      {"inf",
       {OT_DOUBLE,
        "Replace infinities by this amount [default: 1e8]"}},
      {"hpipm",
       {OT_DICT,
        "Options to be passed to hpipm"}}}
  };

  void HpipmInterface::init(const Dict& opts) {
    Conic::init(opts);

    hpipm_mode mode = ROBUST;
    casadi_int struct_cnt=0;

    d_ocp_qp_ipm_arg_set_default(mode, &hpipm_options_);

    inf_ = 1e8;

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
      } else if (op.first=="inf") {
        inf_ = op.second;
      } else if (op.first=="hpipm") {
        Dict hopts = op.second;
        auto it = hopts.find("mode");
        if (it!=hopts.end()) {
          if (it->second=="speed_abs") {
            mode = SPEED_ABS;
          } else if (it->second=="speed") {
            mode = SPEED;
          } else if (it->second=="balance") {
            mode = BALANCE;
          } else if (it->second=="robust") {
            mode = ROBUST;
          } else {
            casadi_error("Unknown mode. Choose from speed_abs, speedm balance, robust.");
          }
        }

        d_ocp_qp_ipm_arg_set_default(mode, &hpipm_options_);

        for (auto&& op : hopts) {
          if (op.first=="mu0") {
            hpipm_options_.mu0 = op.second;
          } else if (op.first=="alpha_min") {
            hpipm_options_.alpha_min = op.second;
          } else if (op.first=="res_g_max") {
            hpipm_options_.res_g_max = op.second;
          } else if (op.first=="res_b_max") {
            hpipm_options_.res_b_max = op.second;
          } else if (op.first=="res_d_max") {
            hpipm_options_.res_d_max = op.second;
          } else if (op.first=="res_m_max") {
            hpipm_options_.res_m_max = op.second;
          } else if (op.first=="iter_max") {
            hpipm_options_.iter_max = op.second;
          } else if (op.first=="stat_max") {
            hpipm_options_.stat_max = op.second;
          } else if (op.first=="pred_corr") {
            hpipm_options_.pred_corr = op.second;
          } else if (op.first=="cond_pred_corr") {
            hpipm_options_.cond_pred_corr = op.second;
          } else if (op.first=="itref_pred_max") {
            hpipm_options_.itref_pred_max = op.second;
          } else if (op.first=="itref_corr_max") {
            hpipm_options_.itref_corr_max = op.second;
          } else if (op.first=="reg_prim") {
            hpipm_options_.reg_prim = op.second;
          } else if (op.first=="lq_fact") {
            hpipm_options_.lq_fact = op.second;
          } else if (op.first=="lam_min") {
            hpipm_options_.lam_min = op.second;
          } else if (op.first=="t_min") {
            hpipm_options_.t_min = op.second;
          } else if (op.first=="warm_start") {
            hpipm_options_.warm_start = op.second;
          } else if (op.first=="abs_form") {
            hpipm_options_.abs_form = op.second;
          } else if (op.first=="comp_dual_sol_eq") {
            hpipm_options_.comp_dual_sol_eq = op.second;
          } else if (op.first=="comp_res_exit") {
            hpipm_options_.comp_res_exit = op.second;
          } else if (op.first=="mode") {
            // already covered
          } else {
            casadi_error("Unknown option.");
          }
        }
      }
    }

    bool detect_structure = struct_cnt==0;
    casadi_assert(struct_cnt==0 || struct_cnt==4,
      "You must either set all of N, nx, nu, ng; "
      "or set none at all (automatic detection).");

    const std::vector<casadi_int>& nx = nxs_;
    const std::vector<casadi_int>& ng = ngs_;
    const std::vector<casadi_int>& nu = nus_;

    if (detect_structure) {
      /* General strategy: look for the xk+1 diagonal part in A
      */

      // Find the right-most column for each row in A -> A_skyline
      // Find the second-to-right-most column -> A_skyline2
      // Find the left-most column -> A_bottomline
      Sparsity AT = A_.T();
      std::vector<casadi_int> A_skyline;
      std::vector<casadi_int> A_skyline2;
      std::vector<casadi_int> A_bottomline;
      for (casadi_int i=0;i<AT.size2();++i) {
        casadi_int pivot = AT.colind()[i+1];
        A_bottomline.push_back(AT.row()[AT.colind()[i]]);
        if (pivot>AT.colind()[i]) {
          A_skyline.push_back(AT.row()[pivot-1]);
          if (pivot>AT.colind()[i]+1) {
            A_skyline2.push_back(AT.row()[pivot-2]);
          } else {
            A_skyline2.push_back(-1);
          }
        } else {
          A_skyline.push_back(-1);
          A_skyline2.push_back(-1);
        }
      }

      /*
      Loop over the right-most columns of A:
      they form the diagonal part due to xk+1 in gap constraints.
      detect when the diagonal pattern is broken -> new stage
      */
      casadi_int pivot = 0; // Current right-most element
      casadi_int start_pivot = pivot; // First right-most element that started the stage
      casadi_int cg = 0; // Counter for non-gap-closing constraints
      for (casadi_int i=0;i<na_;++i) { // Loop over all rows
        bool commit = false; // Set true to jump to the stage
        if (A_skyline[i]>pivot+1) { // Jump to a diagonal in the future
          nus_.push_back(A_skyline[i]-pivot-1); // Size of jump equals number of states
          commit = true;
        } else if (A_skyline[i]==pivot+1) { // Walking the diagonal
          if (A_skyline2[i]<start_pivot) { // Free of below-diagonal entries?
            pivot++;
          } else {
            nus_.push_back(0); // We cannot but conclude that we arrived at a new stage
            commit = true;
          }
        } else { // non-gap-closing constraint detected
          cg++;
        }

        if (commit) {
          nxs_.push_back(pivot-start_pivot+1);
          ngs_.push_back(cg); cg=0;
          start_pivot = A_skyline[i];
          pivot = A_skyline[i];
        }
      }
      nxs_.push_back(pivot-start_pivot+1);

      uout() << nx << nu << ng << std::endl;
      // Correction for k==0
      nxs_[0] = A_skyline[0];
      nus_[0] = 0;
      ngs_.erase(ngs_.begin());
      casadi_int cN=0;
      for (casadi_int i=na_-1;i>=0;--i) {
        if (A_bottomline[i]<start_pivot) break;
        cN++;
      }
      ngs_.push_back(cg-cN);
      ngs_.push_back(cN);

      N_ = nus_.size();
      if (verbose_) {
        casadi_message("Detected structure: N " + str(N_) + ", nx " + str(nx) + ", "
          "nu " + str(nu) + ", ng " + str(ng) + ".");
      }
      nus_.push_back(0);
    }

    uout() << nx << nu << ng << std::endl;

    casadi_assert_dev(nx.size()==N_+1);
    casadi_assert_dev(nu.size()==N_+1);
    casadi_assert_dev(ng.size()==N_+1);

    casadi_assert(nx_ == std::accumulate(nx.begin(), nx.end(), 0) + // NOLINT
      std::accumulate(nu.begin(), nu.end(), 0),
      "sum(nx)+sum(nu) = must equal total size of variables (" + str(nx_) + "). "
      "Structure is: N " + str(N_) + ", nx " + str(nx) + ", "
      "nu " + str(nu) + ", ng " + str(ng) + ".");
    casadi_assert(na_ == std::accumulate(nx.begin()+1, nx.end(), 0) + // NOLINT
      std::accumulate(ng.begin(), ng.end(), 0),
      "sum(nx+1)+sum(ng) = must equal total size of constraints (" + str(na_) + "). "
      "Structure is: N " + str(N_) + ", nx " + str(nx) + ", "
      "nu " + str(nu) + ", ng " + str(ng) + ".");
    // Load library HPIPM when applicable
    std::string searchpath;

    /* Disassemble A input into:
       A B I
       C D
           A B I
           C D
    */
    casadi_int offset_r = 0, offset_c = 0;
    for (casadi_int k=0;k<N_;++k) { // Loop over blocks
      A_blocks.push_back({offset_r,        offset_c,            nx[k+1], nx[k]});
      B_blocks.push_back({offset_r,        offset_c+nx[k],      nx[k+1], nu[k]});
      C_blocks.push_back({offset_r+nx[k+1], offset_c,           ng[k], nx[k]});
      D_blocks.push_back({offset_r+nx[k+1], offset_c+nx[k],     ng[k], nu[k]});

      offset_c+= nx[k]+nu[k];
      if (k+1<N_)
        I_blocks.push_back({offset_r, offset_c, nx[k+1], nx[k+1]});
      else
        I_blocks.push_back({offset_r, offset_c, nx[k+1], nx[k+1]});
      offset_r+= nx[k+1]+ng[k];
    }

    C_blocks.push_back({offset_r, offset_c,            ng[N_], nx[N_]});
    D_blocks.push_back({offset_r, offset_c+nx[N_],     ng[N_], nu[N_]});

    Asp_ = blocksparsity(na_, nx_, A_blocks);
    Bsp_ = blocksparsity(na_, nx_, B_blocks);
    Csp_ = blocksparsity(na_, nx_, C_blocks);
    Dsp_ = blocksparsity(na_, nx_, D_blocks);
    Isp_ = blocksparsity(na_, nx_, I_blocks, true);

    Sparsity total = Asp_ + Bsp_ + Csp_ + Dsp_ + Isp_;

    casadi_assert((A_ + total).nnz() == total.nnz(),
      "HPIPM: specified structure of A does not correspond to what the interface can handle. "
      "Structure is: N " + str(N_) + ", nx " + str(nx) + ", nu " + str(nu) + ", "
      "ng " + str(ng) + ".");
    casadi_assert_dev(total.nnz() == Asp_.nnz() + Bsp_.nnz() + Csp_.nnz() + Dsp_.nnz()
                      + Isp_.nnz());

    /* Disassemble H input into:
       Q S'
       S R
           Q S'
           S R

       Multiply by 2
    */
    casadi_int offset = 0;
    for (casadi_int k=0;k<N_+1;++k) { // Loop over blocks
      R_blocks.push_back({offset+nx[k], offset+nx[k],       nu[k], nu[k]});
      S_blocks.push_back({offset+nx[k], offset,             nu[k], nx[k]});
      Q_blocks.push_back({offset,       offset,             nx[k], nx[k]});
      offset+= nx[k]+nu[k];
    }

    Rsp_ = blocksparsity(nx_, nx_, R_blocks);
    Ssp_ = blocksparsity(nx_, nx_, S_blocks);
    Qsp_ = blocksparsity(nx_, nx_, Q_blocks);

    total = Rsp_ + Ssp_ + Qsp_ + Ssp_.T();
    casadi_assert((H_ + total).nnz() == total.nnz(),
      "HPIPM: specified structure of H does not correspond to what the interface can handle. "
      "Structure is: N " + str(N_) + ", nx " + str(nx) + ", nu " + str(nu) + ", "
      "ng " + str(ng) + ".");
    casadi_assert_dev(total.nnz() == Rsp_.nnz() + 2*Ssp_.nnz() + Qsp_.nnz());

    /* Disassemble LBA/UBA input into:
       b
       lg/ug

       b
       lg/ug
    */
    offset = 0;

    for (casadi_int k=0;k<N_;++k) {
      b_blocks.push_back({offset,   0, nx[k+1], 1}); offset+= nx[k+1];
      lug_blocks.push_back({offset, 0, ng[k], 1}); offset+= ng[k];
    }
    lug_blocks.push_back({offset, 0, ng[N_], 1});

    bsp_ = blocksparsity(na_, 1, b_blocks);
    lugsp_ = blocksparsity(na_, 1, lug_blocks);
    total = bsp_ + lugsp_;
    casadi_assert_dev(total.nnz() == bsp_.nnz() + lugsp_.nnz());
    casadi_assert_dev(total.nnz() == na_);

    /* Disassemble G/X0 input into:
       r/u
       q/x

       r/u
       q/x
    */
    offset = 0;

    for (casadi_int k=0;k<N_+1;++k) {
      x_blocks.push_back({offset, 0, nx[k], 1}); offset+= nx[k];
      u_blocks.push_back({offset, 0, nu[k], 1}); offset+= nu[k];
    }

    usp_ = blocksparsity(nx_, 1, u_blocks);
    xsp_ = blocksparsity(nx_, 1, x_blocks);
    total = usp_ + xsp_;
    casadi_assert_dev(total.nnz() == usp_.nnz() + xsp_.nnz());
    casadi_assert_dev(total.nnz() == nx_);

    std::vector< Block > theirs_u_blocks, theirs_x_blocks;
    offset = 0;

    for (casadi_int k=0;k<N_;++k) {
      theirs_u_blocks.push_back({offset, 0, nu[k], 1}); offset+= nu[k];
      theirs_x_blocks.push_back({offset, 0, nx[k], 1}); offset+= nx[k];
    }
    theirs_x_blocks.push_back({offset, 0, nx[N_], 1});

    theirs_usp_ = blocksparsity(nx_, 1, theirs_u_blocks);
    theirs_xsp_ = blocksparsity(nx_, 1, theirs_x_blocks);
    total = theirs_usp_ + theirs_xsp_;
    casadi_assert_dev(total.nnz() == theirs_usp_.nnz() + theirs_xsp_.nnz());
    casadi_assert_dev(total.nnz() == nx_);

    offset = 0;
    std::vector< Block > lamg_gap_blocks;
    for (casadi_int k=0;k<N_;++k) {
      lamg_gap_blocks.push_back({offset,       0, nx[k+1], 1});offset+= nx[k+1] + ng[k];
    }
    lamg_gapsp_ = blocksparsity(na_, 1, lamg_gap_blocks);
    lamg_csp_ = lamg_gapsp_.pattern_inverse();

    offset = 0;

    for (casadi_int k=0;k<N_;++k) {
      lam_ul_blocks.push_back({offset, 0, nu[k], 1}); offset+= nu[k];
      lam_xl_blocks.push_back({offset, 0, nx[k], 1}); offset+= nx[k];
      lam_uu_blocks.push_back({offset, 0, nu[k], 1}); offset+= nu[k];
      lam_xu_blocks.push_back({offset, 0, nx[k], 1}); offset+= nx[k];
      lam_cl_blocks.push_back({offset, 0, ng[k], 1}); offset+= ng[k];
      lam_cu_blocks.push_back({offset, 0, ng[k], 1}); offset+= ng[k];
    }
    lam_xl_blocks.push_back({offset, 0, nx[N_], 1}); offset+= nx[N_];
    lam_xu_blocks.push_back({offset, 0, nx[N_], 1}); offset+= nx[N_];
    lam_cl_blocks.push_back({offset, 0, ng[N_], 1}); offset+= ng[N_];
    lam_cu_blocks.push_back({offset, 0, ng[N_], 1}); offset+= ng[N_];

    lam_ulsp_ = blocksparsity(offset, 1, lam_ul_blocks);
    lam_uusp_ = blocksparsity(offset, 1, lam_uu_blocks);
    lam_xlsp_ = blocksparsity(offset, 1, lam_xl_blocks);
    lam_xusp_ = blocksparsity(offset, 1, lam_xu_blocks);
    lam_clsp_ = blocksparsity(offset, 1, lam_cl_blocks);
    lam_cusp_ = blocksparsity(offset, 1, lam_cu_blocks);

    pisp_ = Sparsity::dense(std::accumulate(nx.begin()+1, nx.end(), 0), 1);  // NOLINT

    total = lam_ulsp_ + lam_uusp_ + lam_xlsp_ + lam_xusp_ + lam_clsp_ + lam_cusp_;
    casadi_assert_dev(total.nnz() == lam_ulsp_.nnz() + lam_uusp_.nnz() + lam_xlsp_.nnz() +
      lam_xusp_.nnz() + lam_clsp_.nnz() + lam_cusp_.nnz());
    casadi_assert_dev(total.nnz() == offset);

    theirs_Xsp_ = Sparsity::dense(std::accumulate(nx.begin(), nx.end(), 0), 1);  // NOLINT
    theirs_Usp_ = Sparsity::dense(std::accumulate(nu.begin(), nu.end(), 0), 1);  // NOLINT

  }


  int HpipmInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<HpipmMemory*>(mem);

    init_vector(m->nx, nxs_);
    init_vector(m->nu, nus_); m->nu.push_back(0);
    init_vector(m->ng, ngs_);
    m->ns.resize(N_+1, 0);
    m->nsbx.resize(N_+1, 0);
    m->nsbu.resize(N_+1, 0);
    m->nsg.resize(N_+1, 0);

    const std::vector<int>& nx = m->nx;
    const std::vector<int>& nu = m->nu;
    const std::vector<int>& ng = m->ng;
    const std::vector<int>& nbx = m->nbx;
    const std::vector<int>& nbu = m->nbu;
    
    casadi_int offset = 0;

    m->nbx.resize(N_+1);
    for (casadi_int k=0;k<N_+1;++k) m->nbx[k] = nx[k];

    m->nbu.resize(N_+1);
    for (casadi_int k=0;k<N_+1;++k) m->nbu[k] = nu[k];

    uout() <<  "nbx;nbu" << m->nbx <<  m->nbu << std::endl;
    uout() <<  "nx;nu" << m->nx <<  m->nu << std::endl;

    m->A.resize(Asp_.nnz());
    m->B.resize(Bsp_.nnz());
    m->C.resize(Csp_.nnz());
    m->D.resize(Dsp_.nnz());
    m->R.resize(Rsp_.nnz());
    m->I.resize(Isp_.nnz());
    m->S.resize(Ssp_.nnz());
    m->Q.resize(Qsp_.nnz());
    m->b.resize(bsp_.nnz());
    m->b2.resize(bsp_.nnz());
    m->x.resize(xsp_.nnz());m->q.resize(xsp_.nnz());
    m->u.resize(usp_.nnz());m->r.resize(usp_.nnz());
    m->lg.resize(lugsp_.nnz());
    m->ug.resize(lugsp_.nnz());

    m->pi.resize(pisp_.nnz());

    casadi_int nx_total = 0;
    for (casadi_int k=0;k<N_+1;++k) nx_total+=nx[k];
    m->iidxbx.resize(nx_total); // TODO; make range?
    m->lbx.resize(nx_total);
    m->ubx.resize(nx_total);

    casadi_int nu_total = 0;
    for (casadi_int k=0;k<N_+1;++k) nu_total+=nu[k];
    m->iidxbu.resize(nu_total); // TODO; make range?

    m->lbu.resize(nu_total);
    m->ubu.resize(nu_total);

    offset = 0;
    for (casadi_int k=0;k<N_+1;++k) offset+=ng[k]+nx[k];
    m->lam.resize(2*offset);

    // Allocate double* work vectors
    blockptr(m->hA, m->A, A_blocks);
    blockptr(m->hB, m->B, B_blocks);
    blockptr(m->hC, m->C, C_blocks);
    blockptr(m->hD, m->D, D_blocks);
    blockptr(m->hI, m->I, I_blocks, true);
    blockptr(m->hR, m->R, R_blocks);
    blockptr(m->hQ, m->Q, Q_blocks);
    blockptr(m->hS, m->S, S_blocks);
    blockptr(m->hu_guess, m->u, u_blocks);
    blockptr(m->hx_guess, m->x, x_blocks);
    blockptr(m->hr, m->r, u_blocks);
    blockptr(m->hq, m->q, x_blocks);
    blockptr(m->hlg, m->lg, lug_blocks);
    blockptr(m->hug, m->ug, lug_blocks);
    blockptr(m->hb, m->b, b_blocks);

    m->pis.resize(N_);
    offset = 0;
    for (casadi_int k=0;k<N_;++k) {
      m->pis[k] = get_ptr(m->pi)+offset;
      offset+=nx[k+1];
    }

    m->hlbx.resize(N_+1);
    m->hubx.resize(N_+1);
    offset = 0;
    for (casadi_int k=0;k<N_+1;++k) {
      m->hlbx[k] = get_ptr(m->lbx)+offset;
      m->hubx[k] = get_ptr(m->ubx)+offset;
      offset+=nx[k];
    }
    m->hlbu.resize(N_+1);
    m->hubu.resize(N_+1);
    offset = 0;
    for (casadi_int k=0;k<N_+1;++k) {
      m->hlbu[k] = get_ptr(m->lbu)+offset;
      m->hubu[k] = get_ptr(m->ubu)+offset;
      offset+=nu[k];
    }

    m->lams.resize(N_+1);
    offset = 0;
    for (casadi_int k=0;k<N_+1;++k) {
      m->lams[k] = get_ptr(m->lam)+offset;
      offset+=2*(ng[k]+nx[k]);
    }

    m->hidxbx.resize(N_+1);
    m->hidxbu.resize(N_+1);
    offset = 0;
    for (casadi_int k=0;k<N_+1;++k) {
      m->hidxbx[k] = get_ptr(m->iidxbx)+offset;
      for (casadi_int i=0;i<nbx[k];++i) m->hidxbx[k][i] = i;
      offset+=nbx[k];
    }
    offset = 0;
    for (casadi_int k=0;k<N_+1;++k) {
      m->hidxbu[k] = get_ptr(m->iidxbu)+offset;
      for (casadi_int i=0;i<nbx[k];++i) m->hidxbu[k][i] = i;
      offset+=nbu[k];
    }

    m->pv.resize(2*(nx_+na_));

    // We don't interface the slack-variable feature
    m->hZl.resize(N_+1, nullptr);
    m->hZu.resize(N_+1, nullptr);
    m->hzl.resize(N_+1, nullptr);
    m->hzu.resize(N_+1, nullptr);
    m->hlls.resize(N_+1, nullptr);
    m->hlus.resize(N_+1, nullptr);
    m->hidxs.resize(N_+1, 0);

    m->fstats["preprocessing"]  = FStats();
    m->fstats["solver"]         = FStats();
    m->fstats["postprocessing"] = FStats();
    return 0;
  }

  void HpipmInterface::mproject(double factor, const double* x, const casadi_int* sp_x,
                                double* y, const casadi_int* sp_y, double* w) {
    casadi_int ncol_y = sp_y[1];
    const casadi_int *colind_y = sp_y+2;
    casadi_project(x, sp_x, y, sp_y, w);
    casadi_scal(colind_y[ncol_y], factor, y);
  }

  void HpipmInterface::dense_transfer(double factor, const double* x,
                                      const casadi_int* sp_x, double* y,
                                      const casadi_int* sp_y, double* w) {
    CASADI_PREFIX(sparsify)(x, w, sp_x, false);
    casadi_int nrow_y = sp_y[0];
    casadi_int ncol_y = sp_y[1];
    const casadi_int *colind_y = sp_y+2, *row_y = sp_y + 2 + ncol_y+1;
    /* Loop over columns of y */
    casadi_int i, el;
    for (i=0; i<ncol_y; ++i) {
      for (el=colind_y[i]; el<colind_y[i+1]; ++el) y[nrow_y*i + row_y[el]] += factor*(*w++);
    }
  }

  int HpipmInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<HpipmMemory*>(mem);
    // Statistics
    m->fstats.at("preprocessing").tic();

    int dim_size = d_ocp_qp_dim_memsize(N_);

	  void *dim_mem = malloc(dim_size);

	  struct d_ocp_qp_dim dim;
	  d_ocp_qp_dim_create(N_, &dim, dim_mem);

    d_ocp_qp_dim_set_all(
      get_ptr(m->nx), get_ptr(m->nu), get_ptr(m->nbx), get_ptr(m->nbu),
      get_ptr(m->ng), get_ptr(m->nsbx), get_ptr(m->nsbu), get_ptr(m->nsg), &dim);

    int qp_size = d_ocp_qp_memsize(&dim);
    void *qp_mem = malloc(qp_size);

    struct d_ocp_qp qp;
    d_ocp_qp_create(&dim, &qp, qp_mem);


    double* pv =  get_ptr(m->pv);

    // Dissect A matrix
    casadi_project(arg[CONIC_A], A_, get_ptr(m->A), Asp_, pv);
    casadi_project(arg[CONIC_A], A_, get_ptr(m->B), Bsp_, pv);
    casadi_project(arg[CONIC_A], A_, get_ptr(m->C), Csp_, pv);
    casadi_project(arg[CONIC_A], A_, get_ptr(m->D), Dsp_, pv);
    casadi_project(arg[CONIC_A], A_, get_ptr(m->I), Isp_, pv);

    // Dissect H matrix; definition of HPIPM lacks a factor 2
    mproject(0.5, arg[CONIC_H], H_, get_ptr(m->R), Rsp_, pv);
    mproject(0.5, arg[CONIC_H], H_, get_ptr(m->S), Ssp_, pv);
    mproject(0.5, arg[CONIC_H], H_, get_ptr(m->Q), Qsp_, pv);

    // Dissect LBA/UBA
    mproject(-1.0, arg[CONIC_LBA], sparsity_in_.at(CONIC_LBA), get_ptr(m->b), bsp_, pv);
    mproject(-1.0, arg[CONIC_UBA], sparsity_in_.at(CONIC_UBA), get_ptr(m->b2), bsp_, pv);
    casadi_assert_dev(std::equal(m->b.begin(), m->b.end(), m->b2.begin()));
    casadi_project(arg[CONIC_LBA], sparsity_in_.at(CONIC_LBA), get_ptr(m->lg), lugsp_, pv);
    casadi_project(arg[CONIC_UBA], sparsity_in_.at(CONIC_UBA), get_ptr(m->ug), lugsp_, pv);

    // Dissect LBX/UBX input
    std::fill(m->lbu.begin(), m->lbu.end(), 0);
    std::fill(m->ubu.begin(), m->ubu.end(), 0);
    std::fill(m->lbx.begin(), m->lbx.end(), 0);
    std::fill(m->ubx.begin(), m->ubx.end(), 0);

    dense_transfer(1.0, arg[CONIC_LBX], xsp_, get_ptr(m->lbx), theirs_Xsp_, pv);
    dense_transfer(1.0, arg[CONIC_UBX], xsp_, get_ptr(m->ubx), theirs_Xsp_, pv);
    dense_transfer(1.0, arg[CONIC_LBX], usp_, get_ptr(m->lbu), theirs_Usp_, pv);
    dense_transfer(1.0, arg[CONIC_UBX], usp_, get_ptr(m->ubu), theirs_Usp_, pv);

    // Dissect G
    mproject(0.5, arg[CONIC_G], sparsity_in_.at(CONIC_G), get_ptr(m->r), usp_, pv);
    mproject(0.5, arg[CONIC_G], sparsity_in_.at(CONIC_G), get_ptr(m->q), xsp_, pv);

    // Dissect X0
    casadi_project(arg[CONIC_X0], sparsity_in_.at(CONIC_X0), get_ptr(m->u), usp_, pv);
    casadi_project(arg[CONIC_X0], sparsity_in_.at(CONIC_X0), get_ptr(m->x), xsp_, pv);

    m->iter_count = -1;

    // Deal with non-unity I block
    for (casadi_int k=0;k<N_;++k) {
      casadi_int n_row = m->nx[k+1];
      for (casadi_int i=0;i<n_row;++i) {
        double f = -1/m->hI[k][i];
        m->hb[k][i]*=f;
        for (casadi_int j=0;j<m->nx[k];++j) m->hA[k][i+j*n_row]*=f;
        for (casadi_int j=0;j<m->nu[k];++j) m->hB[k][i+j*n_row]*=f;
      }
    }

    // replace infinities
    for (casadi_int i=0;i<m->lbx.size();++i) {
      if (m->lbx[i]==-std::numeric_limits<double>::infinity()) m->lbx[i] = -inf_;
    }
    for (casadi_int i=0;i<m->lbu.size();++i) {
      if (m->lbu[i]==-std::numeric_limits<double>::infinity()) m->lbu[i] = -inf_;
    }
    for (casadi_int i=0;i<m->ubx.size();++i) {
      if (m->ubx[i]==std::numeric_limits<double>::infinity()) m->ubx[i] = inf_;
    }
    for (casadi_int i=0;i<m->ubu.size();++i) {
      if (m->ubu[i]==std::numeric_limits<double>::infinity()) m->ubu[i] = inf_;
    }
    for (casadi_int i=0;i<m->lg.size();++i) {
      if (m->lg[i]==-std::numeric_limits<double>::infinity()) m->lg[i] = -inf_;
    }
    for (casadi_int i=0;i<m->ug.size();++i) {
      if (m->ug[i]==std::numeric_limits<double>::infinity()) m->ug[i] = inf_;
    }

    m->fstats.at("preprocessing").toc();



    std::fill(m->pi.begin(), m->pi.end(), 0);
    std::fill(m->lam.begin(), m->lam.end(), 0);

/*
    if (arg[CONIC_LAM_A0]) {
      dense_transfer(0.5, arg[CONIC_LAM_A0], lamg_gapsp_, get_ptr(m->pi), pisp_, pv);
      // Deal with non-unity I block
      for (casadi_int k=0;k<N_;++k) {
        casadi_int n_row = m->nx[k+1];
        for (casadi_int i=0;i<n_row;++i) {
          double f = -m->hI[k][i];
          m->pis[k][i]*=f;
        }
      }

      dense_transfer(0.5, arg[CONIC_LAM_A0], lamg_csp_, get_ptr(m->lam), lam_cusp_, pv);
    }

    if (arg[CONIC_LAM_X0]) {
      dense_transfer(0.5, arg[CONIC_LAM_X0], usp_, get_ptr(m->lam), lam_uusp_, pv);
      dense_transfer(0.5, arg[CONIC_LAM_X0], xsp_, get_ptr(m->lam), lam_xusp_, pv);
    }*/

    /*
    uout() << m->hA << std::endl;
    uout() << m->hB << std::endl;
    uout() << m->hb << std::endl;
    uout() << m->hQ << std::endl;
    uout() << m->hR << std::endl;
    uout() << m->hS << std::endl;
    uout() << m->hq << std::endl;
    uout() << m->hr << std::endl;
    uout() << m->hlb << std::endl;
    uout() << m->hub << std::endl;
    uout() << m->hC << std::endl;
    uout() << m->hD << std::endl;
    uout() << m->hlg << std::endl;
    uout() << m->hug << std::endl;
    uout() << m->hZl << std::endl;
    uout() << m->hZu << std::endl;
    uout() << m->hzl << std::endl;
    uout() << m->hzu << std::endl;
    uout() << m->hlls << std::endl;
    uout() << m->hlus << std::endl;

    uout() << m->A << std::endl;
    uout() <<  m->B << std::endl;
    uout() <<  m->b << std::endl;
    uout() <<  m->b2 << std::endl;
    uout() <<  m->Q << std::endl;
    uout() <<  m->S << std::endl;
    uout() <<  m->R << std::endl;
    uout() << m->q << std::endl;
    uout() <<  m->r << std::endl;
    uout() <<  m->lb  << std::endl;
    uout() << m->ub << std::endl;
    uout() <<  m->C << std::endl;
    uout() <<  m->D << std::endl;
    uout() << m->lg << std::endl;
    uout() << m->ug << std::endl;
    */

    // https://github.com/giaf/hpipm/commit/590f3b21521ff5784c9497869883df695bc0f441
    std::vector<double*> hI;

    d_ocp_qp_set_all(get_ptr(m->hA), get_ptr(m->hB), get_ptr(m->hb),
      get_ptr(m->hQ), get_ptr(m->hS), get_ptr(m->hR),
      get_ptr(m->hq), get_ptr(m->hr),
      get_ptr(m->hidxbx), get_ptr(m->hlbx), get_ptr(m->hubx),
      get_ptr(m->hidxbu), get_ptr(m->hlbu), get_ptr(m->hubu),
      get_ptr(m->hC), get_ptr(m->hD),
      get_ptr(m->hlg), get_ptr(m->hug),
      get_ptr(m->hZl), get_ptr(m->hZu), get_ptr(m->hzl),
      get_ptr(m->hzu), get_ptr(m->hidxs),
      get_ptr(m->hlls), get_ptr(m->hlus), &qp);

    int qp_sol_size = d_ocp_qp_sol_memsize(&dim);
    void *qp_sol_mem = malloc(qp_sol_size);

    struct d_ocp_qp_sol qp_sol;
    d_ocp_qp_sol_create(&dim, &qp_sol, qp_sol_mem);

    int ipm_arg_size = d_ocp_qp_ipm_arg_memsize(&dim);
    void *ipm_arg_mem = malloc(ipm_arg_size);
    struct d_ocp_qp_ipm_arg myarg;
    d_ocp_qp_ipm_arg_create(&dim, &myarg, ipm_arg_mem);

    memcpy(&myarg, &hpipm_options_, sizeof(d_ocp_qp_ipm_arg));

    int ipm_size = d_ocp_qp_ipm_ws_memsize(&dim, &myarg);
    void *ipm_mem = malloc(ipm_size);

    struct d_ocp_qp_ipm_ws workspace;
    d_ocp_qp_ipm_ws_create(&dim, &myarg, &workspace, ipm_mem);

		// solution guess
		for (casadi_int i=0; i<N_+1; i++)	d_ocp_qp_sol_set_u(i, m->hu_guess[i], &qp_sol);
		for (casadi_int i=0; i<N_+1; i++)	d_ocp_qp_sol_set_x(i, m->hx_guess[i], &qp_sol);

    m->fstats.at("solver").tic();
		// call solver
	  d_ocp_qp_ipm_solve(&qp, &qp_sol, &myarg, &workspace);
    d_ocp_qp_ipm_get_status(&workspace, &m->return_status);
    m->fstats.at("solver").toc();
    m->fstats.at("postprocessing").tic();

    d_ocp_qp_ipm_get_iter(&workspace, &m->iter_count);
    d_ocp_qp_ipm_get_max_res_stat(&workspace, &m->res_stat);
	  d_ocp_qp_ipm_get_max_res_eq(&workspace, &m->res_eq);
	  d_ocp_qp_ipm_get_max_res_ineq(&workspace, &m->res_ineq);
	  d_ocp_qp_ipm_get_max_res_comp(&workspace, &m->res_comp);


	  double *stat;
    d_ocp_qp_ipm_get_stat(&workspace, &stat);
printf("\nalpha_aff\tmu_aff\t\tsigma\t\talpha\t\tmu\n");
    d_print_exp_tran_mat(5, m->iter_count, stat, 5);

    printf("\nHPIPM returned with flag %i.\n", m->return_status);
    if(m->return_status == 0)
		{
        printf("\n -> QP solved!\n");
		}
	else if(m->return_status==1)
		{
        printf("\n -> Solver failed! Maximum number of iterations reached\n");
		}
	else if(m->return_status==2)
		{
        printf("\n -> Solver failed! Minimum step lenght reached\n");
		}
	else if(m->return_status==2)
		{
        printf("\n -> Solver failed! NaN in computations\n");
		}
	else
		{
        printf("\n -> Solver failed! Unknown return flag\n");
		}
	printf("\n\n");


	for(casadi_int i=0; i<N_+1; ++i) d_ocp_qp_sol_get_u(i, &qp_sol, m->hu_guess[i]);
	for(casadi_int i=0; i<N_+1; ++i) d_ocp_qp_sol_get_x(i, &qp_sol, m->hx_guess[i]);
	for(casadi_int i=0; i<N_; ++i) d_ocp_qp_sol_get_pi(i, &qp_sol, m->pis[i]);

  if (res[CONIC_LAM_X]) {
    casadi_int offset = 0;
    for(casadi_int i=0; i<N_+1; ++i) {
      std::vector<double> lam_lb(m->nbx[i]+m->nbu[i], nan), lam_ub(m->nbx[i]+m->nbu[i], nan);
      d_ocp_qp_sol_get_lam_lb(i, &qp_sol, get_ptr(lam_lb));
      d_ocp_qp_sol_get_lam_ub(i, &qp_sol, get_ptr(lam_ub));
      for (casadi_int k=0;k<m->nbx[i];++k) {
        res[CONIC_LAM_X][offset+k] = (lam_ub[k+m->nbu[i]]-lam_lb[k+m->nbu[i]])*2;
      }
      offset += m->nbx[i];
      for (casadi_int k=0;k<m->nbu[i];++k) {
        res[CONIC_LAM_X][offset+k] = (lam_ub[k]-lam_lb[k])*2;
      }
      offset += m->nbu[i];
    }
  }
  if (res[CONIC_LAM_A]) {
    casadi_int offset = 0;
    for(casadi_int i=0; i<N_+1; ++i) {
      std::vector<double> lam_lg(m->ng[i]), lam_ug(m->ng[i]);
      d_ocp_qp_sol_get_lam_lg(i, &qp_sol, get_ptr(lam_lg));
      d_ocp_qp_sol_get_lam_ug(i, &qp_sol, get_ptr(lam_ug));
      for (casadi_int k=0;k<m->ng[i];++k) {
        res[CONIC_LAM_A][offset+(i<N_ ? m->nx[i+1]: 0)+k] = (lam_ug[k]-lam_lg[k])*2;
      }
      offset += m->ng[i]+(i<N_ ? m->nx[i+1]: 0);
    }
  }



      uout() << "HPIPM finished after " << m->iter_count << " iterations." << std::endl;
      uout() << "return status: " << m->return_status << std::endl;
      uout() << "HPIPM residuals: " << m->res_stat << ", " << m->res_eq << ", " << m->res_ineq << ", " << m->res_comp << std::endl;

    m->success = m->return_status==0;

    dense_transfer(1.0, get_ptr(m->x), theirs_Xsp_, res[CONIC_X], xsp_, pv);
    dense_transfer(1.0, get_ptr(m->u), theirs_Usp_, res[CONIC_X], usp_, pv);


    // Deal with non-unity I block
    for (casadi_int k=0;k<N_;++k) {
      casadi_int n_row = m->nx[k+1];
      for (casadi_int i=0;i<n_row;++i) {
        double f = -1/m->hI[k][i];
        m->pis[k][i]*=f;
      }
    }

    dense_transfer(2.0, get_ptr(m->pi), pisp_, res[CONIC_LAM_A], lamg_gapsp_, pv);

    // Construct f
    double f = casadi_dot(nx_, arg[CONIC_G], res[CONIC_X]);
    f += 0.5*casadi_bilin(arg[CONIC_H], H_, res[CONIC_X], res[CONIC_X]);

    if (res[CONIC_COST]) res[CONIC_COST][0] = f;

    m->fstats.at("postprocessing").toc();

    return 0;
  }

  Dict HpipmInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<HpipmMemory*>(mem);
    stats["return_status"] = m->return_status;
    stats["iter_count"] = m->iter_count;
    return stats;
  }

  HpipmMemory::HpipmMemory() {
  }

  HpipmMemory::~HpipmMemory() {

  }

  Sparsity HpipmInterface::blocksparsity(casadi_int rows, casadi_int cols,
      const std::vector<Block>& blocks, bool eye) {
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
  void HpipmInterface::blockptr(std::vector<double *>& vs, std::vector<double>& v,
      const std::vector<Block>& blocks, bool eye) {
    casadi_int N = blocks.size();
    vs.resize(N);
    casadi_int offset=0;
    for (casadi_int k=0;k<N;++k) {
      vs[k] = get_ptr(v)+offset;
      if (eye) {
        casadi_assert_dev(blocks[k].rows==blocks[k].cols);
        offset+=blocks[k].rows;
      } else {
        offset+=blocks[k].rows*blocks[k].cols;
      }
    }
  }


  HpipmInterface::HpipmInterface(DeserializingStream& s) : Conic(s) {
    s.version("HpipmInterface", 1);
  }

  void HpipmInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);

    s.version("HpipmInterface", 1);
  }
} // namespace casadi
