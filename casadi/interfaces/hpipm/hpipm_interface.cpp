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

#include "hpipm_runtime_str.h"

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

    const std::vector<int>& nx = nxs_;
    const std::vector<int>& ng = ngs_;
    const std::vector<int>& nu = nus_;

    Sparsity lamg_csp_, lam_ulsp_, lam_uusp_, lam_xlsp_, lam_xusp_, lam_clsp_;

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

    uout() << "nx,nu,ng" << nx << nu << ng << std::endl;

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

    std::vector< casadi_hpipm_block > theirs_u_blocks, theirs_x_blocks;
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
    std::vector< casadi_hpipm_block > lamg_gap_blocks;
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

    nus_.push_back(0);
    zeros_.resize(N_+1, 0);

    uout() << "nus" << nus_ << std::endl;

    set_hpipm_prob();

    // Allocate memory
    casadi_int sz_arg, sz_res, sz_w, sz_iw;
    casadi_hpipm_work(&p_, &sz_arg, &sz_res, &sz_iw, &sz_w);

    uout() << sz_arg << std::endl;
    uout() << sz_res << std::endl;
    uout() << sz_iw << std::endl;
    uout() << sz_w << std::endl;

    alloc_arg(sz_arg, true);
    alloc_res(sz_res, true);
    alloc_iw(sz_iw, true);
    alloc_w(sz_w, true);



  }

  std::vector<casadi_int> hpipm_blocks_pack(const std::vector<casadi_hpipm_block>& blocks) {
    size_t N = blocks.size();
    std::vector<casadi_int> ret(4*N);
    casadi_int* r = get_ptr(ret);
    for (casadi_int i=0;i<N;++i) {
      *r++ = blocks[i].offset_r;
      *r++ = blocks[i].offset_c;
      *r++ = blocks[i].rows;
      *r++ = blocks[i].cols;
    }
    return ret;
  }


  void codegen_unpack_block(CodeGenerator& g, const std::string& name, const std::vector<casadi_hpipm_block>& blocks) {
      std::string n = "block_" + name + "[" + str(blocks.size()) + "]";
      g.local(n, "static struct casadi_hpipm_block");
      g << "p." << name << " = block_" + name + ";\n";
      g << "casadi_hpipm_unpack_blocks(" << blocks.size()
      << ", p." << name
      << ", " << g.constant(hpipm_blocks_pack(blocks)) << ");\n";
  }

  void codegen_local(CodeGenerator& g, const std::string& name, const std::vector<int>& v) {
    std::string n = name + "[]";
    g.local(n, "static const int");
    std::stringstream init;
    init << "{";
    for (casadi_int i=0;i<v.size();++i) {
      init << v[i];
      if (i<v.size()-1) init << ", ";
    }
    init << "}";
    g.init_local(n, init.str());
  }

  void HpipmInterface::set_hpipm_prob(CodeGenerator& g) const {
    g << "p.qp = &p_qp;\n";
    codegen_local(g, "nx", nxs_);
    codegen_local(g, "nu", nus_);
    codegen_local(g, "ng", ngs_);
    codegen_local(g, "zeros", zeros_);
    g << "p.nx = nx;\n";
    g << "p.nu = nu;\n";
    g << "p.ng = ng;\n";

    g << "p.nbx = nx;\n";
    g << "p.nbu = nu;\n";
    g << "p.ns = zeros;\n";
    g << "p.nsbx = zeros;\n";
    g << "p.nsbu = zeros;\n";
    g << "p.nsg = zeros;\n";

    g << "p.sp_x = " << g.sparsity(sparsity_in(CONIC_X0)) << ";\n";
    g << "p.sp_ba = " << g.sparsity(sparsity_in(CONIC_LBA)) << ";\n";

    g << "p.Asp = " << g.sparsity(Asp_) << ";\n";
    g << "p.Bsp = " << g.sparsity(Bsp_) << ";\n";
    g << "p.Csp = " << g.sparsity(Csp_) << ";\n";
    g << "p.Dsp = " << g.sparsity(Dsp_) << ";\n";

    g << "p.Rsp = " << g.sparsity(Rsp_) << ";\n";
    g << "p.Isp = " << g.sparsity(Isp_) << ";\n";
    g << "p.Ssp = " << g.sparsity(Ssp_) << ";\n";
    g << "p.Qsp = " << g.sparsity(Qsp_) << ";\n";

    g << "p.bsp = " << g.sparsity(bsp_) << ";\n";
    g << "p.xsp = " << g.sparsity(xsp_) << ";\n";
    g << "p.usp = " << g.sparsity(usp_) << ";\n";

    g << "p.pisp = " << g.sparsity(pisp_) << ";\n";

    g << "p.theirs_xsp = " << g.sparsity(theirs_xsp_) << ";\n";
    g << "p.theirs_usp = " << g.sparsity(theirs_usp_) << ";\n";
    g << "p.theirs_Xsp = " << g.sparsity(theirs_Xsp_) << ";\n";
    g << "p.theirs_Usp = " << g.sparsity(theirs_Usp_) << ";\n";

    g << "p.lamg_gapsp = " << g.sparsity(lamg_gapsp_) << ";\n";
    g << "p.lugsp = " << g.sparsity(lugsp_) << ";\n";

    g << "p.N = " << N_ << ";\n";

    g << "p.hpipm_options.mu0 = " << hpipm_options_.mu0 << ";\n";
    g << "p.hpipm_options.alpha_min = " << hpipm_options_.alpha_min << ";\n";
    g << "p.hpipm_options.res_g_max = " << hpipm_options_.res_g_max << ";\n";
    g << "p.hpipm_options.res_b_max = " << hpipm_options_.res_b_max << ";\n";
    g << "p.hpipm_options.res_d_max = " << hpipm_options_.res_d_max << ";\n";
    g << "p.hpipm_options.res_m_max = " << hpipm_options_.res_m_max << ";\n";
    g << "p.hpipm_options.reg_prim = " << hpipm_options_.reg_prim << ";\n";
    g << "p.hpipm_options.lam_min = " << hpipm_options_.lam_min << ";\n";
    g << "p.hpipm_options.t_min = " << hpipm_options_.t_min << ";\n";
    g << "p.hpipm_options.tau_min = " << hpipm_options_.tau_min << ";\n";
    g << "p.hpipm_options.iter_max = " << hpipm_options_.iter_max << ";\n";
    g << "p.hpipm_options.stat_max = " << hpipm_options_.stat_max << ";\n";
    g << "p.hpipm_options.pred_corr = " << hpipm_options_.pred_corr << ";\n";
    g << "p.hpipm_options.cond_pred_corr = " << hpipm_options_.cond_pred_corr << ";\n";
    g << "p.hpipm_options.itref_pred_max = " << hpipm_options_.itref_pred_max << ";\n";
    g << "p.hpipm_options.itref_corr_max = " << hpipm_options_.itref_corr_max << ";\n";
    g << "p.hpipm_options.warm_start = " << hpipm_options_.warm_start << ";\n";
    g << "p.hpipm_options.square_root_alg = " << hpipm_options_.square_root_alg << ";\n";
    g << "p.hpipm_options.lq_fact = " << hpipm_options_.lq_fact << ";\n";
    g << "p.hpipm_options.abs_form = " << hpipm_options_.abs_form << ";\n";
    g << "p.hpipm_options.comp_dual_sol_eq = " << hpipm_options_.comp_dual_sol_eq << ";\n";
    g << "p.hpipm_options.comp_res_exit = " << hpipm_options_.comp_res_exit << ";\n";
    g << "p.hpipm_options.comp_res_pred = " << hpipm_options_.comp_res_pred << ";\n";
    g << "p.hpipm_options.split_step = " << hpipm_options_.split_step << ";\n";
    g << "p.hpipm_options.var_init_scheme = " << hpipm_options_.var_init_scheme << ";\n";
    g << "p.hpipm_options.t_lam_min = " << hpipm_options_.t_lam_min << ";\n";
    g << "p.hpipm_options.mode = " << hpipm_options_.mode << ";\n";
    g << "p.hpipm_options.memsize = " << hpipm_options_.memsize << ";\n";


    g << "p.inf = " << inf_ << ";\n";

    codegen_unpack_block(g, "A", A_blocks);
    codegen_unpack_block(g, "B", B_blocks);
    codegen_unpack_block(g, "C", C_blocks);
    codegen_unpack_block(g, "D", D_blocks);

    codegen_unpack_block(g, "R", R_blocks);
    codegen_unpack_block(g, "I", I_blocks);
    codegen_unpack_block(g, "S", S_blocks);
    codegen_unpack_block(g, "Q", Q_blocks);  
  
    codegen_unpack_block(g, "b", b_blocks);
    codegen_unpack_block(g, "lug", lug_blocks);
    codegen_unpack_block(g, "u", u_blocks);
    codegen_unpack_block(g, "x", x_blocks);   


    codegen_unpack_block(g, "lam_ul", lam_ul_blocks);
    codegen_unpack_block(g, "lam_xl", lam_xl_blocks);
    codegen_unpack_block(g, "lam_uu", lam_uu_blocks);
    codegen_unpack_block(g, "lam_xu", lam_xu_blocks);
    codegen_unpack_block(g, "lam_cl", lam_cl_blocks);
    codegen_unpack_block(g, "lam_cu", lam_cu_blocks);

    g << "casadi_hpipm_setup(&p);\n";

  }

  void HpipmInterface::set_hpipm_prob() {
    p_.qp = &p_qp_;
    p_.nx  = get_ptr(nxs_);
    p_.nu  = get_ptr(nus_);
    uout() << "nu" << nus_.size() << ":" << N_ << std::endl;
    p_.ng  = get_ptr(ngs_);

    p_.nbx = get_ptr(nxs_);
    p_.nbu = get_ptr(nus_);
    p_.ns  = get_ptr(zeros_);
    p_.nsbx = get_ptr(zeros_);
    p_.nsbu = get_ptr(zeros_);
    p_.nsg = get_ptr(zeros_);

    p_.sp_x = sparsity_in(CONIC_X0);
    p_.sp_ba = sparsity_in(CONIC_LBA);

    p_.Asp = Asp_;
    p_.Bsp = Bsp_;
    p_.Csp = Csp_;
    p_.Dsp = Dsp_;

    p_.Rsp = Rsp_;
    p_.Isp = Isp_;
    p_.Ssp = Ssp_;
    p_.Qsp = Qsp_;

    p_.bsp = bsp_;

    p_.xsp = xsp_;
    p_.usp = usp_;

    p_.pisp = pisp_;

    p_.theirs_xsp = theirs_xsp_;
    p_.theirs_usp = theirs_usp_;
    p_.theirs_Xsp = theirs_Xsp_;
    p_.theirs_Usp = theirs_Usp_;

    p_.lamg_gapsp = lamg_gapsp_;
    p_.lugsp = lugsp_;

    p_.N = N_;
    p_.hpipm_options = hpipm_options_;

    p_.inf = inf_;

    p_.A = get_ptr(A_blocks);
    p_.B = get_ptr(B_blocks);
    p_.C = get_ptr(C_blocks);
    p_.D = get_ptr(D_blocks);

    p_.R = get_ptr(R_blocks);
    p_.I = get_ptr(I_blocks);
    p_.S = get_ptr(S_blocks);
    p_.Q = get_ptr(Q_blocks);

    p_.b = get_ptr(b_blocks);
    p_.lug = get_ptr(lug_blocks);
    p_.u = get_ptr(u_blocks);
    p_.x = get_ptr(x_blocks);

    p_.lam_ul = get_ptr(lam_ul_blocks);
    p_.lam_xl = get_ptr(lam_xl_blocks);
    p_.lam_uu = get_ptr(lam_uu_blocks);
    p_.lam_xu = get_ptr(lam_xu_blocks);
    p_.lam_cl = get_ptr(lam_cl_blocks);
    p_.lam_cu = get_ptr(lam_cu_blocks);

    uout() << "p"<< p_.nu[0] << nus_[0] << std::endl;
    casadi_hpipm_setup(&p_);

    uout() << "p"<< p_.nu[0] << std::endl;
  }

  int HpipmInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<HpipmMemory*>(mem);

    m->fstats["preprocessing"]  = FStats();
    m->fstats["solver"]         = FStats();
    m->fstats["postprocessing"] = FStats();
    return 0;
  }

  /** \brief Set the (persistent) work vectors */
  void HpipmInterface::set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const {

    auto m = static_cast<HpipmMemory*>(mem);

    uout() << "set_work" << std::endl;

    uout() << "p"<< p_.nu[0] << nus_[0] << std::endl;
    Conic::set_work(mem, arg, res, iw, w);

    m->d.prob = &p_;
    m->d.qp = &m->d_qp;
    uout() << "p_" << &p_ << p_.nu[0] << std::endl;
    casadi_hpipm_init(&m->d, &arg, &res, &iw, &w);

    m->iter_count = 0;

    uout() << "pset_work"<< p_.nu[0] << std::endl;
  }

  int HpipmInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<HpipmMemory*>(mem);
    // Statistics
    m->fstats.at("preprocessing").tic();

    casadi_hpipm_solve(&m->d, arg, res, iw, w);

    m->success = m->d.return_status==0;

    uout() << "HPIPM finished after " << m->d.iter_count << " iterations." << std::endl;
    uout() << "return status: " << m->d.return_status << std::endl;
    uout() << "HPIPM residuals: " << m->d.res_stat << ", " << m->d.res_eq << ", " << m->d.res_ineq << ", " << m->d.res_comp << std::endl;

    return 0;
  }

  Dict HpipmInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<HpipmMemory*>(mem);

    stats["return_status"] = m->d.return_status;
    stats["iter_count"] = m->d.iter_count;
    return stats;
  }

  HpipmMemory::HpipmMemory() {
  }

  HpipmMemory::~HpipmMemory() {
  }

  Sparsity HpipmInterface::blocksparsity(casadi_int rows, casadi_int cols,
      const std::vector<casadi_hpipm_block>& blocks, bool eye) {
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
      const std::vector<casadi_hpipm_block>& blocks, bool eye) {
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


  void HpipmInterface::codegen_body(CodeGenerator& g) const {
    qp_codegen_body(g);
    g.add_auxiliary(CodeGenerator::AUX_PROJECT);
    g.add_auxiliary(CodeGenerator::AUX_SCAL);
    g.add_auxiliary(CodeGenerator::AUX_SPARSIFY);
    g.add_auxiliary(CodeGenerator::AUX_MAX);
    g.add_auxiliary(CodeGenerator::AUX_SPARSITY);
    g.add_auxiliary(CodeGenerator::AUX_SUM);
    g.add_auxiliary(CodeGenerator::AUX_FILL);
    g.add_auxiliary(CodeGenerator::AUX_CLIP_MIN);
    g.add_auxiliary(CodeGenerator::AUX_CLIP_MAX);
    g.add_auxiliary(CodeGenerator::AUX_DOT);
    g.add_auxiliary(CodeGenerator::AUX_BILIN);
    g.add_include("blasfeo_d_aux_ext_dep.h");
    g.add_include("hpipm_d_ocp_qp_ipm.h");
    g.add_include("hpipm_d_ocp_qp_dim.h");
    g.add_include("hpipm_d_ocp_qp.h");
    g.add_include("hpipm_d_ocp_qp_sol.h");
    g.add_include("hpipm_d_ocp_qp_utils.h");
    g.add_include("stdlib.h");
    g.add_include("string.h");

    g.auxiliaries << g.sanitize_source(hpipm_runtime_str, {"casadi_real"});


    g.local("d", "struct casadi_hpipm_data");
    g.local("p", "struct casadi_hpipm_prob");

    set_hpipm_prob(g);

    // Setup data structure (corresponds to set_work)
    g << "d.prob = &p;\n";
    g << "d.qp = &d_qp;\n";
    g << "casadi_hpipm_init(&d, &arg, &res, &iw, &w);\n";

    g << "casadi_hpipm_solve(&d, arg, res, iw, w);\n";

    g << "if (d.return_status!=0) {\n";
    if (error_on_fail_) {
      g << "return -1000;\n";
    } else {
      g << "return -1;\n";
    }
    g << "}\n";
    g << "return 0;\n";
  }

  HpipmInterface::HpipmInterface(DeserializingStream& s) : Conic(s) {
    s.version("HpipmInterface", 1);
  }

  void HpipmInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);

    s.version("HpipmInterface", 1);
  }

} // namespace casadi
