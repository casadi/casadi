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


#include "smith_lr_dle_internal.hpp"
#include <cassert>
#include "../core/std_vector_tools.hpp"
#include "../core/matrix/matrix_tools.hpp"
#include "../core/mx/mx_tools.hpp"
#include "../core/sx/sx_tools.hpp"
#include "../core/function/mx_function.hpp"
#include "../core/function/sx_function.hpp"
#include <iomanip>

#include <numeric>

INPUTSCHEME(LR_DLEInput)

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LRDLESOLVER_SMITH_EXPORT
  casadi_register_lrdlesolver_smith(LrDleInternal::Plugin* plugin) {
    plugin->creator = SmithLrDleInternal::creator;
    plugin->name = "smith";
    plugin->doc = SmithLrDleInternal::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_LRDLESOLVER_SMITH_EXPORT casadi_load_lrdlesolver_smith() {
    LrDleInternal::registerPlugin(casadi_register_lrdlesolver_smith);
  }

  SmithLrDleInternal::SmithLrDleInternal(
      const LrDleStructure& st, const std::vector<int> &Hs, int nfwd, int nadj) :
      LrDleInternal(st, Hs), nfwd_(nfwd), nadj_(nadj) {

    // set default options
    setOption("name", "unnamed_smith_lr_dle_solver"); // name of the function

    addOption("max_iter", OT_INTEGER, 100,   "Maximum number of iterations for the algorithm");
    addOption("tol",      OT_REAL,    1e-12, "Tolerance for satisfying the Lyapunov equation.");

    addOption("print_iteration", OT_BOOLEAN, false,
              "Print information about each iteration");

  }

  SmithLrDleInternal::~SmithLrDleInternal() {

  }

  void SmithLrDleInternal::smith_prepare() {

    n_ = A_.size1();
    e_ = DMatrix::eye(n_);

    // Allocate helper data-structures
    DMatrix C = DMatrix::zeros(C_);
    DMatrix A = DMatrix::zeros(A_);
    DMatrix V = DMatrix::zeros(V_);
    DMatrix H = DMatrix::zeros(H_);

    D_.clear();    D_.reserve(max_iter_+1);
    VD_.clear();   VD_.reserve(max_iter_+1);
    DTH_.clear();  DTH_.reserve(max_iter_+1);
    VDTH_.clear(); VDTH_.reserve(max_iter_+1);

    D_.push_back(with_C_? C : e_);
    VD_.push_back(with_C_? mul(V, C.T()): V);
    DTH_.push_back(mul(with_C_? C.T() : e_, with_H_?H: e_));
    VDTH_.push_back(mul(V, DTH_[0]));
    DTHv_.push_back(with_H_? horzsplit(DTH_[0], Hi_) : std::vector<DMatrix>(1, DTH_[0]));
    VDTHv_.push_back(with_H_? horzsplit(VDTH_[0], Hi_) : std::vector<DMatrix>(1, VDTH_[0]));

    for (int k=0;k<max_iter_;++k) {
      D_.push_back(mul(A, D_[D_.size()-1]));
      VD_.push_back(mul(V, D_[D_.size()-1].T()));
      DTH_.push_back(mul(D_[D_.size()-1].T(), with_H_?H: e_));
      VDTH_.push_back(mul(V, VDTH_[VDTH_.size()-1]));
      if (with_H_) {
        DTHv_.push_back(horzsplit(DTH_[D_.size()-1], Hi_));
        VDTHv_.push_back(horzsplit(VDTH_[D_.size()-1], Hi_));
      }
    }

    if (with_H_) {
      // Build up indices to effectively slice H_j from H
      for (int j=0;j<max_iter_;++j) {
        DTHvi_.push_back(std::vector<int>(Hs_.size()+1, 0));
        VDTHvi_.push_back(std::vector<int>(Hs_.size()+1, 0));
        for (int i=0;i<Hs_.size();++i) {
          DTHvi_[j][i+1]  = DTHvi_[j][i] + DTHv_[j][i].size();
          VDTHvi_[j][i+1] = VDTHvi_[j][i] + VDTHv_[j][i].size();
        }
      }
    }

    /// Variable to hold the symmetrized DLE_V input
    Vsymm_  = input(LR_DLE_V);

    if (nfwd_>0 || nadj_>0) {
      // Tape for forward mode

      D_ad_.clear();    D_ad_.reserve(max_iter_+1);
      VD_ad_.clear();   VD_ad_.reserve(max_iter_+1);
      DTH_ad_.clear();  DTH_ad_.reserve(max_iter_+1);
      VDTH_ad_.clear(); VDTH_ad_.reserve(max_iter_+1);

      DTH_ad_.insert(DTH_ad_.begin(), DTH_.begin(), DTH_.end());
      VDTH_ad_.insert(VDTH_ad_.begin(), VDTH_.begin(), VDTH_.end());
      DTHv_ad_.insert(DTHv_ad_.begin(), DTHv_.begin(), DTHv_.end());
      VDTHv_ad_.insert(VDTHv_ad_.begin(), VDTHv_.begin(), VDTHv_.end());
      D_ad_.insert(D_ad_.begin(), D_.begin(), D_.end());
      VD_ad_.insert(VD_ad_.begin(), VD_.begin(), VD_.end());
    }

    if (nfwd_>0) {
      Vdsymm_ = input(LR_DLE_V);
    }



    // Allocate work vectors for norm_inf_mul_tt
    Dwork_norm_.resize(n_);
    Iwork_norm_.resize(n_+1+n_);


    em_ = DMatrix::eye(V_.size1())*0.5;

  }

  void SmithLrDleInternal::init() {
    max_iter_  = getOption("max_iter");
    tol_  = getOption("tol");

    LrDleInternal::init();

    // Adapt the inputs/outputs if AD is requested
    if (nfwd_>0 || nadj_>0) {
      setNumInputs(LR_DLE_NUM_IN*(nfwd_+1)+LR_DLE_NUM_OUT*nadj_);
      setNumOutputs(LR_DLE_NUM_IN*nadj_+LR_DLE_NUM_OUT*(nfwd_+1));

      for (int i=0;i<nfwd_;++i) {
        for (int j=0;j<LR_DLE_NUM_IN;++j) {
          input(LR_DLE_NUM_IN*(i+1)+j) = input(j);
        }
        for (int j=0;j<LR_DLE_NUM_OUT;++j) {
          output(LR_DLE_NUM_OUT*(i+1)+j) = output(j);
        }
      }

      for (int i=0;i<nadj_;++i) {
        for (int j=0;j<LR_DLE_NUM_IN;++j) {
          output(LR_DLE_NUM_OUT*(nfwd_+1)+LR_DLE_NUM_IN*i+j) = input(j);
        }
        for (int j=0;j<LR_DLE_NUM_OUT;++j) {
          input(LR_DLE_NUM_IN*(nfwd_+1)+LR_DLE_NUM_OUT*i+j) = output(j);
        }
      }
      input_.scheme = IOScheme();
    }

    /// Prepare intermediate variables and AD tape
    smith_prepare();

    casadi_assert_message(!pos_def_,
      "pos_def option set to True: Solver only handles the indefinite case.");

    print_iteration_ = getOption("print_iteration");
  }

  void SmithLrDleInternal::printIteration(std::ostream &stream) {
    stream << setw(5) << "iter";
    stream << setw(10) << "norm_inf";
    stream << std::endl;
    stream.unsetf(std::ios::floatfield);
  }

  void SmithLrDleInternal::printIteration(std::ostream &stream, int iter, double norm_inf) {
    stream << setw(5) << iter;
    stream << setw(10) << scientific << setprecision(2) << norm_inf;

    stream << fixed;
    stream << std::endl;
    stream.unsetf(std::ios::floatfield);
  }

  void SmithLrDleInternal::evaluate() {

    // Set-up aliases for inputs
    DMatrix &A = input(LR_DLE_A);
    DMatrix &C = input(LR_DLE_C);
    DMatrix &V = input(LR_DLE_V);
    DMatrix &H = input(LR_DLE_H);

    // Symmetrize V
    Vsymm_.setAll(0);
    DMatrix::mul_no_alloc_nn(V, em_, Vsymm_);
    DMatrix::mul_no_alloc_tn(V, em_, Vsymm_);

    // Initializing D
    D_[0].set(with_C_? C : e_);

    // Computing VC^T
    VD_[0].setAll(0); // Clear variable
    DMatrix::mul_no_alloc_nt(Vsymm_, with_C_? C : e_, VD_[0]);

    // Iteraton counter
    int i = 0;
    while (true) {

      // Clear intermediate variables
      D_[i+1].setAll(0);
      VD_[i+1].setAll(0);

      // The heart of the algorithm; obtaining a new D-entry
      DMatrix::mul_no_alloc_nn(A, D_[i], D_[i+1]);

      // Calculate the stopping criterion
      DMatrix::mul_no_alloc_nt(Vsymm_, D_[i+1], VD_[i+1]);
      double norm = norm_inf_mul_nn(D_[i+1], VD_[i+1], Dwork_norm_, Iwork_norm_);

      if (print_iteration_) {
        // Only print iteration header oncein a while
        if (i % 10==0) {
          printIteration(std::cout);
        }

        // Print iteration information
        printIteration(std::cout, i, norm);
      }

      // Break if stopping criterion reached
      if (norm < tol_) {
        break;
      }

      // Abort after too much iterations
      if (i>=max_iter_-1) {
        casadi_error("Maximum number of iterations reached");
      }

      i+=1;
    }

    // The final number of iterations
    int r = i;


    if (with_H_) {
      // Reset Outputs
      for (int k=0;k<Hs_.size();++k) {
        Pv_[k].setAll(0);
      }
    } else {
      output(LR_DLE_Y).setAll(0);
    }

    // Calculating the outputs
    for (int j=0;j<r;++j) {

      // Clear intermediate variables
      DTH_[j].setAll(0);
      VDTH_[j].setAll(0);

      // Obtain H_j^T P H_j without ever calculating P:
      // Instead work with products of D and H_j
      if (with_H_) {
        DMatrix::mul_no_alloc_tn(D_[j], H, DTH_[j]);
      } else {
        DMatrix::mul_no_alloc_tn(D_[j], e_, DTH_[j]);
      }

      DMatrix::mul_no_alloc_nn(V, DTH_[j], VDTH_[j]);

      if (with_H_) {
        for (int k=0;k<Hs_.size();++k) {
          // We need the equivalent of a slice to obtain H_j from H
          std::copy(DTH_[j].begin()+DTHvi_[j][k],
                    DTH_[j].begin()+DTHvi_[j][k+1],
                    DTHv_[j][k].begin());
          std::copy(VDTH_[j].begin()+VDTHvi_[j][k],
                    VDTH_[j].begin()+VDTHvi_[j][k+1],
                    VDTHv_[j][k].begin());

          DMatrix::mul_no_alloc_tn(DTHv_[j][k], VDTHv_[j][k], Pv_[k]);
          std::copy(Pv_[k].begin(),
                    Pv_[k].end(),
                    output(LR_DLE_Y).begin()+Pi_[k]);
        }
      } else {
        DMatrix::mul_no_alloc_tn(DTH_[j], VDTH_[j], output(LR_DLE_Y));
      }
    }

    // Forward sweeps
    for (int d=0;d<nfwd_;++d) {

      // Set-up aliases for forward seeds
      DMatrix & Vd = input(LR_DLE_NUM_IN*(d+1)+LR_DLE_V);
      DMatrix & Cd = input(LR_DLE_NUM_IN*(d+1)+LR_DLE_C);
      DMatrix & Ad = input(LR_DLE_NUM_IN*(d+1)+LR_DLE_A);
      DMatrix & Hd = input(LR_DLE_NUM_IN*(d+1)+LR_DLE_H);

      // Symmetrize Vd
      Vdsymm_.setAll(0);
      DMatrix::mul_no_alloc_nn(Vd, em_, Vdsymm_);
      DMatrix::mul_no_alloc_tn(Vd, em_, Vdsymm_);

      if (with_C_) {
        std::copy(Cd.begin(), Cd.end(), D_ad_[0].begin());
      } else {
        std::fill(D_ad_[0].begin(), D_ad_[0].end(), 0);
        //std::copy(e_.begin(), e_.end(), D_ad_[0].begin());
      }

      VD_ad_[0].setAll(0); // Clear variable

      DMatrix::mul_no_alloc_nt(Vdsymm_, with_C_ ? C : e_, VD_ad_[0]);
      if (with_C_) {
        DMatrix::mul_no_alloc_nt(V, Cd, VD_ad_[0]);
      }

      // Replay main-loop
      for (int i=0;i<r;++i) {
       D_ad_[i+1].setAll(0);
       DMatrix::mul_no_alloc_nn(A, D_ad_[i], D_ad_[i+1]);
       DMatrix::mul_no_alloc_nn(Ad, D_[i], D_ad_[i+1]);
      }

      if (with_H_) {
        for (int k=0;k<Hs_.size();++k) {
          Pv_[k].setAll(0);
        }
      } else {
        output(LR_DLE_NUM_OUT*(d+1)+LR_DLE_Y).setAll(0);
      }

      // Calculating the outputs
      for (int j=0;j<r;++j) {

        // Clear intermediate variables
        DTH_ad_[j].setAll(0);
        VDTH_ad_[j].setAll(0);

        if (with_H_) {
          DMatrix::mul_no_alloc_tn(D_ad_[j], H, DTH_ad_[j]);
          DMatrix::mul_no_alloc_tn(D_[j], Hd, DTH_ad_[j]);
        } else {
          DMatrix::mul_no_alloc_tn(D_ad_[j], e_, DTH_ad_[j]);
        }
        DMatrix::mul_no_alloc_nn(Vdsymm_, DTH_[j], VDTH_ad_[j]);
        DMatrix::mul_no_alloc_nn(V, DTH_ad_[j], VDTH_ad_[j]);

        if (with_H_) {
          for (int k=0;k<Hs_.size();++k) {
            std::copy(DTH_ad_[j].begin()+DTHvi_[j][k],
              DTH_ad_[j].begin()+DTHvi_[j][k+1], DTHv_ad_[j][k].begin());
            std::copy(VDTH_ad_[j].begin()+VDTHvi_[j][k],
              VDTH_ad_[j].begin()+VDTHvi_[j][k+1], VDTHv_ad_[j][k].begin());
            DMatrix::mul_no_alloc_tn(DTHv_ad_[j][k], VDTHv_[j][k], Pv_[k]);
            DMatrix::mul_no_alloc_tn(DTHv_[j][k], VDTHv_ad_[j][k], Pv_[k]);

          }
        } else {
          DMatrix::mul_no_alloc_tn(DTH_ad_[j], VDTH_[j], output(LR_DLE_NUM_OUT*(d+1)+LR_DLE_Y));
          DMatrix::mul_no_alloc_tn(DTH_[j], VDTH_ad_[j], output(LR_DLE_NUM_OUT*(d+1)+LR_DLE_Y));
        }
      }

      if (with_H_) {
        for (int k=0;k<Hs_.size();++k) {
          std::copy(Pv_[k].begin(),
          Pv_[k].end(),
          output(LR_DLE_NUM_OUT*(d+1)+LR_DLE_Y).begin()+Pi_[k]);
        }
      }
    }

    // Adjoint sweeps
    for (int d=0;d<nadj_;++d) {

      /// Reset bar quantities
      for (int j=0;j<r;++j) {
        std::fill(DTH_ad_[j].begin(), DTH_ad_[j].end(), 0);
        std::fill(VDTH_ad_[j].begin(), VDTH_ad_[j].end(), 0);
        std::fill(D_ad_[j].begin(), D_ad_[j].end(), 0);
        std::fill(VD_ad_[j].begin(), VD_ad_[j].end(), 0);

        if (with_H_) {
          for (int i=0;i<Hs_.size();++i) {
            std::fill(DTHv_ad_[j][i].begin(), DTHv_ad_[j][i].end(), 0);
            std::fill(VDTHv_ad_[j][i].begin(), VDTHv_ad_[j][i].end(), 0);
          }
        }
      }

      DMatrix & Vb = output(LR_DLE_NUM_OUT*(nfwd_+1)+LR_DLE_NUM_IN*d+LR_DLE_V);
      DMatrix & Cb = output(LR_DLE_NUM_OUT*(nfwd_+1)+LR_DLE_NUM_IN*d+LR_DLE_C);
      DMatrix & Ab = output(LR_DLE_NUM_OUT*(nfwd_+1)+LR_DLE_NUM_IN*d+LR_DLE_A);
      DMatrix & Hb = output(LR_DLE_NUM_OUT*(nfwd_+1)+LR_DLE_NUM_IN*d+LR_DLE_H);

      Vb.setAll(0);
      Cb.setAll(0);
      Ab.setAll(0);
      Hb.setAll(0);

      if (with_H_) {
        for (int k=0;k<Hs_.size();++k) {
           std::copy(input(LR_DLE_NUM_IN*(nfwd_+1)+LR_DLE_NUM_OUT*d+LR_DLE_Y).begin()+Pi_[k],
                     input(LR_DLE_NUM_IN*(nfwd_+1)+LR_DLE_NUM_OUT*d+LR_DLE_Y).begin()+Pi_[k+1],
                     Pv_[k].begin());
        }
      }

      for (int j=0;j<r;++j) {

        if (with_H_) {
          for (int k=0;k<Hs_.size();++k) {
            DMatrix & Yb = Pv_[k];

            DMatrix::mul_no_alloc_nt(VDTHv_[j][k], Yb, DTHv_ad_[j][k]);
            DMatrix::mul_no_alloc_nt(DTHv_[j][k], Yb, VDTHv_ad_[j][k]);

            DMatrix::mul_no_alloc_nn(VDTHv_[j][k], Yb, DTHv_ad_[j][k]);
            DMatrix::mul_no_alloc_nn(DTHv_[j][k], Yb, VDTHv_ad_[j][k]);

            for (int kk=0;kk<DTHv_ad_[j][k].size();++kk) DTHv_ad_[j][k][kk]*=0.5;
            for (int kk=0;kk<VDTHv_ad_[j][k].size();++kk) VDTHv_ad_[j][k][kk]*=0.5;

            std::copy(DTHv_ad_[j][k].begin(),
                     DTHv_ad_[j][k].end(),
                     DTH_ad_[j].begin()+DTHvi_[j][k]);
            std::copy(VDTHv_ad_[j][k].begin(),
                     VDTHv_ad_[j][k].end(),
                     VDTH_ad_[j].begin()+DTHvi_[j][k]);
          }
        } else {
          DMatrix & Yb = input(LR_DLE_NUM_IN*(nfwd_+1)+LR_DLE_NUM_OUT*d+LR_DLE_Y);
          DMatrix::mul_no_alloc_nt(VDTH_[j], Yb, DTH_ad_[j]);
          DMatrix::mul_no_alloc_nt(DTH_[j], Yb, VDTH_ad_[j]);

          DMatrix::mul_no_alloc_nn(VDTH_[j], Yb, DTH_ad_[j]);
          DMatrix::mul_no_alloc_nn(DTH_[j], Yb, VDTH_ad_[j]);

          for (int kk=0;kk<DTH_ad_[j].size();++kk) DTH_ad_[j][kk]*=0.5;
          for (int kk=0;kk<VDTH_ad_[j].size();++kk) VDTH_ad_[j][kk]*=0.5;

        }
        DMatrix::mul_no_alloc_nt(VDTH_ad_[j], DTH_[j], Vb);
        DMatrix::mul_no_alloc_tn(V, VDTH_ad_[j], DTH_ad_[j]);

        if (with_H_)  {
          DMatrix::mul_no_alloc_nn(D_[j], DTH_ad_[j], Hb);
        }
        DMatrix::mul_no_alloc_nt(with_H_? H : e_, DTH_ad_[j], D_ad_[j]);

      }
      // Replay main-loop in reverse order
      for (int i=r;i>=0;--i) {
       DMatrix::mul_no_alloc_nt(D_ad_[i+1], D_[i], Ab);
       DMatrix::mul_no_alloc_tn(A, D_ad_[i+1], D_ad_[i]);
      }

      if (with_C_) {
        DMatrix::mul_no_alloc_tn(VD_ad_[0], V, Cb);
      }
      DMatrix::mul_no_alloc_nn(VD_ad_[0], with_C_? C : e_, Vb);

      if (with_C_) {
        DMatrix::mul_no_alloc_nn(e_, D_ad_[0], Cb);
      }

      // Clear adjoint seeds
      input(LR_DLE_NUM_IN*(nfwd_+1)+LR_DLE_NUM_OUT*d+LR_DLE_Y).setAll(0);
    }

    if (gather_stats_) {
      stats_["iter_count"] = r;
    }
  }

  Function SmithLrDleInternal::getDerivative(int nfwd, int nadj) {
    casadi_assert_message((nfwd_==0 && nadj_==0) || (nfwd==0 && nadj==0),
      "SmithLrDleInternal::second order derivatives are not supported");
    // Return a deep copy
    SmithLrDleInternal* node = new SmithLrDleInternal(st_, Hs_, nfwd, nadj);
    node->setOption(dictionary());
    node->init();
    return node->shared_from_this<Function>();
  }

  void SmithLrDleInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    LrDleInternal::deepCopyMembers(already_copied);
  }

  SmithLrDleInternal* SmithLrDleInternal::clone() const {
    // Return a deep copy
    SmithLrDleInternal* node = new SmithLrDleInternal(st_, Hs_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


