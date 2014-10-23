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


#include "dsdp_interface.hpp"

#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/function/mx_function.hpp"
/**
   Some implementation details
   "Multiple cones can be created for the same solver, but it is usually more efficient to group
   all blocks into the same conic structure." user manual
*/
using namespace std;
namespace casadi {

  extern "C"
  int CASADI_SDPSOLVER_DSDP_EXPORT
  casadi_register_sdpsolver_dsdp(SdpSolverInternal::Plugin* plugin) {
    plugin->creator = DsdpInterface::creator;
    plugin->name = "dsdp";
    plugin->doc = DsdpInterface::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_SDPSOLVER_DSDP_EXPORT casadi_load_sdpsolver_dsdp() {
    SdpSolverInternal::registerPlugin(casadi_register_sdpsolver_dsdp);
  }

  DsdpInterface* DsdpInterface::clone() const {
    // Return a deep copy
    DsdpInterface* node = new DsdpInterface(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  DsdpInterface::DsdpInterface(const std::vector<Sparsity> &st) : SdpSolverInternal(st) {
    casadi_assert_message(
      static_cast<double>(m_)*(static_cast<double>(m_)+1)/2 < std::numeric_limits<int>::max(),
      "Your problem size m is too large to be handled by DSDP.");

    addOption("gapTol", OT_REAL, 1e-8,
              "Convergence criterion based on distance between primal and dual objective");
    addOption("maxIter", OT_INTEGER, 500, "Maximum number of iterations");
    addOption("dualTol", OT_REAL, 1e-4,
              "Tolerance for dual infeasibility "
              "(translates to primal infeasibility in dsdp terms)");
    addOption("primalTol", OT_REAL, 1e-4,
              "Tolerance for primal infeasibility "
              "(translates to dual infeasibility in dsdp terms)");
    addOption("stepTol", OT_REAL, 5e-2,
              "Terminate the solver if the step length in the primal is below this tolerance. ");
    addOption("infinity", OT_REAL, 1e30, "Treat numbers higher than this as infinity");
    addOption("_use_penalty", OT_BOOLEAN, true,
              "Modifies the algorithm to use a penality gamma on r.");
    addOption("_penalty", OT_REAL, 1e5,
              "Penality parameter lambda. Must exceed the trace of Y. This "
              "parameter heavily influences the ability of DSDP to treat "
              "linear equalities. The DSDP standard default (1e8) will make "
              "a problem with linear equality return unusable solutions.");
    addOption("_rho", OT_REAL, 4.0, "Potential parameter. Must be >=1");
    addOption("_zbar", OT_REAL, 1e10, "Initial upper bound on the objective of the dual problem.");
    addOption("_reuse", OT_INTEGER, 4,
              "Maximum on the number of times the Schur complement matrix is reused");
    addOption("_printlevel", OT_INTEGER, 1,
              "A printlevel of zero will disable all output. "
              "Another number indicates how often a line is printed.");
    addOption("_loglevel", OT_INTEGER, 0,
              "An integer that specifies how much logging is done on stdout.");

    // Set DSDP memory blocks to null
    dsdp_ = 0;
    sdpcone_ = 0;
  }

  DsdpInterface::~DsdpInterface() {
    if (dsdp_!=0) {
      DSDPDestroy(dsdp_);
      dsdp_ = 0;
    }
  }

  const char* DsdpInterface::terminationReason(int flag) {
    switch (flag) {
    case DSDP_CONVERGED: return "DSDP_CONVERGED";
    case DSDP_MAX_IT: return "DSDP_MAX_IT";
    case DSDP_INFEASIBLE_START: return "DSDP_INFEASIBLE_START";
    case DSDP_INDEFINITE_SCHUR_MATRIX: return "DSDP_INDEFINITE SCHUR";
    case DSDP_SMALL_STEPS: return "DSDP_SMALL_STEPS";
    case DSDP_NUMERICAL_ERROR: return "DSDP_NUMERICAL_ERROR";
    case DSDP_UPPERBOUND: return "DSDP_UPPERBOUND";
    case DSDP_USER_TERMINATION: return "DSDP_USER_TERMINATION";
    case CONTINUE_ITERATING: return "CONTINUE_ITERATING";
    default: return "N/A";
    }
  }

  const char* DsdpInterface::solutionType(int flag) {
    switch (flag) {
    case DSDP_PDFEASIBLE: return  "DSDP_PDFEASIBLE";
    case DSDP_UNBOUNDED: return  "DSDP_UNBOUNDED";
    case DSDP_INFEASIBLE: return  "DSDP_INFEASIBLE";
    case DSDP_PDUNKNOWN: return  "DSDP_PDUNKNOWN";
    default: return "N/A";
    }
  }

  void DsdpInterface::init() {
    // Initialize the base classes
    SdpSolverInternal::init();
    log("DsdpInterface::init", "Enter");

    // Fill the data structures that hold DSDP-style sparse symmetric matrix
    pattern_.resize(n_+1);
    values_.resize(n_+1);

    for (int i=0;i<n_+1;++i) {
      pattern_[i].resize(nb_);
      values_[i].resize(nb_);
      for (int j=0;j<nb_;++j) {
        Sparsity CAij = mapping_.output(i*nb_+j).sparsity();
        pattern_[i][j].resize(CAij.sizeU());
        values_[i][j].resize(pattern_[i][j].size());
        int nz=0;
        const vector<int>& colind = CAij.colind();
        const vector<int>& row = CAij.row();
        for (int cc=0; cc<colind.size()-1; ++cc) {
          int rr;
          // upper triangular part (= lower triangular part for row-major)
          for (int el=colind[cc]; el<colind[cc+1] && (rr=row[el])<=cc; ++el) {
            pattern_[i][j][nz++] = cc*(cc + 1)/2 + rr; // DSDP is row-major --> indices swapped
          }
        }
        mapping_.output(i*nb_+j).get(values_[i][j], SPARSESYM);
      }
    }

    if (nc_>0) {
      // Fill in the linear program structure
      MX A = MX::sym("A", input(SDP_SOLVER_A).sparsity());
      MX LBA = MX::sym("LBA", input(SDP_SOLVER_LBA).sparsity());
      MX UBA = MX::sym("UBA", input(SDP_SOLVER_UBA).sparsity());

      std::vector< MX >  syms;
      syms.push_back(A);
      syms.push_back(LBA);
      syms.push_back(UBA);

      // DSDP has no infinities -- replace by a big number
      // There is a way to deal with this properly, but requires modifying the dsdp source
      // We already did this for the variable bounds
      double dsdp_inf = getOption("infinity");
      MX lba = horzcat(fmax(fmin(LBA, dsdp_inf), -dsdp_inf), A).T();
      MX uba = horzcat(fmax(fmin(UBA, dsdp_inf), -dsdp_inf), A).T();

      mappingA_ = MXFunction(syms, horzcat(-lba, uba).T());
      mappingA_.init();
    }

    if (calc_dual_) {
      store_X_.resize(nb_);
      for (int j=0;j<nb_;++j) {
        store_X_[j].resize(block_sizes_[j]*(block_sizes_[j]+1)/2);
      }
    }
    if (calc_p_) {
      store_P_.resize(nb_);
      for (int j=0;j<nb_;++j) {
        store_P_[j].resize(block_sizes_[j]*(block_sizes_[j]+1)/2);
      }
    }


  }

  void DsdpInterface::evaluate() {
    if (inputs_check_) checkInputs();
    if (print_problem_) printProblem();

    // TODO(Joris): Return flag from DSDP functions not checked!
    // Seems unavoidable to do DSDPCreate here

    // Destroy existing DSDP instance if already allocated
    if (dsdp_!=0) {
      DSDPDestroy(dsdp_);
      dsdp_ = 0;
    }

    // Allocate DSDP solver memory
    DSDPCreate(n_, &dsdp_);
    DSDPSetGapTolerance(dsdp_, getOption("gapTol"));
    DSDPSetMaxIts(dsdp_, getOption("maxIter"));
    DSDPSetPTolerance(dsdp_, getOption("dualTol"));
    DSDPSetRTolerance(dsdp_, getOption("primalTol"));
    DSDPSetStepTolerance(dsdp_, getOption("stepTol"));
    DSDPSetStandardMonitor(dsdp_, getOption("_printlevel"));

    DSDPCreateSDPCone(dsdp_, nb_, &sdpcone_);
    for (int j=0; j<nb_; ++j) {
      log("DsdpInterface::init", "Setting");
      SDPConeSetBlockSize(sdpcone_, j, block_sizes_[j]);
      SDPConeSetSparsity(sdpcone_, j, block_sizes_[j]);
    }
    if (nc_>0) {
      DSDPCreateLPCone(dsdp_, &lpcone_);
      LPConeSetData(lpcone_, nc_*2, getPtr(mappingA_.output().colind()),
                    getPtr(mappingA_.output().row()), mappingA_.output().ptr());
    }

    DSDPCreateBCone(dsdp_, &bcone_);
    BConeAllocateBounds(bcone_, n_);

    DSDPUsePenalty(dsdp_, getOption("_use_penalty") ? 1: 0);
    DSDPSetPenaltyParameter(dsdp_, getOption("_penalty"));
    DSDPSetPotentialParameter(dsdp_, getOption("_rho"));
    DSDPSetZBar(dsdp_, getOption("_zbar"));
    DSDPReuseMatrix(dsdp_, getOption("_reuse") ? 1: 0);
    DSDPLogInfoAllow(getOption("_loglevel"), 0);

    // Copy bounds
    for (int i=0;i<n_;++i) {
      if (input(SDP_SOLVER_LBX).at(i)==-std::numeric_limits< double >::infinity()) {
        BConeSetUnboundedLower(bcone_, i+1);
      } else {
        BConeSetLowerBound(bcone_, i+1, input(SDP_SOLVER_LBX).at(i));
      }
      if (input(SDP_SOLVER_UBX).at(i)==std::numeric_limits< double >::infinity()) {
        BConeSetUnboundedUpper(bcone_, i+1);
      } else {
        BConeSetUpperBound(bcone_, i+1, input(SDP_SOLVER_UBX).at(i));
      }
    }

    if (nc_>0) {
      // Copy linear constraints
      mappingA_.setInput(input(SDP_SOLVER_A), 0);
      mappingA_.setInput(input(SDP_SOLVER_LBA), 1);
      mappingA_.setInput(input(SDP_SOLVER_UBA), 2);
      mappingA_.evaluate();

      // TODO(Joris): this can be made non-allocating bu hacking into DSDP source code
      LPConeSetDataC(lpcone_, nc_*2, mappingA_.output().ptr());

    }

    // Copy b vector
    for (int i=0;i<n_;++i) {
      DSDPSetDualObjective(dsdp_, i+1, -input(SDP_SOLVER_C).at(i));
    }


    if (nb_>0) {
      for (int i=0;i<n_+1;++i) {
        for (int j=0;j<nb_;++j) {
          SDPConeSetASparseVecMat(sdpcone_, j, i, block_sizes_[j], 1, 0, getPtr(pattern_[i][j]),
                                  getPtr(values_[i][j]), pattern_[i][j].size());
        }
      }
    }


    if (nb_>0) {
      // Get Ai from supplied A
      mapping_.setInput(input(SDP_SOLVER_G), 0);
      mapping_.setInput(input(SDP_SOLVER_F), 1);
      mapping_.evaluate();

      for (int i=0;i<n_+1;++i) {
        for (int j=0;j<nb_;++j) {
          mapping_.output(i*nb_+j).get(values_[i][j], SPARSESYM);
        }
      }
    }

    // Manual: Do not set data into the cones after calling this routine.
    // Manual: Should be called only once for each DSDP solver created
    DSDPSetup(dsdp_);
    int info = DSDPSolve(dsdp_);

    casadi_assert_message(info==0, "DsdpInterface failed");

    // Get termination reason
    DSDPTerminationReason reason;
    DSDPStopReason(dsdp_, &reason);
    stats_["termination_reason"] = terminationReason(reason);

    // Get solution type
    DSDPSolutionType pdfeasible;
    DSDPGetSolutionType(dsdp_, &pdfeasible);
    stats_["solution_type"] =  solutionType(pdfeasible);

    info = DSDPGetY(dsdp_, &output(SDP_SOLVER_X).at(0), n_);

    double temp;
    DSDPGetDDObjective(dsdp_, &temp);
    output(SDP_SOLVER_COST).set(-temp);
    DSDPGetPPObjective(dsdp_, &temp);
    output(SDP_SOLVER_DUAL_COST).set(-temp);

    if (calc_dual_) {
      for (int j=0;j<nb_;++j) {
        info = SDPConeComputeX(sdpcone_, j, block_sizes_[j],
                               getPtr(store_X_[j]), store_X_[j].size());
        Pmapper_.input(j).set(store_X_[j], SPARSESYM);
      }
      Pmapper_.evaluate();
      std::copy(Pmapper_.output().data().begin(),
                Pmapper_.output().data().end(),
                output(SDP_SOLVER_DUAL).data().begin());
    }

    if (calc_p_) {
      for (int j=0;j<nb_;++j) {
        info = SDPConeComputeS(sdpcone_, j, 1.0,  output(SDP_SOLVER_X).ptr(), n_, 0,
                               block_sizes_[j] , getPtr(store_P_[j]), store_P_[j].size());
        Pmapper_.input(j).set(store_P_[j], SPARSESYM);
      }
      Pmapper_.evaluate();
      std::copy(Pmapper_.output().data().begin(),
                Pmapper_.output().data().end(),
                output(SDP_SOLVER_P).data().begin());
    }

    DSDPComputeX(dsdp_);

    info = BConeCopyXSingle(bcone_, output(SDP_SOLVER_LAM_X).ptr(), n_);

    if (nc_>0) {
      int dummy;
      double *lam;

      info = LPConeGetXArray(lpcone_, &lam, &dummy);
      std::transform(lam + nc_, lam + 2*nc_, lam,
                      output(SDP_SOLVER_LAM_A).ptr(), std::minus<double>());
    }

  }

} // namespace casadi
