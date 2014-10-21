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


#include "ooqp_interface.hpp"
#include "casadi/core/function/qp_solver.hpp"
#include "casadi/core/matrix/sparsity_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/matrix/generic_matrix_tools.hpp"
#include "casadi/core/std_vector_tools.hpp"

// OOQP headers
#include <cQpGenSparse.h>
#include <Status.h>
#include <GondzioSolver.h>

// A variable that controls the printlevel of OOQP
// This is the only possible way to access it using the C++ interface
extern int gOoqpPrintLevel;

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_QPSOLVER_OOQP_EXPORT
  casadi_register_qpsolver_ooqp(QpSolverInternal::Plugin* plugin) {
    plugin->creator = OoqpInterface::creator;
    plugin->name = "ooqp";
    plugin->doc = OoqpInterface::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_QPSOLVER_OOQP_EXPORT casadi_load_qpsolver_ooqp() {
    QpSolverInternal::registerPlugin(casadi_register_qpsolver_ooqp);
  }

  OoqpInterface::OoqpInterface(const std::vector<Sparsity>& st) : QpSolverInternal(st) {
    addOption("print_level", OT_INTEGER, 0,
              "Print level. OOQP listens to print_level 0, 10 and 100");
    addOption("mutol", OT_REAL, 1e-8, "tolerance as provided with setMuTol to OOQP");
    addOption("artol", OT_REAL, 1e-8, "tolerance as provided with setArTol to OOQP");
  }

  OoqpInterface::~OoqpInterface() {
  }

  void OoqpInterface::init() {
    // Initialize the base classes
    QpSolverInternal::init();

    // Read options
    print_level_ = getOption("print_level");
    mutol_ = getOption("mutol");
    artol_ = getOption("artol");

    // Allocate memory for problem
    c_.resize(n_);
    xlow_.resize(n_);
    ixlow_.resize(n_);
    xupp_.resize(n_);
    ixupp_.resize(n_);
    bA_.resize(nc_);
    clow_.resize(nc_);
    iclow_.resize(nc_);
    cupp_.resize(nc_);
    icupp_.resize(nc_);
    dQ_.resize(input(QP_SOLVER_H).sizeU());
    irowQ_.resize(dQ_.size());
    jcolQ_.resize(dQ_.size());
    dA_.resize(input(QP_SOLVER_A).size());
    irowA_.resize(dA_.size());
    jcolA_.resize(dA_.size());
    dC_.resize(dA_.size());
    irowC_.resize(dC_.size());
    jcolC_.resize(dC_.size());
    x_index_.resize(n_);
    c_index_.resize(nc_);
    p_.resize(n_);
    AT_ = DMatrix::zeros(input(QP_SOLVER_A).sparsity().transpose());
    AT_tmp_.resize(nc_);

    // Solution
    x_.resize(n_);
    gamma_.resize(n_);
    phi_.resize(n_);
    y_.resize(nc_);
    z_.resize(nc_);
    lambda_.resize(nc_);
    pi_.resize(nc_);
  }

  void OoqpInterface::evaluate() {
    // Check inputs for consistency
    if (inputs_check_) checkInputs();

    // Get problem data
    const vector<double>& lbx = input(QP_SOLVER_LBX).data();
    const vector<double>& ubx = input(QP_SOLVER_UBX).data();
    const vector<double>& lba = input(QP_SOLVER_LBA).data();
    const vector<double>& uba = input(QP_SOLVER_UBA).data();
    const vector<double>& g = input(QP_SOLVER_G).data();

    // Parameter contribution to the objective
    double objParam = 0;

    // Get the number of free variables and their types
    int nx = 0, np=0;
    for (int i=0; i<n_; ++i) {
      if (lbx[i]==ubx[i]) {
        // Save parameter
        p_[np] = lbx[i];

        // Add contribution to objective
        objParam += g[i]*p_[np];

        // Save index
        x_index_[i] = -1-np++;

      } else {
        // True free variable
        if (lbx[i]==-numeric_limits<double>::infinity()) {
          xlow_[nx] = 0;
          ixlow_[nx] = 0;
        } else {
          xlow_[nx] = lbx[i];
          ixlow_[nx] = 1;
        }
        if (ubx[i]==numeric_limits<double>::infinity()) {
          xupp_[nx] = 0;
          ixupp_[nx] = 0;
        } else {
          xupp_[nx] = ubx[i];
          ixupp_[nx] = 1;
        }
        c_[nx] = g[i];
        x_index_[i] = nx++;
      }
    }

    // Get quadratic term
    const vector<double>& H = input(QP_SOLVER_H).data();
    const vector<int>& H_colind = input(QP_SOLVER_H).colind();
    const vector<int>& H_row = input(QP_SOLVER_H).row();
    int nnzQ = 0;
    // Loop over the columns of the quadratic term
    for (int cc=0; cc<n_; ++cc) {

      // Loop over nonzero elements of the column
      for (int el=H_colind[cc]; el<H_colind[cc+1]; ++el) {

        // Only upper triangular part
        int rr=H_row[el];
        if (rr>cc) break;

        // Get variable types
        int icc=x_index_[cc];
        int irr=x_index_[rr];

        if (icc<0) {
          if (irr<0) {
            // Add contribution to objective
            objParam += icc==irr ? H[el]*sq(p_[-1-icc])/2 : H[el]*p_[-1-irr]*p_[-1-icc];
          } else {
            // Add contribution to gradient term
            c_[irr] += H[el]*p_[-1-icc];
          }
        } else {
          if (irr<0) {
            // Add contribution to gradient term
            c_[icc] += H[el]*p_[-1-irr];
          } else {
            // Add to sparsity pattern
            irowQ_[nnzQ] = icc; // row-major --> indices swapped
            jcolQ_[nnzQ] = irr; // row-major --> indices swapped
            dQ_[nnzQ++] = H[el];
          }
        }
      }
    }

    // Get the transpose of the sparsity pattern to be able to loop over the constraints
    const vector<double>& A = input(QP_SOLVER_A).data();
    const vector<int>& A_colind = input(QP_SOLVER_A).colind();
    const vector<int>& A_row = input(QP_SOLVER_A).row();
    vector<double>& AT = AT_.data();
    const vector<int>& AT_colind = AT_.colind();
    const vector<int>& AT_row = AT_.row();
    std::copy(AT_colind.begin(), AT_colind.begin()+nc_, AT_tmp_.begin());
    for (int cc=0; cc<n_; ++cc) {
      for (int el=A_colind[cc]; el<A_colind[cc+1]; ++el) {
        int rr=A_row[el];
        int elT = AT_tmp_[rr]++;
        AT[elT] = A[el];
      }
    }

    // Loop over constraints
    int nA=0, nC=0, /*mz=0, */ nnzA=0, nnzC=0;
    for (int j=0; j<nc_; ++j) {
      if (lba[j] == -numeric_limits<double>::infinity() &&
          uba[j] ==  numeric_limits<double>::infinity()) {
        // Redundant constraint
        c_index_[j] = 0;
      } else if (lba[j]==uba[j]) {
        // Equality constraint
        bA_[nA] = lba[j];

        // Add to A
        for (int el=AT_colind[j]; el<AT_colind[j+1]; ++el) {
          int i=AT_row[el];
          if (x_index_[i]<0) {
            // Parameter
            bA_[nA] -= AT[el]*p_[-x_index_[i]-1];
          } else {
            // Free variable
            irowA_[nnzA] = nA;
            jcolA_[nnzA] = x_index_[i];
            dA_[nnzA++] = AT[el];
          }
        }
        c_index_[j] = -1-nA++;
      } else {
        // Inequality constraint
        if (lba[j]==-numeric_limits<double>::infinity()) {
          clow_[nC] = 0;
          iclow_[nC] = 0;
        } else {
          clow_[nC] = lba[j];
          iclow_[nC] = 1;
        }
        if (uba[j]==numeric_limits<double>::infinity()) {
          cupp_[nC] = 0;
          icupp_[nC] = 0;
        } else {
          cupp_[nC] = uba[j];
          icupp_[nC] = 1;
        }

        // Add to C
        for (int el=AT_colind[j]; el<AT_colind[j+1]; ++el) {
          int i=AT_row[el];
          if (x_index_[i]<0) {
            // Parameter
            if (iclow_[nC]==1) clow_[nC] -= AT[el]*p_[-x_index_[i]-1];
            if (icupp_[nC]==1) cupp_[nC] -= AT[el]*p_[-x_index_[i]-1];
          } else {
            // Free variable
            irowC_[nnzC] = nC;
            jcolC_[nnzC] = x_index_[i];
            dC_[nnzC++] = AT[el];
          }
        }
        c_index_[j] = 1+nC++;
      }
    }

    // Print problem
    if (verbose()) {
      cout << "Attempting to solve the QP:" << endl;
      cout << "   minimize    1/2*x'*Q*x + c'*x" << endl;
      cout << "   subject to  A*x = b,  d <= C*x <= f, l <= x <= u" << endl;
      cout << "with" << endl;
      cout << "Q = " <<
          tril2symm(DMatrix::triplet(vector<int>(irowQ_.begin(), irowQ_.begin()+nnzQ),
                                     vector<int>(jcolQ_.begin(), jcolQ_.begin()+nnzQ),
                                     vector<double>(dQ_.begin(), dQ_.begin()+nnzQ),
                                     nx, nx)) << endl;
      cout << "c = " << vector<double>(c_.begin(), c_.begin()+nx) << endl;
      cout << "A = " << DMatrix::triplet(vector<int>(irowA_.begin(),
                                                     irowA_.begin()+nnzA),
                                        vector<int>(jcolA_.begin(),
                                                    jcolA_.begin()+nnzA),
                                        vector<double>(dA_.begin(),
                                                       dA_.begin()+nnzA),
                                         nA, nx) << endl;
      cout << "b = " << vector<double>(bA_.begin(), bA_.begin()+nA) << endl;
      cout << "C = " << DMatrix::triplet(vector<int>(irowC_.begin(),
                                                     irowC_.begin()+nnzC),
                                        vector<int>(jcolC_.begin(),
                                                    jcolC_.begin()+nnzC),
                                        vector<double>(dC_.begin(),
                                                       dC_.begin()+nnzC),
                                         nC, nx) << endl;
      cout << "d = " << printBounds(clow_, iclow_, nC, "-") << endl;
      cout << "f = " << printBounds(cupp_, icupp_, nC, "") << endl;
      cout << "l = " << printBounds(xlow_, ixlow_, nx, "-") << endl;
      cout << "u = " << printBounds(xupp_, ixupp_, nx, "") << endl;
    }

    // Reset the solution
    fill(x_.begin(), x_.end(), 0);
    fill(gamma_.begin(), gamma_.end(), 0);
    fill(phi_.begin(), phi_.end(), 0);
    fill(y_.begin(), y_.end(), 0);
    fill(z_.begin(), z_.end(), 0);
    fill(lambda_.begin(), lambda_.end(), 0);
    fill(pi_.begin(), pi_.end(), 0);

    // Solve the QP
    double objectiveValue;

    int ierr;
    if (false) { // Use C interface
      // TODO(jgillis): Change to qpsolvehb, see OOQP users guide
      qpsolvesp(getPtr(c_), nx,
                getPtr(irowQ_),  nnzQ, getPtr(jcolQ_), getPtr(dQ_),
                getPtr(xlow_), getPtr(ixlow_),
                getPtr(xupp_), getPtr(ixupp_),
                getPtr(irowA_), nnzA, getPtr(jcolA_), getPtr(dA_),
                getPtr(bA_), nA,
                getPtr(irowC_), nnzC, getPtr(jcolC_), getPtr(dC_),
                getPtr(clow_), nC, getPtr(iclow_),
                getPtr(cupp_), getPtr(icupp_),
                getPtr(x_), getPtr(gamma_), getPtr(phi_),
                getPtr(y_),
                getPtr(z_), getPtr(lambda_), getPtr(pi_),
                &objectiveValue,
                print_level_, &ierr);
    } else { // Use C++ interface
      ierr=0;
      // All OOQP related allocations in evaluate

      std::vector<int> krowQ(nx+1);
      std::vector<int> krowA(nA+1);
      std::vector<int> krowC(nC+1);

      //int status_code = 0;
      makehb(getPtr(irowQ_), nnzQ, getPtr(krowQ), nx, &ierr);
      if (ierr == 0) makehb(getPtr(irowA_), nnzA, getPtr(krowA), nA, &ierr);
      if (ierr == 0) makehb(getPtr(irowC_), nnzC, getPtr(krowC), nC, &ierr);

      if (ierr == 0) {
        QpGenContext ctx;

        QpGenHbGondzioSetup(getPtr(c_), nx, getPtr(krowQ), getPtr(jcolQ_), getPtr(dQ_),
        getPtr(xlow_), getPtr(ixlow_), getPtr(xupp_), getPtr(ixupp_),
        getPtr(krowA), nA, getPtr(jcolA_), getPtr(dA_), getPtr(bA_),
        getPtr(krowC), nC, getPtr(jcolC_), getPtr(dC_),
        getPtr(clow_), getPtr(iclow_), getPtr(cupp_), getPtr(icupp_), &ctx,
        &ierr);
        if (ierr == 0) {
          Solver* solver = static_cast<Solver *>(ctx.solver);
          gOoqpPrintLevel = print_level_;
          solver->monitorSelf();
          solver->setMuTol(mutol_);
          solver->setMuTol(mutol_);

          QpGenFinish(&ctx, getPtr(x_), getPtr(gamma_), getPtr(phi_),
              getPtr(y_), getPtr(z_), getPtr(lambda_), getPtr(pi_),
              &objectiveValue, &ierr);
        }

        QpGenCleanup(&ctx);
      }
    }

    if (ierr>0) {
      casadi_warning("Unable to solve problem: " << errFlag(ierr));
    } else if (ierr<0) {
      casadi_error("Fatal error: " << errFlag(ierr));
    }

    // Save optimal cost
    output(QP_SOLVER_COST).set(objectiveValue + objParam);

    // Save primal solution
    vector<double>& x = output(QP_SOLVER_X).data();
    for (int i=0; i<n_; ++i) {
      int ii = x_index_[i];
      if (ii<0) {
        x[i] = p_[-1-ii];
      } else {
        x[i] = x_[ii];
      }
    }

    // Save dual solution (linear bounds)
    vector<double>& lam_a = output(QP_SOLVER_LAM_A).data();
    for (int j=0; j<nc_; ++j) {
      int jj = c_index_[j];
      if (jj==0) {
        lam_a[j] = 0;
      } else if (jj<0) {
        lam_a[j] = -y_[-1-jj];
      } else {
        lam_a[j] = pi_[-1+jj]-lambda_[-1+jj];
      }
    }

    // Save dual solution (simple bounds)
    vector<double>& lam_x = output(QP_SOLVER_LAM_X).data();
    for (int i=0; i<n_; ++i) {
      int ii = x_index_[i];
      if (ii<0) {
        // The dual solution for the fixed parameters follows from the KKT conditions
        lam_x[i] = -g[i];
        for (int el=H_colind[i]; el<H_colind[i+1]; ++el) {
          int j=H_row[el];
          lam_x[i] -= H[el]*x[j];
        }
        for (int el=A_colind[i]; el<A_colind[i+1]; ++el) {
          int j=A_row[el];
          lam_x[i] -= A[el]*lam_a[j];
        }
      } else {
        lam_x[i] = phi_[ii]-gamma_[ii];
      }
    }
  }

  const char* OoqpInterface::errFlag(int flag) {
    // Find the error
    //const char* msg;
    switch (flag) {
    case SUCCESSFUL_TERMINATION: return  "SUCCESSFUL_TERMINATION";
    case NOT_FINISHED:           return  "NOT_FINISHED";
    case MAX_ITS_EXCEEDED:       return  "MAX_ITS_EXCEEDED";
    case INFEASIBLE:             return  "INFEASIBLE";
    case UNKNOWN:                return  "UNKNOWN";
    default:                     return  "N/A";
    }
  }

  std::string OoqpInterface::printBounds(const std::vector<double>& b,
                                        const std::vector<char>& ib, int n, const char *sign) {
    stringstream ss;
    ss << "[";
    for (int i=0; i<n; ++i) {
      if (i!=0) ss << ", ";
      if (ib[i]==0) {
        ss << sign << "inf";
      } else {
        ss << b[i];
      }
    }
    ss << "]";
    return ss.str();
  }


} // namespace casadi

