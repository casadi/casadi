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
#include "casadi/core/function/qpsol.hpp"
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
  int CASADI_QPSOL_OOQP_EXPORT
  casadi_register_qpsol_ooqp(Qpsol::Plugin* plugin) {
    plugin->creator = OoqpInterface::creator;
    plugin->name = "ooqp";
    plugin->doc = OoqpInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_QPSOL_OOQP_EXPORT casadi_load_qpsol_ooqp() {
    Qpsol::registerPlugin(casadi_register_qpsol_ooqp);
  }

  OoqpInterface::OoqpInterface(const std::string& name,
                               const std::map<std::string, Sparsity>& st)
    : Qpsol(name, st) {

    addOption("print_level", OT_INT, 0,
              "Print level. OOQP listens to print_level 0, 10 and 100");
    addOption("mutol", OT_DOUBLE, 1e-8, "tolerance as provided with setMuTol to OOQP");
    addOption("artol", OT_DOUBLE, 1e-8, "tolerance as provided with setArTol to OOQP");
  }

  OoqpInterface::~OoqpInterface() {
  }

  void OoqpInterface::init(const Dict& opts) {
    // Initialize the base classes
    Qpsol::init(opts);

    // Read options
    print_level_ = option("print_level");
    mutol_ = option("mutol");
    artol_ = option("artol");

    // Allocate memory for problem
    nQ_ = sparsity_in(QPSOL_H).nnz_upper();
    nA_ = nnz_in(QPSOL_A);
    nH_ = nnz_in(QPSOL_H);
    spAT_ = sparsity_in(QPSOL_A).T();

    // Allocate work vectors
    alloc_w(n_, true); // g
    alloc_w(n_, true); // lbx
    alloc_w(n_, true); // ubx
    alloc_w(nc_, true); // lba
    alloc_w(nc_, true); // uba
    alloc_w(nH_, true); // H
    alloc_w(nA_, true); // A
    alloc_w(n_, true); // c_
    alloc_w(nc_, true); // bA_
    alloc_w(n_, true); // xlow_
    alloc_w(n_, true); // xupp_
    alloc_w(nc_, true); // clow_
    alloc_w(nc_, true); // cupp_
    alloc_w(n_, true); // x_
    alloc_w(n_, true); // gamma_
    alloc_w(n_, true); // phi_
    alloc_w(nc_, true); // y_
    alloc_w(nc_, true); // z_
    alloc_w(nc_, true); // lambda_
    alloc_w(nc_, true); // pi_
    alloc_iw(n_, true); // ixlow_
    alloc_iw(n_, true); // ixupp_
    alloc_iw(nc_, true); // iclow_
    alloc_iw(nc_, true); // icupp_
    alloc_w(nQ_, true); // dQ_
    alloc_w(nA_, true); // dA_
    alloc_w(nA_, true); // dC_
    alloc_iw(nQ_, true); // irowQ_
    alloc_iw(nQ_, true); // jcolQ_
    alloc_iw(nA_, true); // irowA_
    alloc_iw(nA_, true); // jcolA_
    alloc_iw(nA_, true); // irowC_
    alloc_iw(nA_, true); // jcolC_
    alloc_iw(n_, true); // x_index_
    alloc_iw(nc_, true); // c_index_
    alloc_w(n_, true); // p_
    alloc_w(nA_, true); // AT
    alloc_iw(nc_); // casadi_trans
  }

  void OoqpInterface::
  eval(Memory& mem, const double** arg, double** res, int* iw, double* w) const {
    if (inputs_check_) {
      checkInputs(arg[QPSOL_LBX], arg[QPSOL_UBX], arg[QPSOL_LBA], arg[QPSOL_UBA]);
    }

    // Get problem data
    double* g=w; w += n_;
    casadi_copy(arg[QPSOL_G], n_, g);
    double* lbx=w; w += n_;
    casadi_copy(arg[QPSOL_LBX], n_, lbx);
    double* ubx=w; w += n_;
    casadi_copy(arg[QPSOL_UBX], n_, ubx);
    double* lba=w; w += nc_;
    casadi_copy(arg[QPSOL_LBA], nc_, lba);
    double* uba=w; w += nc_;
    casadi_copy(arg[QPSOL_UBA], nc_, uba);
    double* H=w; w += nnz_in(QPSOL_H);
    casadi_copy(arg[QPSOL_H], nnz_in(QPSOL_H), H);
    double* A=w; w += nnz_in(QPSOL_A);
    casadi_copy(arg[QPSOL_A], nnz_in(QPSOL_A), A);

    // Temporary memory
    double* c_ = w; w += n_;
    double* bA_ = w; w += nc_;
    double* xlow_ = w; w += n_;
    double* xupp_ = w; w += n_;
    double* clow_ = w; w += nc_;
    double* cupp_ = w; w += nc_;
    double* x_ = w; w += n_;
    double* gamma_ = w; w += n_;
    double* phi_ = w; w += n_;
    double* y_ = w; w += nc_;
    double* z_ = w; w += nc_;
    double* lambda_ = w; w += nc_;
    double* pi_ = w; w += nc_;
    char* ixlow_ = reinterpret_cast<char*>(iw); iw += n_;
    char* ixupp_ = reinterpret_cast<char*>(iw); iw += n_;
    char* iclow_ = reinterpret_cast<char*>(iw); iw += nc_;
    char* icupp_ = reinterpret_cast<char*>(iw); iw += nc_;
    double* dQ_ = w; w += nQ_;
    double* dA_ = w; w += nA_;
    double* dC_ = w; w += nA_;
    int* irowQ_ = iw; iw += nQ_;
    int* jcolQ_ = iw; iw += nQ_;
    int* irowA_ = iw; iw += nA_;
    int* jcolA_ = iw; iw += nA_;
    int* irowC_ = iw; iw += nA_;
    int* jcolC_ = iw; iw += nA_;
    int* x_index_ = iw; iw += n_;
    int* c_index_ = iw; iw += nc_;
    double* p_ = w; w += n_;
    double* AT = w; w += nA_;

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
    const int* H_colind = sparsity_in(QPSOL_H).colind();
    const int* H_row = sparsity_in(QPSOL_H).row();
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
    casadi_trans(A, sparsity_in(QPSOL_A), AT, spAT_, iw);

    // Loop over constraints
    const int* A_colind = sparsity_in(QPSOL_A).colind();
    const int* A_row = sparsity_in(QPSOL_A).row();
    const int* AT_colind = spAT_.colind();
    const int* AT_row = spAT_.row();
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

    // Reset the solution
    casadi_fill(x_, n_, 0.);
    casadi_fill(gamma_, n_, 0.);
    casadi_fill(phi_, n_, 0.);
    casadi_fill(y_, nc_, 0.);
    casadi_fill(z_, nc_, 0.);
    casadi_fill(lambda_, nc_, 0.);
    casadi_fill(pi_, nc_, 0.);

    // Solve the QP
    double objectiveValue;

    int ierr;
    if (false) { // Use C interface
      // TODO(jgillis): Change to qpsolvehb, see OOQP users guide
      qpsolvesp(c_, nx,
                irowQ_,  nnzQ, jcolQ_, dQ_,
                xlow_, ixlow_,
                xupp_, ixupp_,
                irowA_, nnzA, jcolA_, dA_,
                bA_, nA,
                irowC_, nnzC, jcolC_, dC_,
                clow_, nC, iclow_,
                cupp_, icupp_,
                x_, gamma_, phi_,
                y_,
                z_, lambda_, pi_,
                &objectiveValue,
                print_level_, &ierr);
    } else { // Use C++ interface
      ierr=0;
      // All OOQP related allocations in evaluate

      std::vector<int> krowQ(nx+1);
      std::vector<int> krowA(nA+1);
      std::vector<int> krowC(nC+1);

      //int status_code = 0;
      makehb(irowQ_, nnzQ, get_ptr(krowQ), nx, &ierr);
      if (ierr == 0) makehb(irowA_, nnzA, get_ptr(krowA), nA, &ierr);
      if (ierr == 0) makehb(irowC_, nnzC, get_ptr(krowC), nC, &ierr);

      if (ierr == 0) {
        QpGenContext ctx;

        QpGenHbGondzioSetup(c_, nx, get_ptr(krowQ), jcolQ_, dQ_,
                            xlow_, ixlow_, xupp_, ixupp_,
                            get_ptr(krowA), nA, jcolA_, dA_, bA_,
                            get_ptr(krowC), nC, jcolC_, dC_,
                            clow_, iclow_, cupp_, icupp_, &ctx,
                            &ierr);
        if (ierr == 0) {
          Solver* solver = static_cast<Solver *>(ctx.solver);
          gOoqpPrintLevel = print_level_;
          solver->monitorSelf();
          solver->setMuTol(mutol_);
          solver->setMuTol(mutol_);

          QpGenFinish(&ctx, x_, gamma_, phi_,
                      y_, z_, lambda_, pi_,
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

    // Retrieve eliminated decision variables
    for (int i=n_-1; i>=0; --i) {
      int ii = x_index_[i];
      if (ii<0) {
        x_[i] = p_[-1-ii];
      } else {
        x_[i] = x_[ii];
      }
    }

    // Retreive eliminated dual variables (linear bounds)
    for (int j=nc_-1; j>=0; --j) {
      int jj = c_index_[j];
      if (jj==0) {
        lambda_[j] = 0;
      } else if (jj<0) {
        lambda_[j] = -y_[-1-jj];
      } else {
        lambda_[j] = pi_[-1+jj]-lambda_[-1+jj];
      }
    }

    // Retreive eliminated dual variables (simple bounds)
    for (int i=n_-1; i>=0; --i) {
      int ii = x_index_[i];
      if (ii<0) {
        // The dual solution for the fixed parameters follows from the KKT conditions
        gamma_[i] = -g[i];
        for (int el=H_colind[i]; el<H_colind[i+1]; ++el) {
          int j=H_row[el];
          gamma_[i] -= H[el]*x_[j];
        }
        for (int el=A_colind[i]; el<A_colind[i+1]; ++el) {
          int j=A_row[el];
          gamma_[i] -= A[el]*lambda_[j];
        }
      } else {
        gamma_[i] = phi_[ii]-gamma_[ii];
      }
    }

    // Save optimal cost
    if (res[QPSOL_COST]) *res[QPSOL_COST] = objectiveValue + objParam;

    // Save primal solution
    casadi_copy(x_, n_, res[QPSOL_X]);

    // Save dual solution (linear bounds)
    casadi_copy(lambda_, nc_, res[QPSOL_LAM_A]);

    // Save dual solution (simple bounds)
    casadi_copy(gamma_, n_, res[QPSOL_LAM_X]);
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
