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


#include "symbolic_qr.hpp"

#ifdef WITH_DL
#include <cstdlib>
#endif // WITH_DL

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINEARSOLVER_SYMBOLICQR_EXPORT
  casadi_register_linearsolver_symbolicqr(LinearSolverInternal::Plugin* plugin) {
    plugin->creator = SymbolicQr::creator;
    plugin->name = "symbolicqr";
    plugin->doc = SymbolicQr::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_LINEARSOLVER_SYMBOLICQR_EXPORT casadi_load_linearsolver_symbolicqr() {
    LinearSolverInternal::registerPlugin(casadi_register_linearsolver_symbolicqr);
  }

  SymbolicQr::SymbolicQr(const std::string& name, const Sparsity& sparsity, int nrhs) :
    LinearSolverInternal(name, sparsity, nrhs) {
    addOption("codegen",           OT_BOOLEAN,  false,               "C-code generation");
    addOption("compiler",          OT_STRING,    "gcc -fPIC -O2",
              "Compiler command to be used for compiling generated code");
  }

  SymbolicQr::~SymbolicQr() {
  }

  void SymbolicQr::init() {
    // Call the base class initializer
    LinearSolverInternal::init();

    // Read options
    bool codegen = option("codegen");
    string compiler = option("compiler");

    // Make sure that command processor is available
    if (codegen) {
#ifdef WITH_DL
      int flag = system(static_cast<const char*>(0));
      casadi_assert_message(flag!=0, "No command procesor available");
#else // WITH_DL
      casadi_error("Codegen requires CasADi to be compiled with option \"WITH_DL\" enabled");
#endif // WITH_DL
    }

    // Symbolic expression for A
    SX A = SX::sym("A", input(0).sparsity());

    // Get the inverted column permutation
    std::vector<int> inv_colperm(colperm_.size());
    for (int k=0; k<colperm_.size(); ++k)
      inv_colperm[colperm_[k]] = k;

    // Get the inverted row permutation
    std::vector<int> inv_rowperm(rowperm_.size());
    for (int k=0; k<rowperm_.size(); ++k)
      inv_rowperm[rowperm_[k]] = k;

    // Permute the linear system
    SX Aperm = A(rowperm_, colperm_);

    // Generate the QR factorization function
    vector<SX> QR(2);
    qr(Aperm, QR[0], QR[1]);
    Function fact_fcn("QR_fact", {A}, QR);

    // Optionally generate c code and load as DLL
    if (codegen) {
      stringstream ss;
      ss << "symbolic_qr_fact_fcn_" << this;
      fact_fcn_ = dynamicCompilation(fact_fcn, ss.str(),
                                     "Symbolic QR factorization function", compiler);
    } else {
      fact_fcn_ = fact_fcn;
    }

    // Symbolic expressions for solve function
    SX Q = SX::sym("Q", QR[0].sparsity());
    SX R = SX::sym("R", QR[1].sparsity());
    SX b = SX::sym("b", input(1).size1(), 1);

    // Solve non-transposed
    // We have Pb' * Q * R * Px * x = b <=> x = Px' * inv(R) * Q' * Pb * b

    // Permute the right hand sides
    SX bperm = b(rowperm_, ALL);

    // Solve the factorized system
    SX xperm = R.zz_solve(mul(Q.T(), bperm));

    // Permute back the solution
    SX x = xperm(inv_colperm, ALL);

    // Generate the QR solve function
    vector<SX> solv_in = {Q, R, b};
    Function solv_fcn("QR_solv", solv_in, {x});

    // Optionally generate c code and load as DLL
    if (codegen) {
      stringstream ss;
      ss << "symbolic_qr_solv_fcn_N_" << this;
      solv_fcn_N_ = dynamicCompilation(solv_fcn, ss.str(), "QR_solv_N", compiler);
    } else {
      solv_fcn_N_ = solv_fcn;
    }

    // Solve transposed
    // We have (Pb' * Q * R * Px)' * x = b
    // <=> Px' * R' * Q' * Pb * x = b
    // <=> x = Pb' * Q * inv(R') * Px * b

    // Permute the right hand side
    bperm = b(colperm_, ALL);

    // Solve the factorized system
    xperm = mul(Q, R.T().zz_solve(bperm));

    // Permute back the solution
    x = xperm(inv_rowperm, ALL);

    // Mofify the QR solve function
    solv_fcn = Function("QR_solv_T", solv_in, {x});

    // Optionally generate c code and load as DLL
    if (codegen) {
      stringstream ss;
      ss << "symbolic_qr_solv_fcn_T_" << this;
      solv_fcn_T_ = dynamicCompilation(solv_fcn, ss.str(), "QR_solv_T", compiler);
    } else {
      solv_fcn_T_ = solv_fcn;
    }

    // Allocate storage for QR factorization
    Q_ = DMatrix::zeros(Q.sparsity());
    R_ = DMatrix::zeros(R.sparsity());
  }

  void SymbolicQr::linsol_prepare() {
    // Factorize
    fact_fcn_.setInput(input(LINSOL_A));
    fact_fcn_.evaluate();
    fact_fcn_.getOutput(Q_, 0);
    fact_fcn_.getOutput(R_, 1);
  }

  void SymbolicQr::linsol_solve(double* x, int nrhs, bool transpose) {
    // Select solve function
    Function& solv = transpose ? solv_fcn_T_ : solv_fcn_N_;

    // Pass QR factorization
    solv.setInput(Q_, 0);
    solv.setInput(R_, 1);

    // Solve for all right hand sides
    for (int i=0; i<nrhs; ++i) {
      solv.setInputNZ(x, 2);
      solv.evaluate();
      solv.getOutputNZ(x);
      x += solv.nnz_out(0);
    }
  }

  void SymbolicQr::linsol_evalSX(void* mem, const SXElem** arg, SXElem** res,
                                 int* iw, SXElem* w, bool tr, int nrhs) {
    casadi_assert(arg[0]!=0);
    casadi_assert(arg[1]!=0);
    casadi_assert(res[0]!=0);

    // Get A and factorize it
    SX A = SX::zeros(input(LINSOL_A).sparsity());
    copy(arg[1], arg[1]+A.nnz(), A->begin());
    vector<SX> v = fact_fcn_(A);

    // Select solve function
    Function& solv = tr ? solv_fcn_T_ : solv_fcn_N_;

    // Solve for every right hand side
    v.push_back(SX::zeros(A.size1()));
    const SXElem* a=arg[0];
    SXElem* r=res[0];
    for (int i=0; i<nrhs; ++i) {
      copy(a, a+v[2].nnz(), v[2]->begin());
      SX rr = solv(v).at(0);
      copy(rr->begin(), rr->end(), r);
      r += rr.nnz();
    }
  }

  void SymbolicQr::generateDeclarations(CodeGenerator& g) const {

    // Generate code for the embedded functions
    fact_fcn_->addDependency(g);
    solv_fcn_N_->addDependency(g);
    solv_fcn_T_->addDependency(g);
  }

  void SymbolicQr::generateBody(CodeGenerator& g) const {
    casadi_error("Code generation for SymbolicQR still experimental");
#if 0
    // Data structures to hold A, Q and R
    g.body
       << "  static int prepared = 0;" << endl
       << "  static d A[" << input(LINSOL_A).nnz() << "];" << endl
       << "  static d Q[" << Q_.nnz() << "];" << endl
       << "  static d R[" << R_.nnz() << "];" << endl;

    // Check if the factorization is up-to-date
    g.body
      << "  int i;" << endl
      << "  for (i=0; prepared && i<" << input(LINSOL_A).nnz()
      << "; ++i) prepared=A[i]!=x0[i];" << endl;

    // Factorize if needed
    int fact_ind = g.getDependency(fact_fcn_);
    g.body
      << "  if (!prepared) {" << endl
      << "    for (i=0; i<" << input(LINSOL_A).nnz() << "; ++i) A[i]=x0[i];" << endl
      << "    f" << fact_ind << "(A, Q, R);" << endl
      << "    prepared = 1;" << endl
      << "  }" << endl;

    // Solve
    int solv_ind_N = g.getDependency(solv_fcn_N_);
    int neq = input(LINSOL_B).size1();
    int nrhs = input(LINSOL_B).size2();
    g.body
      << "  for (i=0; i<" << nrhs << "; ++i) {" << endl
      << "    f" << solv_ind_N << "(Q, R, x1, r0);" << endl
      << "    x1+=" << neq << "; r0+=" << neq << ";" << endl
      << "  }" << endl;
#endif
  }

} // namespace casadi
