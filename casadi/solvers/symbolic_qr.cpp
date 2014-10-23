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
#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/function/sx_function.hpp"

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
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_LINEARSOLVER_SYMBOLICQR_EXPORT casadi_load_linearsolver_symbolicqr() {
    LinearSolverInternal::registerPlugin(casadi_register_linearsolver_symbolicqr);
  }

  SymbolicQr::SymbolicQr(const Sparsity& sparsity, int nrhs) :
      LinearSolverInternal(sparsity, nrhs) {
    addOption("codegen",           OT_BOOLEAN,  false,               "C-code generation");
    addOption("compiler",          OT_STRING,    "gcc -fPIC -O2",
              "Compiler command to be used for compiling generated code");
  }

  SymbolicQr::~SymbolicQr() {
  }

  void SymbolicQr::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    LinearSolverInternal::deepCopyMembers(already_copied);
    fact_fcn_ = deepcopy(fact_fcn_, already_copied);
    solv_fcn_N_ = deepcopy(solv_fcn_N_, already_copied);
    solv_fcn_T_ = deepcopy(solv_fcn_T_, already_copied);
  }

  void SymbolicQr::init() {
    // Call the base class initializer
    LinearSolverInternal::init();

    // Read options
    bool codegen = getOption("codegen");
    string compiler = getOption("compiler");

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
    SXFunction fact_fcn(A, QR);

    // Optionally generate c code and load as DLL
    if (codegen) {
      stringstream ss;
      ss << "symbolic_qr_fact_fcn_" << this;
      fact_fcn_ = dynamicCompilation(fact_fcn, ss.str(),
                                     "Symbolic QR factorization function", compiler);
    } else {
      fact_fcn_ = fact_fcn;
    }

    // Initialize factorization function
    fact_fcn_.setOption("name", "QR_fact");
    fact_fcn_.init();

    // Symbolic expressions for solve function
    SX Q = SX::sym("Q", QR[0].sparsity());
    SX R = SX::sym("R", QR[1].sparsity());
    SX b = SX::sym("b", input(1).size1(), 1);

    // Solve non-transposed
    // We have Pb' * Q * R * Px * x = b <=> x = Px' * inv(R) * Q' * Pb * b

    // Permute the right hand sides
    SX bperm = b(rowperm_, ALL);

    // Solve the factorized system
    SX xperm = casadi::solve(R, mul(Q.T(), bperm));

    // Permute back the solution
    SX x = xperm(inv_colperm, ALL);

    // Generate the QR solve function
    vector<SX> solv_in(3);
    solv_in[0] = Q;
    solv_in[1] = R;
    solv_in[2] = b;
    SXFunction solv_fcn(solv_in, x);

    // Optionally generate c code and load as DLL
    if (codegen) {
      stringstream ss;
      ss << "symbolic_qr_solv_fcn_N_" << this;
      solv_fcn_N_ = dynamicCompilation(solv_fcn, ss.str(), "QR_solv_N", compiler);
    } else {
      solv_fcn_N_ = solv_fcn;
    }

    // Initialize solve function
    solv_fcn_N_.setOption("name", "QR_solv");
    solv_fcn_N_.init();

    // Solve transposed
    // We have (Pb' * Q * R * Px)' * x = b
    // <=> Px' * R' * Q' * Pb * x = b
    // <=> x = Pb' * Q * inv(R') * Px * b

    // Permute the right hand side
    bperm = b(colperm_, ALL);

    // Solve the factorized system
    xperm = mul(Q, casadi::solve(R.T(), bperm));

    // Permute back the solution
    x = xperm(inv_rowperm, ALL);

    // Mofify the QR solve function
    solv_fcn = SXFunction(solv_in, x);

    // Optionally generate c code and load as DLL
    if (codegen) {
      stringstream ss;
      ss << "symbolic_qr_solv_fcn_T_" << this;
      solv_fcn_T_ = dynamicCompilation(solv_fcn, ss.str(), "QR_solv_T", compiler);
    } else {
      solv_fcn_T_ = solv_fcn;
    }

    // Initialize solve function
    solv_fcn_T_.setOption("name", "QR_solv_T");
    solv_fcn_T_.init();

    // Allocate storage for QR factorization
    Q_ = DMatrix::zeros(Q.sparsity());
    R_ = DMatrix::zeros(R.sparsity());
  }

  void SymbolicQr::prepare() {
    // Factorize
    fact_fcn_.setInput(input(LINSOL_A));
    fact_fcn_.evaluate();
    fact_fcn_.getOutput(Q_, 0);
    fact_fcn_.getOutput(R_, 1);
    prepared_ = true;
  }

  void SymbolicQr::solve(double* x, int nrhs, bool transpose) {
    // Select solve function
    Function& solv = transpose ? solv_fcn_T_ : solv_fcn_N_;

    // Pass QR factorization
    solv.setInput(Q_, 0);
    solv.setInput(R_, 1);

    // Solve for all right hand sides
    for (int i=0; i<nrhs; ++i) {
      solv.setInput(x, 2);
      solv.evaluate();
      solv.getOutput(x);
      x += solv.output().size();
    }
  }

  void SymbolicQr::evaluateSXGen(const SXPtrV& input, SXPtrV& output, bool tr) {
    // Get arguments
    casadi_assert(input.at(0)!=0);
    SX r = *input.at(0);
    casadi_assert(input.at(1)!=0);
    SX A = *input.at(1);

    // Number of right hand sides
    int nrhs = r.size2();

    // Factorize A
    vector<SX> v = fact_fcn_(A);

    // Select solve function
    Function& solv = tr ? solv_fcn_T_ : solv_fcn_N_;

    // Solve for every right hand side
    vector<SX> resv;
    v.resize(3);
    for (int i=0; i<nrhs; ++i) {
      v[2] = r(Slice(), i);
      resv.push_back(solv(v).at(0));
    }

    // Collect the right hand sides
    casadi_assert(output[0]!=0);
    *output.at(0) = horzcat(resv);
  }

  void SymbolicQr::generateDeclarations(std::ostream &stream, const std::string& type,
                                                CodeGenerator& gen) const {

    // Generate code for the embedded functions
    gen.addDependency(fact_fcn_);
    gen.addDependency(solv_fcn_N_);
    gen.addDependency(solv_fcn_T_);
  }

  void SymbolicQr::generateBody(std::ostream &stream, const std::string& type,
                                        CodeGenerator& gen) const {
    casadi_warning("Code generation for SymbolicQR still experimental");

    // Data structures to hold A, Q and R
    stream << "  static int prepared = 0;" << endl;
    stream << "  static d A[" << input(LINSOL_A).size() << "];" << endl;
    stream << "  static d Q[" << Q_.size() << "];" << endl;
    stream << "  static d R[" << R_.size() << "];" << endl;

    // Check if the factorization is up-to-date
    stream << "  int i;" << endl;
    stream << "  for (i=0; prepared && i<" << input(LINSOL_A).size()
           << "; ++i) prepared=A[i]!=x0[i];" << endl;

    // Factorize if needed
    int fact_ind = gen.getDependency(fact_fcn_);
    stream << "  if (!prepared) {" << endl;
    stream << "    for (i=0; i<" << input(LINSOL_A).size() << "; ++i) A[i]=x0[i];" << endl;
    stream << "    f" << fact_ind << "(A, Q, R);" << endl;
    stream << "    prepared = 1;" << endl;
    stream << "  }" << endl;

    // Solve
    int solv_ind_N = gen.getDependency(solv_fcn_N_);
    int neq = input(LINSOL_B).size1();
    int nrhs = input(LINSOL_B).size2();
    stream << "  for (i=0; i<" << nrhs << "; ++i) {" << endl;
    stream << "    f" << solv_ind_N << "(Q, R, x1, r0);" << endl;
    stream << "    x1+=" << neq << "; r0+=" << neq << ";" << endl;
    stream << "  }" << endl;
  }

} // namespace casadi
