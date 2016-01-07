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
  int CASADI_LINSOL_SYMBOLICQR_EXPORT
  casadi_register_linsol_symbolicqr(Linsol::Plugin* plugin) {
    plugin->creator = SymbolicQr::creator;
    plugin->name = "symbolicqr";
    plugin->doc = SymbolicQr::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_SYMBOLICQR_EXPORT casadi_load_linsol_symbolicqr() {
    Linsol::registerPlugin(casadi_register_linsol_symbolicqr);
  }

  SymbolicQr::SymbolicQr(const std::string& name, const Sparsity& sparsity, int nrhs) :
    Linsol(name, sparsity, nrhs) {
    addOption("codegen",           OT_BOOLEAN,  false,               "C-code generation");
    addOption("compiler",          OT_STRING,    "gcc -fPIC -O2",
              "Compiler command to be used for compiling generated code");
  }

  SymbolicQr::~SymbolicQr() {
  }

  void SymbolicQr::init() {
    // Call the base class initializer
    Linsol::init();

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
    SX A = SX::sym("A", sparsity_in(LINSOL_A));

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
    SX Q1, R1;
    qr(Aperm, Q1, R1);
    Function fact_fcn("QR_fact", {A}, {Q1, R1});

    // Optionally generate c code and load as DLL
    if (codegen) {
      stringstream ss;
      ss << "symbolic_qr_fact_fcn_" << this;
      fact_fcn_ = dynamicCompilation(fact_fcn, ss.str(),
                                     "Symbolic QR factorization function", compiler);
    } else {
      fact_fcn_ = fact_fcn;
    }
    alloc(fact_fcn_);

    // Symbolic expressions for solve function
    SX Q = SX::sym("Q", Q1.sparsity());
    SX R = SX::sym("R", R1.sparsity());
    SX b = SX::sym("b", size1_in(LINSOL_B), 1);

    // Solve non-transposed
    // We have Pb' * Q * R * Px * x = b <=> x = Px' * inv(R) * Q' * Pb * b

    // Permute the right hand sides
    SX bperm = b(rowperm_, ALL);

    // Solve the factorized system
    SX xperm = SX::solve(R, mtimes(Q.T(), bperm));

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
    alloc(solv_fcn_N_);

    // Solve transposed
    // We have (Pb' * Q * R * Px)' * x = b
    // <=> Px' * R' * Q' * Pb * x = b
    // <=> x = Pb' * Q * inv(R') * Px * b

    // Permute the right hand side
    bperm = b(colperm_, ALL);

    // Solve the factorized system
    xperm = mtimes(Q, SX::solve(R.T(), bperm));

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
    alloc(solv_fcn_T_);

    // Temporary storage
    alloc_w(neq_, true);
  }

  Memory* SymbolicQr::memory() const {
    SymbolicQrMemory* m = new SymbolicQrMemory();
    try {
      // Allocate storage for QR factorization
      m->q.resize(fact_fcn_.nnz_out(0));
      m->r.resize(fact_fcn_.nnz_out(1));
      return m;
    } catch (...) {
      delete m;
      return 0;
    }
  }

  void SymbolicQr::set_temp(Memory& mem, const double** arg, double** res,
                            int* iw, double* w) const {
    SymbolicQrMemory& m = dynamic_cast<SymbolicQrMemory&>(mem);
    m.arg = arg;
    m.res = res;
    m.iw = iw;
    m.w = w;
  }

  void SymbolicQr::linsol_factorize(Memory& mem, const double* A) const {
    SymbolicQrMemory& m = dynamic_cast<SymbolicQrMemory&>(mem);

    // Factorize
    fill_n(m.arg, fact_fcn_.n_in(), nullptr);
    m.arg[0] = A;
    fill_n(m.res, fact_fcn_.n_out(), nullptr);
    m.res[0] = getPtr(m.q);
    m.res[1] = getPtr(m.r);
    fact_fcn_(m.arg, m.res, m.iw, m.w);
  }

  void SymbolicQr::linsol_solve(Memory& mem, double* x, int nrhs, bool tr) const {
    SymbolicQrMemory& m = dynamic_cast<SymbolicQrMemory&>(mem);

    // Select solve function
    const Function& solv = tr ? solv_fcn_T_ : solv_fcn_N_;

    // Solve for all right hand sides
    fill_n(m.arg, solv.n_in(), nullptr);
    m.arg[0] = getPtr(m.q);
    m.arg[1] = getPtr(m.r);
    fill_n(m.res, solv.n_out(), nullptr);
    for (int i=0; i<nrhs; ++i) {
      copy_n(x, neq_, m.w); // Copy x to a temporary
      m.arg[2] = m.w;
      m.res[0] = x;
      solv(m.arg, m.res, m.iw, m.w+neq_, 0);
      x += neq_;
    }
  }

  void SymbolicQr::linsol_eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, void* mem,
                                 bool tr, int nrhs) {
    casadi_assert(arg[0]!=0);
    casadi_assert(arg[1]!=0);
    casadi_assert(res[0]!=0);

    // Get A and factorize it
    SX A = SX::zeros(sparsity_in(LINSOL_A));
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
  }

} // namespace casadi
