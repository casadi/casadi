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

#include "clp_interface.hpp"

namespace casadi {

  using namespace std;

  extern "C"
  int CASADI_CONIC_CLP_EXPORT
  casadi_register_conic_clp(Conic::Plugin* plugin) {
    plugin->creator = ClpInterface::creator;
    plugin->name = "clp";
    plugin->doc = ClpInterface::meta_doc.c_str();
    plugin->version = 30;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_CLP_EXPORT casadi_load_conic_clp() {
    Conic::registerPlugin(casadi_register_conic_clp);
  }

  ClpInterface::ClpInterface(const std::string& name,
                             const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
  }

  void ClpInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Conic::init(opts);

    // Default options
    casadi_assert_message(H_.nnz()==0, "Not an LP");

    // Read options
    //for (auto&& op : opts) {
    //}

    // Allocate work vectors
    alloc_w(nx_, true); // g
    alloc_w(nx_, true); // lbx
    alloc_w(nx_, true); // ubx
    alloc_w(na_, true); // lba
    alloc_w(na_, true); // uba
    alloc_w(nnz_in(CONIC_H), true); // H
    alloc_w(nnz_in(CONIC_A), true); // A
  }

  void ClpInterface::init_memory(void* mem) const {
    //auto m = static_cast<ClpMemory*>(mem);
  }

  void ClpInterface::
  eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    //auto m = static_cast<ClpMemory*>(mem);

    if (inputs_check_) {
      checkInputs(arg[CONIC_LBX], arg[CONIC_UBX], arg[CONIC_LBA], arg[CONIC_UBA]);
    }

    // Get inputs
    double* g=w; w += nx_;
    casadi_copy(arg[CONIC_G], nx_, g);
    double* lbx=w; w += nx_;
    casadi_copy(arg[CONIC_LBX], nx_, lbx);
    double* ubx=w; w += nx_;
    casadi_copy(arg[CONIC_UBX], nx_, ubx);
    double* lba=w; w += na_;
    casadi_copy(arg[CONIC_LBA], na_, lba);
    double* uba=w; w += na_;
    casadi_copy(arg[CONIC_UBA], na_, uba);
    double* H=w; w += nnz_in(CONIC_H);
    casadi_copy(arg[CONIC_H], nnz_in(CONIC_H), H);
    double* A=w; w += nnz_in(CONIC_A);
    casadi_copy(arg[CONIC_A], nnz_in(CONIC_A), A);

    // Create model
    ClpSimplex model;
    model.loadProblem(A_.size2(), A_.size1(), A_.colind(), A_.row(), A,
                      lbx, ubx, g, lba, uba, nullptr);

    // Solve the problem using the primal simplex algorithm
    model.primal();

    // Primal solution
    double* x = model.primalColumnSolution();
    casadi_copy(x, nx_, res[CONIC_X]);

    // Dual solution (x)
    double* minus_lam_x = model.dualColumnSolution();
    if (res[CONIC_LAM_X]) {
      casadi_copy(minus_lam_x, nx_, res[CONIC_LAM_X]);
      casadi_scal(nx_, -1., res[CONIC_LAM_X]);
    }

    // Dual solution (A)
    double* minus_lam_a = model.dualRowSolution();
    if (res[CONIC_LAM_A]) {
      casadi_copy(minus_lam_a, na_, res[CONIC_LAM_A]);
      casadi_scal(na_, -1., res[CONIC_LAM_A]);
    }

    // Optimal cost
    double f = model.rawObjectiveValue();
    if (res[CONIC_COST]) *res[CONIC_COST] = f;
  }

  ClpInterface::~ClpInterface() {
    clear_memory();
  }

  ClpMemory::ClpMemory() {
  }

  ClpMemory::~ClpMemory() {
  }

} // end namespace casadi
