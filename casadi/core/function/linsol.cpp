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


#include "linsol_internal.hpp"

using namespace std;
namespace casadi {

  Linsol::Linsol() {
  }

  Linsol::Linsol(const std::string& name, const std::string& solver,
         const Sparsity& sp, const Dict& opts) {
    assignNode(LinsolInternal::getPlugin(solver).creator(name, sp));
    (*this)->construct(opts);
  }

  LinsolInternal* Linsol::operator->() {
    return static_cast<LinsolInternal*>(SharedObject::operator->());
  }

  const LinsolInternal* Linsol::operator->() const {
    return static_cast<const LinsolInternal*>(SharedObject::operator->());
  }

  bool Linsol::test_cast(const SharedObjectNode* ptr) {
    return dynamic_cast<const LinsolInternal*>(ptr)!=0;
  }

  bool Linsol::has_plugin(const std::string& name) {
    return LinsolInternal::has_plugin(name);
  }

  void Linsol::load_plugin(const std::string& name) {
    LinsolInternal::load_plugin(name);
  }

  std::string Linsol::doc(const std::string& name) {
    return LinsolInternal::getPlugin(name).doc;
  }

  std::string Linsol::plugin_name() const {
    return (*this)->plugin_name();
  }

  DM Linsol::solve(const DM& A, const DM& B, bool tr) {
    // Factorize
    factorize(A.ptr());

    // Solve
    DM x = densify(B);
    solve(x.ptr(), x.size2());
    return x;
  }

  MX Linsol::solve(const MX& A, const MX& B, bool tr) {
    return (*this)->linsol_solve(A, B, tr);
  }

  void Linsol::solveL(double* x, int nrhs, bool tr) const {
    (*this)->linsol_solveL((*this)->memory(0), x, nrhs, tr);
  }

  void Linsol::factorize(const double* A) const {
    (*this)->linsol_factorize((*this)->memory(0), A);
  }

  void Linsol::solve(double* x, int nrhs, bool tr) const {
    (*this)->linsol_solve((*this)->memory(0), x, nrhs, tr);
  }

  Sparsity Linsol::cholesky_sparsity(bool tr) const {
    return (*this)->linsol_cholesky_sparsity((*this)->memory(0), tr);
  }

  DM Linsol::cholesky(bool tr) const {
    return (*this)->linsol_cholesky((*this)->memory(0), tr);
  }

} // namespace casadi
