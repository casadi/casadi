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


#include "interpolant_impl.hpp"
#include "../std_vector_tools.hpp"
#include "../mx/mx_node.hpp"
#include <typeinfo>

using namespace std;
namespace casadi {

  bool has_interpolant(const string& name) {
    return Interpolant::hasPlugin(name);
  }

  void load_interpolant(const string& name) {
    Interpolant::loadPlugin(name);
  }

  string doc_interpolant(const string& name) {
    return Interpolant::getPlugin(name).doc;
  }

  Function interpolant(const std::string& name,
                       const std::string& solver,
                       const std::vector<std::vector<double> >& grid,
                       const std::vector<double>& values,
                       const Dict& opts) {
    Function ret;
    ret.assignNode(Interpolant::getPlugin(solver).creator(name, grid, values));
    ret->construct(opts);
    return ret;
  }

  Interpolant::
  Interpolant(const std::string& name,
              const std::vector<std::vector<double> >& grid,
              const std::vector<double>& values)
              : FunctionInternal(name), grid_(grid), values_(values) {
  }

  Interpolant::~Interpolant() {
  }

  Sparsity Interpolant::get_sparsity_in(int i) {
    casadi_assert(i==0);
    return Sparsity::dense(grid_.size());
  }

  Sparsity Interpolant::get_sparsity_out(int i) {
    casadi_assert(i==0);
    return Sparsity::scalar();
  }

  std::string Interpolant::get_name_in(int i) {
    casadi_assert(i==0);
    return "x";
  }

  std::string Interpolant::get_name_out(int i) {
    casadi_assert(i==0);
    return "f";
  }

  void Interpolant::init(const Dict& opts) {
    // Call the base class initializer
    FunctionInternal::init(opts);

  }

  void Interpolant::eval(void* mem, const double** arg, double** res,
                    int* iw, double* w) const {
    casadi_error("not implemented");
  }

  std::map<std::string, Interpolant::Plugin> Interpolant::solvers_;

  const std::string Interpolant::infix_ = "interpolant";

} // namespace casadi
