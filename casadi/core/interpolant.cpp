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
#include "std_vector_tools.hpp"
#include "mx_node.hpp"
#include <typeinfo>

using namespace std;
namespace casadi {

  bool has_interpolant(const string& name) {
    return Interpolant::has_plugin(name);
  }

  void load_interpolant(const string& name) {
    Interpolant::load_plugin(name);
  }

  string doc_interpolant(const string& name) {
    return Interpolant::getPlugin(name).doc;
  }

  Function interpolant(const std::string& name,
                       const std::string& solver,
                       const std::vector<std::vector<double> >& grid,
                       const std::vector<double>& values,
                       const Dict& opts) {

    // Dimension at least 1
    casadi_assert_message(grid.size()>0, "At least one input required");

    // Consistency check, number of elements
    unsigned int nel=1;
    for (auto&& g : grid) {
      casadi_assert_message(g.size()>=2, "Need at least two grid points for every input")
      nel *= g.size();
    }
    casadi_assert_message(nel==values.size(), "Inconsistent number of elements");

    // Grid must be strictly increasing
    for (auto&& g : grid) {
      double last = -inf;
      for (auto&& e : g) {
        casadi_assert_message(!isinf(e) && e>last,
          "Gridpoints must be finite and strictly increasing");
        last = e;
      }
    }

    // Get offset for each input dimension
    vector<int> offset;
    offset.reserve(grid.size()+1);
    offset.push_back(0);
    for (auto&& g : grid) offset.push_back(offset.back()+g.size());

    // Stack input grids
    vector<double> stacked;
    stacked.reserve(offset.back());
    for (auto&& g : grid) stacked.insert(stacked.end(), g.begin(), g.end());


    Function ret;
    ret.own(Interpolant::getPlugin(solver).creator(name, stacked, offset, values));
    ret->construct(opts);
    return ret;
  }

  Interpolant::
  Interpolant(const std::string& name,
              const std::vector<double>& grid,
              const std::vector<int>& offset,
              const std::vector<double>& values)
              : FunctionInternal(name), grid_(grid), offset_(offset), values_(values) {
    // Number of grid points
    ndim_ = offset_.size()-1;
  }

  Interpolant::~Interpolant() {
  }

  Sparsity Interpolant::get_sparsity_in(int i) {
    casadi_assert(i==0);
    return Sparsity::dense(ndim_);
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

  std::map<std::string, Interpolant::Plugin> Interpolant::solvers_;

  const std::string Interpolant::infix_ = "interpolant";

} // namespace casadi
