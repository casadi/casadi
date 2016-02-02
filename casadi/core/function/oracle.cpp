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


#include "oracle.hpp"

using namespace std;

namespace casadi {

  XProblem::XProblem(const SXProblem& d) : sx_p(new SXProblem(d)), is_sx(true) {
  }

  XProblem::XProblem(const MXProblem& d) : mx_p(new MXProblem(d)), is_sx(false) {
  }

  XProblem::~XProblem() {
    if (is_sx) {
      delete sx_p;
    } else {
      delete mx_p;
    }
  }

  XProblem::XProblem(const XProblem& d) : is_sx(d.is_sx) {
    if (d.is_sx) {
      sx_p = new SXProblem(*d.sx_p);
    } else {
      mx_p = new MXProblem(*d.mx_p);
    }
  }

  XProblem& XProblem::operator=(const XProblem& d) {
    if (&d!=this) {
      // Delete the previous object
      if (is_sx) {
        delete sx_p;
      } else {
        delete mx_p;
      }
      // Assign
      is_sx = d.is_sx;
      if (is_sx) {
        sx_p = new SXProblem(*d.sx_p);
      } else {
        mx_p = new MXProblem(*d.mx_p);
      }
    }
    return *this;
  }

  XProblem::operator const SXProblem&() const {
    casadi_assert(is_sx);
    return *sx_p;
  }

  XProblem::operator const MXProblem&() const {
    casadi_assert(!is_sx);
    return *mx_p;
  }

  const Sparsity& XProblem::sparsity_in(int i) const {
    return is_sx ? sx_p->sparsity_in(i) : mx_p->sparsity_in(i);
  }

  const Sparsity& XProblem::sparsity_out(int i) const {
    return is_sx ? sx_p->sparsity_out(i) : mx_p->sparsity_out(i);
  }

  template<typename XType>
  Function Problem<XType>::create(const std::string& fname,
                                  const std::vector<std::string>& s_in,
                                  const std::vector<std::string>& s_out,
                                  const Dict& opts) const {
    // Assemble inputs
    std::vector<XType> x_in(s_in.size());
    for (int i=0; i<x_in.size(); ++i) {
      x_in[i] = in.at(std::find(ischeme.begin(), ischeme.end(), s_in[i])-ischeme.begin());
    }

    // Assemble outputs
    std::vector<XType> x_out(s_out.size());
    for (int i=0; i<x_out.size(); ++i) {
      x_out[i] = out.at(std::find(oscheme.begin(), oscheme.end(), s_out[i])-oscheme.begin());
    }

    // Create function and return
    return Function(fname, x_in, x_out, s_in, s_out, opts);
  }

  Function XProblem::create(const std::string& fname,
                            const std::vector<std::string>& s_in,
                            const std::vector<std::string>& s_out,
                            const Dict& opts) const {
    if (is_sx) {
      return sx_p->create(fname, s_in, s_out, opts);
    } else {
      return mx_p->create(fname, s_in, s_out, opts);
    }
  }

} // namespace casadi

