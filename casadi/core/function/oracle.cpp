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
                                  const std::vector<LinComb>& lincomb,
                                  const Dict& opts) const {
    // Inputs and outputs of function being assembled
    vector<XType> ret_in(s_in.size());
    vector<XType> ret_out(s_out.size());

    // Which input and output expressions have been assigned so far
    vector<bool> assigned_in(s_in.size(), false);
    vector<bool> assigned_out(s_out.size(), false);

    // List of valid attributes
    enum Attributes {
      ATTR_TRANSPOSE,
      ATTR_TRIU,
      ATTR_TRIL,
      ATTR_DENSIFY};

    // Separarate attributes
    vector<vector<Attributes> > attr(s_out.size());
    vector<string> s_out_noatt = s_out;
    for (int i=0; i<s_out_noatt.size(); ++i) {
      // Currently processed string
      string& s = s_out_noatt[i];

      // Loop over attributes
      while (true) {
        // Find the first underscore separator
        size_t pos = s.find('_');
        if (pos>=s.size()) break; // No more underscore
        string a = s.substr(0, pos);

        // Abort if not attribute
        if (a=="transpose") {
          attr[i].push_back(ATTR_TRANSPOSE);
        } else if (a=="triu") {
          attr[i].push_back(ATTR_TRIU);
        } else if (a=="tril") {
          attr[i].push_back(ATTR_TRIL);
        } else if (a=="densify") {
          attr[i].push_back(ATTR_DENSIFY);
        } else {
          // No more attribute
          break;
        }

        // Strip attribute
        s = s.substr(pos+1, string::npos);
      }
    }

    // Assign non-differentiated inputs
    for (int i=0; i<ret_in.size(); ++i) {
      // Locate it in the input scheme
      auto it = std::find(ischeme.begin(), ischeme.end(), s_in[i]);
      if (it!=ischeme.end()) {
        ret_in[i] = in.at(it-ischeme.begin());
        assigned_in[i] = true;
      }
    }

    // Make sure all inputs have been treated
    for (int i=0; i<ret_in.size(); ++i) {
      casadi_assert_message(assigned_in[i], "Cannot treat input \"" + s_in[i] + "\"");
    }

    // Assign non-differentiated outputs
    for (int i=0; i<ret_out.size(); ++i) {
      // Locate it in the input scheme
      auto it = std::find(oscheme.begin(), oscheme.end(), s_out_noatt[i]);
      if (it!=oscheme.end()) {
        ret_out[i] = out.at(it-oscheme.begin());
        assigned_out[i] = true;
      }
    }

    // Make sure all outputs have been treated
    for (int i=0; i<ret_out.size(); ++i) {
      casadi_assert_message(assigned_out[i], "Cannot treat output \"" + s_out[i] + "\"");
    }

    // Apply attributes (starting from the right-most one)
    for (int i=0; i<s_out.size(); ++i) {
      XType& r = ret_out[i];
      for (auto a=attr[i].rbegin(); a!=attr[i].rend(); ++a) {
        switch (*a) {
        case ATTR_TRANSPOSE: r = r.T(); break;
        case ATTR_TRIU: r = triu(r); break;
        case ATTR_TRIL: r = tril(r); break;
        case ATTR_DENSIFY: r = densify(r); break;
        }
      }
    }

    // Create function and return
    return Function(fname, ret_in, ret_out, s_in, s_out, opts);
  }

  Function XProblem::create(const std::string& fname,
                            const std::vector<std::string>& s_in,
                            const std::vector<std::string>& s_out,
                            const std::vector<LinComb>& lincomb,
                            const Dict& opts) const {
    if (is_sx) {
      return sx_p->create(fname, s_in, s_out, lincomb, opts);
    } else {
      return mx_p->create(fname, s_in, s_out, lincomb, opts);
    }
  }

} // namespace casadi

