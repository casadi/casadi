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
    // Inputs of function being assembled
    vector<XType> ret_in(s_in.size());

    // Handle inputs
    for (int i=0; i<ret_in.size(); ++i) {
      // Locate it in the input scheme
      auto it = std::find(ischeme.begin(), ischeme.end(), s_in[i]);
      if (it!=ischeme.end()) {
        ret_in[i] = in.at(it-ischeme.begin());
        continue;
      }

      // Not found
      casadi_error("Cannot treat input \"" + s_in[i] + "\"");
    }

    // Outputs of function being assembled
    vector<XType> ret_out(s_out.size());

    // List of valid attributes
    const vector<string> all_attr = {"transpose", "triu", "tril", "densify"};

    // Handle outputs
    for (int i=0; i<ret_out.size(); ++i) {
      XType& r = ret_out[i];

      // Output name
      string s = s_out[i];

      // Separarate attributes
      vector<string> attr;
      while (true) {
        // Find the first underscore separator
        size_t pos = s.find('_');
        if (pos>=s.size()) break; // No more underscore
        string a = s.substr(0, pos);

        // Try to match with attributes
        auto it = std::find(all_attr.begin(), all_attr.end(), a);
        if (it==all_attr.end()) break; // No more attribute

        // Save attribute and strip from s
        attr.push_back(*it);
        s = s.substr(pos+1, string::npos);
      }

      // Try to locate in list of outputs
      auto it = std::find(oscheme.begin(), oscheme.end(), s);
      if (it!=oscheme.end()) {
        // Non-differentiated output
        r = out.at(it-oscheme.begin());
      } else {
        // Must be an operator
        size_t pos = s.find('_');
        casadi_assert_message(pos<s.size(), s_out[i] + " is not an output or operator");
        string op = s.substr(0, pos);
        s = s.substr(pos+1);

        // Handle different types of operators
        if (op=="grad") {
          // Output
          pos = s.find('_');
          casadi_assert_message(pos<s.size(), s_out[i] + " is ill-posed");
          string res = s.substr(0, pos);
          int res_i = std::find(oscheme.begin(), oscheme.end(), res) - oscheme.begin();
          casadi_assert_message(res_i<oscheme.size(),
                                "Unrecognized output " + res + " in " + s_out[i]);
          s = s.substr(pos+1);

          // Input
          pos = s.find('_');
          casadi_assert_message(pos<s.size(), s_out[i] + " is ill-posed");
          string arg = s.substr(0, pos);
          int arg_i = std::find(ischeme.begin(), ischeme.end(), arg) - ischeme.begin();
          casadi_assert_message(arg_i<ischeme.size(),
                                "Unrecognized input " + arg + " in " + s_out[i]);
          s = s.substr(pos+1);
          casadi_assert_message(s.empty(), s_out[i] + " is ill-posed");

          // Calculate gradient
          r = gradient(out[res_i], in[arg_i]);

        } else {
          casadi_error("Unknown operator: " + op);
        }
      }

      // Apply attributes (starting from the right-most one)
      for (auto a=attr.rbegin(); a!=attr.rend(); ++a) {
        if (*a=="transpose") {
          r = r = r.T();
        } else if (*a=="triu") {
          r = triu(r);
        } else if (*a=="tril") {
          r = tril(r);
        } else if (*a=="densify") {
          r = densify(r);
        } else {
          casadi_assert(0);
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

