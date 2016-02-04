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

  XProblem::XProblem(const SXProblem& d) : p(new SXProblem(d)) {
  }

  XProblem::XProblem(const MXProblem& d) : p(new MXProblem(d)) {
  }

  XProblem::~XProblem() {
    delete p;
  }

  XProblem::XProblem(const XProblem& d) {
    p = d.p->clone();
  }

  XProblem& XProblem::operator=(const XProblem& d) {
    if (&d!=this) {
      delete p;

      // Assign
      p = d.p->clone();
    }
    return *this;
  }

  const Sparsity& XProblem::sparsity_in(int i) const {
    return p->sparsity_in(i);
  }

  const Sparsity& XProblem::sparsity_out(int i) const {
    return p->sparsity_out(i);
  }

  template<typename XType>
  Function Problem<XType>::create(const std::string& fname,
                                  const std::vector<std::string>& s_in,
                                  const std::vector<std::string>& s_out,
                                  const std::vector<LinComb>& lincomb,
                                  const Dict& opts) const {
    // Create maps for all inputs and outputs
    map<string, XType> map_in, map_out;

    // Non-differentiated inputs
    for (int i=0; i<in.size(); ++i) {
      map_in[ischeme[i]] = in[i];
    }

    // Non-differentiated outputs
    for (int i=0; i<out.size(); ++i) {
      map_out[oscheme[i]] = out[i];
    }

    // Dual variables
    for (int i=0; i<out.size(); ++i) {
      string dual_name = "lam_" + oscheme[i];
      map_in[dual_name] = XType::sym(dual_name, out[i].sparsity());
    }

    // Add linear combinations
    for (auto i : lincomb) {
      XType lc = 0;
      for (auto j : i.second) {
        lc += dot(map_in.at("lam_" + j), map_out.at(j));
      }
      map_out[i.first] = lc;
    }

    // Inputs of function being assembled
    vector<XType> ret_in(s_in.size());

    // Handle inputs
    for (int i=0; i<ret_in.size(); ++i) {
      ret_in[i] = map_in.at(s_in[i]);
    }

    // Outputs of function being assembled
    vector<XType> ret_out(s_out.size());

    // List of valid attributes
    const vector<string> all_attr = {"transpose", "triu", "tril", "densify", "sym", "withdiag"};

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
      auto it = map_out.find(s);
      if (it!=map_out.end()) {
        // Non-differentiated output
        r = it->second;
      } else {
        // Must be an operator
        size_t pos = s.find('_');
        casadi_assert_message(pos<s.size(), s_out[i] + " is not an output or operator");
        string op = s.substr(0, pos);
        s = s.substr(pos+1);

        // Handle different types of operators
        if (op=="grad" || op=="jac" || op=="hess") {
          // Output
          pos = s.find('_');
          casadi_assert_message(pos<s.size(), s_out[i] + " is ill-posed");
          string res = s.substr(0, pos);
          auto res_i = map_out.find(res);
          casadi_assert_message(res_i!=map_out.end(),
                                "Unrecognized output " + res + " in " + s_out[i]);
          s = s.substr(pos+1);

          // Input
          string arg;
          if (op=="hess") {
            pos = s.find('_');
            casadi_assert_message(pos<s.size(), s_out[i] + " is ill-posed");
            arg = s.substr(0, pos);
            s = s.substr(pos+1);
            casadi_assert_message(s==arg, "Mixed Hessian terms not supported");
          } else {
            arg = s;
          }
          auto arg_i = map_in.find(arg);
          casadi_assert_message(arg_i!=map_in.end(),
                                "Unrecognized input " + arg + " in " + s_out[i]);

          // Calculate gradient or Jacobian
          if (op=="jac") {
            r = jacobian(res_i->second, arg_i->second);
          } else if (op=="grad") {
            r = project(gradient(res_i->second, arg_i->second), arg_i->second.sparsity());
          } else {
            casadi_assert(op=="hess");
            r = triu(hessian(res_i->second, arg_i->second));
          }
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
        } else if (*a=="sym") {
          r = triu2symm(r);
        } else if (*a=="withdiag") {
          r = project(r, r.sparsity() + Sparsity::diag(r.size1()));
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
    return p->create(fname, s_in, s_out, lincomb, opts);
  }

  template<typename XType>
  std::vector<bool> Problem<XType>::nl_var(const std::string& s_in,
                                           const std::vector<std::string>& s_out) {
    // Input arguments
    auto it = find(ischeme.begin(), ischeme.end(), s_in);
    casadi_assert(it!=ischeme.end());
    XType arg = in.at(it-ischeme.begin());

    // Output arguments
    vector<XType> res;
    for (auto&& s : s_out) {
      it = find(oscheme.begin(), oscheme.end(), s);
      casadi_assert(it!=oscheme.end());
      res.push_back(out.at(it-oscheme.begin()));
    }

    // Extract variables entering nonlinearly
    return XType::nl_var(veccat(res), arg);
  }

  std::vector<bool> XProblem::nl_var(const std::string& s_in,
                                     const std::vector<std::string>& s_out) {
    return p->nl_var(s_in, s_out);
  }

  Oracle* Oracle::create(const std::vector<SX>& in,
                         const std::vector<SX>& out,
                         const std::vector<std::string>& ischeme,
                         const std::vector<std::string>& oscheme) {
    return new Problem<SX>(in, out, ischeme, oscheme);
  }

  Oracle* Oracle::create(const std::vector<MX>& in,
                         const std::vector<MX>& out,
                         const std::vector<std::string>& ischeme,
                         const std::vector<std::string>& oscheme) {
    return new Problem<MX>(in, out, ischeme, oscheme);
  }

} // namespace casadi

