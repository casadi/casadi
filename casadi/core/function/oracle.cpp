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
#include "external.hpp"

using namespace std;

namespace casadi {

  template<typename XType>
  Function XOracle<XType>::create(const std::string& fname,
                                  const std::vector<std::string>& s_in,
                                  const std::vector<std::string>& s_out,
                                  const std::vector<LinComb>& lincomb,
                                  const Dict& opts) const {
    // Create maps for all inputs and outputs
    map<string, XType> map_in, map_out;

    // Non-differentiated inputs
    for (int i=0; i<in_.size(); ++i) {
      map_in[ischeme_[i]] = in_[i];
    }

    // Non-differentiated outputs
    for (int i=0; i<out_.size(); ++i) {
      map_out[oscheme_[i]] = out_[i];
    }

    // Dual variables
    for (int i=0; i<out_.size(); ++i) {
      string dual_name = "lam_" + oscheme_[i];
      map_in[dual_name] = XType::sym(dual_name, out_[i].sparsity());
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
    Function ret(fname, ret_in, ret_out, s_in, s_out, opts);
    if (ret.has_free()) {
      stringstream ss;
      ss << "Cannot generate " << fname << " as the expressions contain free variables: ";
      if (ret.is_a("sxfunction")) {
        ss << ret.free_sx();
      } else {
        casadi_assert(ret.is_a("mxfunction"));
        ss << ret.free_mx();
      }
      casadi_error(ss.str());
    }
    return ret;
  }

  template<typename XType>
  std::vector<bool>
  XOracle<XType>::nl_var(const std::string& s_in,
                         const std::vector<std::string>& s_out) const {
    // Input arguments
    auto it = find(ischeme_.begin(), ischeme_.end(), s_in);
    casadi_assert(it!=ischeme_.end());
    XType arg = in_.at(it-ischeme_.begin());

    // Output arguments
    vector<XType> res;
    for (auto&& s : s_out) {
      it = find(oscheme_.begin(), oscheme_.end(), s);
      casadi_assert(it!=oscheme_.end());
      res.push_back(out_.at(it-oscheme_.begin()));
    }

    // Extract variables entering nonlinearly
    return XType::nl_var(veccat(res), arg);
  }

  Oracle* Oracle::construct(const std::vector<SX>& in,
                            const std::vector<SX>& out,
                            const std::vector<std::string>& ischeme,
                            const std::vector<std::string>& oscheme) {
    return new XOracle<SX>(in, out, ischeme, oscheme);
  }

  Oracle* Oracle::construct(const std::vector<MX>& in,
                            const std::vector<MX>& out,
                            const std::vector<std::string>& ischeme,
                            const std::vector<std::string>& oscheme) {
    return new XOracle<MX>(in, out, ischeme, oscheme);
  }

  Oracle* Oracle::construct(const Importer& compiler, const std::string& all_io) {
    return new LibOracle<Importer>(compiler, all_io);
  }

  Oracle* Oracle::construct(const std::string& fname, const std::string& all_io) {
    // If fname ends with .c, JIT
    if (fname.size()>2 && fname.compare(fname.size()-2, fname.size(), ".c")==0) {
      Importer compiler(fname, "clang");
      return construct(compiler, all_io);
    } else {
      return new LibOracle<std::string>(fname, all_io);
    }
  }

  std::vector<bool> Oracle::nl_var(const std::string& s_in,
                                   const std::vector<std::string>& s_out) const {
    casadi_error("'nl_var' not defined for " + type_name());
  }

  std::string Oracle::name_in(int i) const {
    casadi_error("'name_in' not defined for " + type_name());
  }

  std::string Oracle::name_out(int i) const {
    casadi_error("'name_out' not defined for " + type_name());
  }

  template<typename XType>
  std::string XOracle<XType>::type_name() const {
    return XType::type_name() + "Oracle";
  }

  Function Oracle::all_io(const std::string& fname, const Dict& opts) const {
    casadi_error("'all_io' not defined for " + type_name());
  }

  template<typename XType>
  Function XOracle<XType>::all_io(const std::string& fname, const Dict& opts) const {
    return Function(fname, in_, out_, ischeme_, oscheme_, opts);
  }

  template<typename LibType>
  LibOracle<LibType>::LibOracle(const LibType& libtype, const std::string& all_io)
    : libtype_(libtype) {
    all_io_ = external(all_io, libtype);
  }

  template<typename LibType>
  const Sparsity& LibOracle<LibType>::sparsity_in(int i) const {
    return all_io_.sparsity_in(i);
  }

  template<typename LibType>
  const Sparsity& LibOracle<LibType>::sparsity_out(int i) const {
    return all_io_.sparsity_out(i);
  }

  template<typename LibType>
  Function LibOracle<LibType>::create(const std::string& fname,
                             const std::vector<std::string>& s_in,
                             const std::vector<std::string>& s_out,
                             const std::vector<LinComb>& lincomb,
                             const Dict& opts) const {
    Function ret = external(fname, libtype_, opts);
    casadi_assert(s_in.size()==ret.n_in());
    casadi_assert(s_out.size()==ret.n_out());
    return ret;
  }

  template<typename LibType>
  std::string LibOracle<LibType>::type_name() const {
    return "LibOracle";
  }

  Oracle* Oracle::expand() const {
    // Get a MX function from all inputs to all outputs
    Function f = all_io();

    // Expand to SX
    if (!f.is_a("sxfunction")) {
      f = f.expand();
    }

    // Get expressions
    vector<SX> arg = f.sx_in();
    vector<SX> res = f(arg);

    // Construct SX oracle
    return construct(arg, res, f.name_in(), f.name_out());
  }

} // namespace casadi

