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


#ifndef CASADI_FACTORY_HPP
#define CASADI_FACTORY_HPP

#include "function.hpp"

/// \cond INTERNAL

namespace casadi {

  // A Jacobian or gradient block
  struct Block {
    std::string ex, arg;
    Block(const std::string& s) {
      size_t pos = s.find(':');
      if (pos<s.size()) {
        this->ex = s.substr(0, pos);
        this->arg = s.substr(pos+1, std::string::npos);
      }
    }
  };

  // A Hessian block
  struct HBlock {
    std::string ex, arg1, arg2;
    HBlock(const std::string& s) {
      size_t pos1 = s.find(':');
      if (pos1<s.size()) {
        size_t pos2 = s.find(':', pos1+1);
        if (pos2<s.size()) {
          this->ex = s.substr(0, pos1);
          this->arg1 = s.substr(pos1+1, pos2-pos1-1);
          this->arg2 = s.substr(pos2+1, std::string::npos);
        }
      }
    }
  };

  // Helper class for generating new functions
  template<typename MatType>
  class Factory {
  public:

    // All auxiliary outputs
    const Function::AuxOut& aux_;

    // All input and output expressions created so far
    std::vector<std::pair<std::string, MatType> > in_, out_;

    // Forward mode directional derivatives
    std::vector<std::string> fwd_in_, fwd_out_;

    // Reverse mode directional derivatives
    std::vector<std::string> adj_in_, adj_out_;

    // Jacobian/gradient blocks
    std::vector<Block> jac_, grad_;

    // Hessian blocks
    std::vector<HBlock> hess_;

    // Constructor
    Factory(const Function::AuxOut& aux) : aux_(aux) {}

    // Add an input expression
    void add_input(const std::string& s, const MatType& e);

    // Add an output expression
    void add_output(const std::string& s, const MatType& e);

    // Request a factory input
    std::string request_input(const std::string& s);

    // Request a factory output
    std::string request_output(const std::string& s);

    // Calculate requested outputs
    void calculate();

    // Retrieve an input
    MatType get_input(const std::string& s);

    // Retrieve an output
    MatType get_output(const std::string& s);

    // Helper function
    static bool has_prefix(const std::string& s);

    // Split prefix
    static std::pair<std::string, std::string> split_prefix(const std::string& s);

    // Check if input exists
    bool has_in(const std::string& s) const;

    // Get an input
    std::string get_in(const std::string& s) const;

    // Check if out exists
    bool has_out(const std::string& s) const;

    // Get an output
    std::string get_out(const std::string& s) const;

    // Get input scheme
    std::vector<std::string> name_in() const;

    // Get output scheme
    std::vector<std::string> name_out() const;

  };

  template<typename MatType>
  void Factory<MatType>::
  add_input(const std::string& s, const MatType& e) {
    casadi_assert_message(!has_in(s), "Duplicate input expression \"" + s + "\"");
    in_.push_back(make_pair(s, e));
  }

  template<typename MatType>
  void Factory<MatType>::
  add_output(const std::string& s, const MatType& e) {
    casadi_assert_message(!has_out(s), "Duplicate output expression \"" + s + "\"");
    out_.push_back(make_pair(s, e));
  }

  template<typename MatType>
  std::string Factory<MatType>::
  request_input(const std::string& s) {
    using namespace std;

    // Quick return if already available
    if (has_in(s)) return get_in(s);

    // Get prefix
    casadi_assert_message(has_prefix(s), "Cannot process \"" + s + "\" as input."
                                         " Available: " + join(name_in()) + ".");
    pair<string, string> ss = split_prefix(s);

    if (ss.first=="fwd") {
      // Forward mode directional derivative
      fwd_in_.push_back(get_in(ss.second));
      return "fwd_" + fwd_in_.back();
    } else if (ss.first=="adj") {
      // Reverse mode directional derivative
      adj_in_.push_back(get_out(ss.second));
      return "adj_" + adj_in_.back();
    }
  }

  template<typename MatType>
  std::string Factory<MatType>::
  request_output(const std::string& s) {
    using namespace std;

    // Quick return if already available
    if (has_out(s)) return get_out(s);

    // Get prefix
    casadi_assert_message(has_prefix(s), "Cannot process \"" + s + "\" as output."
                                         " Available: " + join(name_out()) + ".");
    pair<string, string> ss = split_prefix(s);

    if (ss.first=="fwd") {
      // Forward mode directional derivative
      ss.second = get_out(ss.second);
      fwd_out_.push_back(ss.second);
    } else if (ss.first=="adj") {
      // Reverse mode directional derivative
      ss.second = get_in(ss.second);
      adj_out_.push_back(ss.second);
    } else if (ss.first=="jac") {
      jac_.push_back(ss.second);
      jac_.back().ex = get_out(jac_.back().ex);
      jac_.back().arg = get_in(jac_.back().arg);
      ss.second = jac_.back().ex + "_" + jac_.back().arg;
    } else if (ss.first=="grad") {
      grad_.push_back(ss.second);
      grad_.back().ex = get_out(grad_.back().ex);
      grad_.back().arg = get_in(grad_.back().arg);
      ss.second = grad_.back().ex + "_" + grad_.back().arg;
    } else if (ss.first=="hess") {
      hess_.push_back(ss.second);
      hess_.back().ex = get_out(hess_.back().ex);
      hess_.back().arg1 = get_in(hess_.back().arg1);
      hess_.back().arg2 = get_in(hess_.back().arg2);
      ss.second = hess_.back().ex + "_" + hess_.back().arg1 + "_" + hess_.back().arg2;
    } else {
      // Assume attribute
      return request_output(ss.second);
    }

    // Replace colons with underscore
    return ss.first + "_" + ss.second;
  }

  template<typename MatType>
  void Factory<MatType>::calculate() {
    using namespace std;

    // Dual variables
    for (auto&& e : out_) {
      in_.push_back(make_pair("lam:" + e.first, MatType::sym("lam_" + e.first, e.second.sparsity())));
    }

    // Forward mode directional derivatives
    if (!fwd_out_.empty()) {
      casadi_assert(!fwd_in_.empty());

      vector<MatType> arg, res;
      vector<vector<MatType>> seed(1), sens(1);
      // Inputs and forward mode seeds
      for (const string& s : fwd_in_) {
        arg.push_back(get_input(s));
        seed[0].push_back(MatType::sym("fwd_" + s, arg.back().sparsity()));
        in_.push_back(make_pair("fwd:" + s, seed[0].back()));
      }
      // Outputs
      for (const string& s : fwd_out_) {
        res.push_back(get_output(s));
      }
      // Calculate directional derivatives
      Dict opts = {{"always_inline", true}};
      sens = forward(res, arg, seed, opts);

      // Get directional derivatives
      for (int i=0; i<fwd_out_.size(); ++i) {
        out_.push_back(make_pair("fwd:" + fwd_out_[i],
				 project(sens[0].at(i), res.at(i).sparsity())));
      }
    }

    // Reverse mode directional derivatives
    if (!adj_out_.empty()) {
      casadi_assert(!adj_in_.empty());
      vector<MatType> arg, res;
      vector<vector<MatType>> seed(1), sens(1);
      // Inputs
      for (const string& s : adj_out_) {
        arg.push_back(get_input(s));
      }
      // Outputs and reverse mode seeds
      for (const string& s : adj_in_) {
        res.push_back(get_output(s));
        seed[0].push_back(MatType::sym("adj_" + s, res.back().sparsity()));
        in_.push_back(make_pair("adj:" + s, seed[0].back()));
      }
      // Calculate directional derivatives
      Dict opts = {{"always_inline", true}};
      sens = reverse(res, arg, seed, opts);

      // Get directional derivatives
      for (int i=0; i<adj_out_.size(); ++i) {
        out_.push_back(make_pair("adj:" + adj_out_[i],
				 project(sens[0].at(i), arg.at(i).sparsity())));
      }
    }

    // Add linear combinations
    for (auto i : aux_) {
      MatType lc = 0;
      for (auto j : i.second) {
        lc += dot(get_input("lam:" + j), get_output(j));
      }
      out_.push_back(make_pair(i.first, lc));
    }

    // Jacobian blocks
    for (auto &&b : jac_) {
      const MatType& ex = get_output(b.ex);
      const MatType& arg = get_input(b.arg);
      out_.push_back(make_pair("jac:" + b.ex + ":" + b.arg,
			       MatType::jacobian(ex, arg)));
    }

    // Gradient blocks
    for (auto &&b : grad_) {
      const MatType& ex = get_output(b.ex);
      const MatType& arg = get_input(b.arg);
      out_.push_back(make_pair("grad:" + b.ex + ":" + b.arg,
			       project(gradient(ex, arg), arg.sparsity())));
    }

    // Hessian blocks
    for (auto &&b : hess_) {
      const MatType& ex = get_output(b.ex);
      casadi_assert_message(b.arg1==b.arg2, "Mixed Hessian terms not supported");
      const MatType& arg1 = get_input(b.arg1);
      //const MatType& arg2 = get_input(b.arg2);
      out_.push_back(make_pair("hess:" + b.ex + ":" + b.arg1 + ":" + b.arg2,
			       triu(hessian(ex, arg1))));
    }
  }

  template<typename MatType>
  MatType Factory<MatType>::get_input(const std::string& s) {
    for (int i=0; i<in_.size(); ++i) {
      if (s==in_[i].first || s==str(i)) return in_[i].second;
    }
    casadi_error("Cannot retrieve \"" + s + "\"");
  }

  template<typename MatType>
  MatType Factory<MatType>::get_output(const std::string& s) {
    using namespace std;

    // Quick return if output
    for (int i=0; i<out_.size(); ++i) {
      if (s==out_[i].first || s==str(i)) return out_[i].second;
    }
    
    // Assume attribute
    casadi_assert_message(has_prefix(s), "Cannot process \"" + s + "\"");
    pair<string, string> ss = split_prefix(s);
    string a = ss.first;
    MatType r = get_output(ss.second);

    // Process attributes
    if (a=="transpose") {
      return r.T();
    } else if (a=="triu") {
      return triu(r);
    } else if (a=="tril") {
      return tril(r);
    } else if (a=="densify") {
      return densify(r);
    } else if (a=="sym") {
      return triu2symm(r);
    } else if (a=="withdiag") {
      return project(r, r.sparsity() + Sparsity::diag(r.size1()));
    } else {
      casadi_error("Cannot process attribute \"" + a + "\"");
      return MatType();
    }
  }

  template<typename MatType>
  bool Factory<MatType>::
  has_prefix(const std::string& s) {
    return s.find(':') < s.size();
  }

  template<typename MatType>
  std::pair<std::string, std::string> Factory<MatType>::
  split_prefix(const std::string& s) {
    // Get prefix
    casadi_assert(!s.empty());
    size_t pos = s.find(':');
    casadi_assert_message(pos<s.size(), "Cannot process \"" + s + "\"");
    return make_pair(s.substr(0, pos), s.substr(pos+1, std::string::npos));
  }

  template<typename MatType>
  std::vector<std::string> Factory<MatType>::name_in() const {
    std::vector<std::string> ret;
    for (auto i : in_) {
      ret.push_back(i.first);
    }
    return ret;
  }

  template<typename MatType>
  std::vector<std::string> Factory<MatType>::name_out() const {
    std::vector<std::string> ret;
    for (auto i : out_) {
      ret.push_back(i.first);
    }
    return ret;
  }

  template<typename MatType>
  bool Factory<MatType>::has_in(const std::string& s) const {
    if (s.empty()) return false;
    for (int i=0; i<in_.size(); ++i) {
      if (s==in_[i].first || s==str(i)) return true;
    }
  }

  template<typename MatType>
  std::string Factory<MatType>::get_in(const std::string& s) const {
    if (!has_in(s)) {
      casadi_error("Cannot process \"" + s + "\" as input. Available: " + str(name_in()) + ".");
    }
    if (isdigit(s.at(0))) {
      for (int i=0; i<in_.size(); ++i) {
	if (s==str(i)) return in_[i].first;
      }
    }
    return s;
  }

  template<typename MatType>
  bool Factory<MatType>::has_out(const std::string& s) const {
    if (s.empty()) return false;
    // Standard output
    for (int i=0; i<out_.size(); ++i) {
      if (s==out_[i].first || s==str(i)) return true;
    }
    // Auxiliary output?
    return aux_.find(s)!=aux_.end();
  }

  template<typename MatType>
  std::string Factory<MatType>::get_out(const std::string& s) const {
    if (!has_out(s)) {
      std::vector<std::string> all_out = name_out();
      for (auto i : aux_) all_out.push_back(i.first);
      casadi_error("Cannot process \"" + s + "\" as output. Available: " + str(all_out) + ".");
    }
    if (isdigit(s.at(0))) {
      for (int i=0; i<out_.size(); ++i) {
	if (s==str(i)) return out_[i].first;
      }
    }
    return s;
  }

} // namespace casadi
/// \endcond

#endif // CASADI_FACTORY_HPP
