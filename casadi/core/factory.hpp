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
    std::map<std::string, MatType> in_, out_;
    std::map<std::string, bool> is_diff_in_, is_diff_out_;

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
    void add_input(const std::string& s, const MatType& e, bool is_diff);

    // Add an output expression
    void add_output(const std::string& s, const MatType& e, bool is_diff);

    // Request a factory input
    std::string request_input(const std::string& s);

    // Request a factory output
    std::string request_output(const std::string& s);

    // Calculate forward mode directional derivatives
    void calculate_fwd(const Dict& opts);

    // Calculate reverse mode directional derivatives
    void calculate_adj(const Dict& opts);

    // Calculate Jacobian blocks
    void calculate_jac(const Dict& opts);

    // Calculate gradient blocks
    void calculate_grad(const Dict& opts);

    // Calculate Hessian blocks
    void calculate_hess(const Dict& opts);

    // Calculate requested outputs
    void calculate(const Dict& opts = Dict());

    // Retrieve an input
    MatType get_input(const std::string& s);

    // Retrieve an output
    MatType get_output(const std::string& s);

    // Helper function
    static bool has_prefix(const std::string& s);

    // Split prefix
    static std::pair<std::string, std::string> split_prefix(const std::string& s);

    // Check if input exists
    bool has_in(const std::string& s) const { return in_.find(s)!=in_.end();}

    // Check if out exists
    bool has_out(const std::string& s) const {
      // Standard output
      if (out_.find(s)!=out_.end()) return true;
      // Auxiliary output?
      return aux_.find(s)!=aux_.end();
    }

    // Get input scheme
    std::vector<std::string> name_in() const;

    // Get output scheme
    std::vector<std::string> name_out() const;
  };

  template<typename MatType>
  void Factory<MatType>::
  add_input(const std::string& s, const MatType& e, bool is_diff) {
    auto it = in_.insert(make_pair(s, e));
    casadi_assert(it.second, "Duplicate input expression \"" + s + "\"");
    is_diff_in_.insert(make_pair(s, is_diff));
  }

  template<typename MatType>
  void Factory<MatType>::
  add_output(const std::string& s, const MatType& e, bool is_diff) {
    auto it = out_.insert(make_pair(s, e));
    casadi_assert(it.second, "Duplicate output expression \"" + s + "\"");
    is_diff_out_.insert(make_pair(s, is_diff));
  }

  template<typename MatType>
  std::string Factory<MatType>::
  request_input(const std::string& s) {
    using namespace std;

    // Quick return if already available
    if (has_in(s)) return s;

    // Get prefix
    casadi_assert(has_prefix(s), "Cannot process \"" + s + "\" as input."
                                         " Available: " + join(name_in()) + ".");
    pair<string, string> ss = split_prefix(s);

    if (ss.first=="fwd") {
      // Forward mode directional derivative
      casadi_assert(has_in(ss.second), "Cannot process \"" + ss.second + "\""
                                               " (from \"" + s + "\") as input."
                                               " Available: " + join(name_in()) + ".");
      fwd_in_.push_back(ss.second);
    } else if (ss.first=="adj") {
      // Reverse mode directional derivative
      casadi_assert(has_out(ss.second), "Cannot process \"" + ss.second + "\""
                                                " (from \"" + s + "\") as output."
                                                " Available: " + join(name_out()) + ".");
      adj_in_.push_back(ss.second);
    }

    // Replace colons with underscore
    string ret = s;
    replace(ret.begin(), ret.end(), ':', '_');
    return ret;
  }

  template<typename MatType>
  std::string Factory<MatType>::
  request_output(const std::string& s) {
    using namespace std;

    // Quick return if already available
    if (has_out(s)) return s;

    // Get prefix
    casadi_assert(has_prefix(s), "Cannot process \"" + s + "\" as output."
                                         " Available: " + join(name_out()) + ".");
    pair<string, string> ss = split_prefix(s);

    if (ss.first=="fwd") {
      // Forward mode directional derivative
      casadi_assert(has_out(ss.second), "Cannot process \"" + ss.second + "\""
                                                " (from \"" + s + "\") as output."
                                                " Available: " + join(name_out()) + ".");
      fwd_out_.push_back(ss.second);
    } else if (ss.first=="adj") {
      // Reverse mode directional derivative
      casadi_assert(has_in(ss.second),
        "Cannot process \"" + ss.second + "\" (from \"" + s + "\") as input. "
        "Available: " + join(name_in()) + ".");
      adj_out_.push_back(ss.second);
    } else if (ss.first=="jac") {
      jac_.push_back(ss.second);
      casadi_assert(has_out(jac_.back().ex),
        "Cannot process \"" + jac_.back().ex + "\" (from \"" + s + "\") as output. "
        "Available: " + join(name_out()) + ".");
      casadi_assert(has_in(jac_.back().arg),
        "Cannot process \"" + jac_.back().arg + "\" (from \"" + s + "\") as input. "
        "Available: " + join(name_in()) + ".");
    } else if (ss.first=="grad") {
      grad_.push_back(ss.second);
      casadi_assert(has_out(grad_.back().ex),
        "Cannot process \"" + grad_.back().ex + "\" (from \"" + s + "\") as output. "
        "Available: " + join(name_out()) + ".");
      casadi_assert(has_in(grad_.back().arg),
        "Cannot process \"" + grad_.back().arg + "\" (from \"" + s + "\") as input. "
        "Available: " + join(name_in()) + ".");
    } else if (ss.first=="hess") {
      hess_.push_back(ss.second);
      casadi_assert(has_out(hess_.back().ex),
        "Cannot process \"" + hess_.back().ex + "\" (from \"" + s + "\") as output. "
        "Available: " + join(name_out()) + ".");
      casadi_assert(has_in(hess_.back().arg1),
        "Cannot process \"" + hess_.back().arg1 + "\" (from \"" + s + "\") as input. "
        "Available: " + join(name_in()) + ".");
      casadi_assert(has_in(hess_.back().arg2),
        "Cannot process \"" + hess_.back().arg2 + "\" (from \"" + s + "\") as input. "
        "Available: " + join(name_in()) + ".");
    } else {
      // Assume attribute
      request_output(ss.second);
    }

    // Replace colons with underscore
    string ret = s;
    replace(ret.begin(), ret.end(), ':', '_');
    return ret;
  }

  template<typename MatType>
  void Factory<MatType>::calculate_fwd(const Dict& opts) {
    if (fwd_out_.empty()) return;
    casadi_assert_dev(!fwd_in_.empty());

    std::vector<MatType> arg, res;
    std::vector<std::vector<MatType>> seed(1), sens(1);
    // Inputs and forward mode seeds
    for (const std::string& s : fwd_in_) {
      arg.push_back(in_[s]);
      Sparsity sp = is_diff_in_[s] ? arg.back().sparsity() : Sparsity(arg.back().size());
      seed[0].push_back(MatType::sym("fwd_" + s, sp));
      in_["fwd:" + s] = seed[0].back();
    }
    // Outputs
    for (const std::string& s : fwd_out_) {
      res.push_back(out_[s]);
    }
    // Calculate directional derivatives
    Dict local_opts = opts;
    local_opts["always_inline"] = true;
    sens = forward(res, arg, seed, local_opts);

    // Get directional derivatives
    for (casadi_int i=0; i<fwd_out_.size(); ++i) {
      std::string s = fwd_out_[i];
      Sparsity sp = is_diff_out_[s] ? res.at(i).sparsity() : Sparsity(res.at(i).size());
      out_["fwd:" + s] = project(sens[0].at(i), sp);
      is_diff_out_["fwd:" + s] = is_diff_out_[s];
    }
  }

  template<typename MatType>
  void Factory<MatType>::calculate_adj(const Dict& opts) {
    if (adj_out_.empty()) return;
    casadi_assert_dev(!adj_in_.empty());
    std::vector<MatType> arg, res;
    std::vector<std::vector<MatType>> seed(1), sens(1);
    // Inputs
    for (const std::string& s : adj_out_) {
      arg.push_back(in_[s]);
    }
    // Outputs and reverse mode seeds
    for (const std::string& s : adj_in_) {
      res.push_back(out_[s]);
      Sparsity sp = is_diff_out_[s] ? res.back().sparsity() : Sparsity(res.back().size());
      seed[0].push_back(MatType::sym("adj_" + s, sp));
      in_["adj:" + s] = seed[0].back();
    }
    // Calculate directional derivatives
    Dict local_opts;
    local_opts["always_inline"] = true;
    sens = reverse(res, arg, seed, local_opts);

    // Get directional derivatives
    for (casadi_int i=0; i<adj_out_.size(); ++i) {
      std::string s = adj_out_[i];
      Sparsity sp = is_diff_in_[s] ? arg.at(i).sparsity() : Sparsity(arg.at(i).size());
      out_["adj:" + s] = project(sens[0].at(i), sp);
      is_diff_in_["adj:" + s] = is_diff_in_[s];
    }
  }

  template<typename MatType>
  void Factory<MatType>::calculate_jac(const Dict& opts) {
    // Calculate blocks for all non-differentiable inputs and outputs
    for (auto &&b : jac_) {
      std::string s = "jac:" + b.ex + ":" + b.arg;
      if (!is_diff_out_.at(b.ex) || !is_diff_in_.at(b.arg)) {
        out_[s] = MatType(out_.at(b.ex).numel(), in_.at(b.arg).numel());
        is_diff_out_[s] = false;
      } else {
        is_diff_out_[s] = true;
      }
    }
    // Calculate regular blocks
    for (auto &&b : jac_) {
      // Get block name, skip if already calculated
      std::string s = "jac:" + b.ex + ":" + b.arg;
      if (out_.find(s) != out_.end()) continue;
      // Find other blocks with the same input, but different (not yet calculated) outputs
      std::vector<MatType> ex;
      std::vector<std::string> all_ex;
      for (auto &&b1 : jac_) {
        // Check if same input
        if (b1.arg != b.arg) continue;
        // Check if already calculated
        std::string s1 = "jac:" + b1.ex + ":" + b1.arg;
        if (out_.find(s1) != out_.end()) continue;
        // Collect expressions
        all_ex.push_back(b1.ex);
        ex.push_back(out_.at(b1.ex));
      }
      // Now find other blocks with *all* the same outputs, but different inputs
      std::vector<MatType> arg{in_.at(b.arg)};
      std::vector<std::string> all_arg{b.arg};
      for (auto &&b1 : jac_) {
        // Candidate b1.arg: Check if already added
        bool skip = false;
        for (const std::string& a : all_arg) {
          if (a == b1.arg) {
            skip = true;
            break;
          }
        }
        if (skip) continue;
        // Check if all blocks corresponding to the same input are needed
        for (const std::string& e : all_ex) {
          std::string s1 = "jac:" + e + ":" + b1.arg;
          if (is_diff_out_.find(s1) == is_diff_out_.end() || out_.find(s1) != out_.end()) {
            // Block is not requested or has already been calculated
            skip = true;
            break;
          }
        }
        if (skip) continue;
        // Keep candidate
        arg.push_back(in_.at(b1.arg));
        all_arg.push_back(b1.arg);
      }
      try {
        // Calculate Jacobian block(s)
        if (ex.size() == 1 && arg.size() == 1) {
          out_[s] = MatType::jacobian(ex[0], arg[0], opts);
        } else {
          // Calculate Jacobian of all outputs with respect to all inputs
          MatType J = MatType::jacobian(vertcat(ex), vertcat(arg), opts);
          // Split Jacobian into blocks
          std::vector<std::vector<MatType>> J_all = blocksplit(J, offset(ex), offset(arg));
          // Save blocks
          for (size_t i = 0; i < all_ex.size(); ++i) {
            for (size_t j = 0; j < all_arg.size(); ++j) {
              out_["jac:" + all_ex[i] + ":" + all_arg[j]] = J_all.at(i).at(j);
            }
          }
        }
      } catch (std::exception& e) {
        std::stringstream ss;
        ss << "Calculating Jacobian of " << all_ex << " w.r.t. " << all_arg << ": " << e.what();
        casadi_error(ss.str());
      }
    }
  }

  template<typename MatType>
  void Factory<MatType>::calculate_grad(const Dict& opts) {
    for (auto &&b : grad_) {
      const MatType& ex = out_.at(b.ex);
      const MatType& arg = in_.at(b.arg);
      if (is_diff_out_.at(b.ex) && is_diff_in_.at(b.arg)) {
        out_["grad:" + b.ex + ":" + b.arg] = project(gradient(ex, arg, opts), arg.sparsity());
        is_diff_out_["grad:" + b.ex + ":" + b.arg] = true;
      } else {
        casadi_assert(ex.is_scalar(), "Can only take gradient of scalar expression.");
        out_["grad:" + b.ex + ":" + b.arg] = MatType(1, arg.numel());
        is_diff_out_["grad:" + b.ex + ":" + b.arg] = false;
      }
    }
  }

  template<typename MatType>
  void Factory<MatType>::calculate_hess(const Dict& opts) {
    for (auto &&b : hess_) {
      const MatType& ex = out_.at(b.ex);
      casadi_assert(b.arg1==b.arg2, "Mixed Hessian terms not supported");
      const MatType& arg1 = in_.at(b.arg1);
      //const MatType& arg2 = in_.at(b.arg2);
      if (is_diff_out_.at(b.ex) && is_diff_in_.at(b.arg1)) {
        out_["hess:" + b.ex + ":" + b.arg1 + ":" + b.arg2] = triu(hessian(ex, arg1, opts));
        is_diff_out_["hess:" + b.ex + ":" + b.arg1 + ":" + b.arg2] = true;
      } else {
        casadi_assert(ex.is_scalar(), "Can only take Hessian of scalar expression.");
        out_["hess:" + b.ex + ":" + b.arg1 + ":" + b.arg2] = MatType(arg1.numel(), arg1.numel());
        is_diff_out_["hess:" + b.ex + ":" + b.arg1 + ":" + b.arg2] = false;
      }
    }
  }

  template<typename MatType>
  void Factory<MatType>::calculate(const Dict& opts) {
    // Dual variables
    for (auto&& e : out_) {
      Sparsity sp = is_diff_out_[e.first] ? e.second.sparsity() : Sparsity(e.second.size());
      in_["lam:" + e.first] = MatType::sym("lam_" + e.first, sp);
    }

    // Forward mode directional derivatives
    try {
      calculate_fwd(opts);
    } catch (std::exception& e) {
      casadi_error("Forward mode AD failed:\n" + str(e.what()));
    }

    // Reverse mode directional derivatives
    try {
      calculate_adj(opts);
    } catch (std::exception& e) {
      casadi_error("Reverse mode AD failed:\n" + str(e.what()));
    }

    // Add linear combinations
    for (auto i : aux_) {
      MatType lc = 0;
      for (auto j : i.second) {
        lc += dot(in_.at("lam:" + j), out_.at(j));
      }
      out_[i.first] = lc;
      is_diff_out_[i.first] = true;
    }

    // Jacobian blocks
    try {
      calculate_jac(opts);
    } catch (std::exception& e) {
      casadi_error("Jacobian generation failed:\n" + str(e.what()));
    }

    // Gradient blocks
    try {
      calculate_grad(opts);
    } catch (std::exception& e) {
      casadi_error("Gradient generation failed:\n" + str(e.what()));
    }

    // Hessian blocks
    try {
      calculate_hess(opts);
    } catch (std::exception& e) {
      casadi_error("Hessian generation failed:\n" + str(e.what()));
    }
  }

  template<typename MatType>
  MatType Factory<MatType>::get_input(const std::string& s) {
    auto it = in_.find(s);
    casadi_assert(it!=in_.end(), "Cannot retrieve \"" + s + "\"");
    return it->second;
  }

  template<typename MatType>
  MatType Factory<MatType>::get_output(const std::string& s) {
    using namespace std;

    // Quick return if output
    auto it = out_.find(s);
    if (it!=out_.end()) return it->second;

    // Assume attribute
    casadi_assert(has_prefix(s), "Cannot process \"" + s + "\"");
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
    casadi_assert_dev(!s.empty());
    size_t pos = s.find(':');
    casadi_assert(pos<s.size(), "Cannot process \"" + s + "\"");
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

} // namespace casadi
/// \endcond

#endif // CASADI_FACTORY_HPP
