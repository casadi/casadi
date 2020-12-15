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
    size_t f, x;
  };

  // A Hessian block
  struct HBlock {
    size_t f, x1, x2;
    std::string s;
    bool calculated;
  };

  // Helper class for generating new functions
  template<typename MatType>
  class Factory {
  public:

    // All input and output expressions
    std::vector<MatType> in_, out_;

    // Names of inputs and outputs
    std::vector<std::string> iname_, oname_;

    // All input and output expressions created so far
    std::map<std::string, size_t> imap_, omap_;
    std::map<std::string, bool> is_diff_imap_, is_diff_omap_;

    // Forward mode directional derivatives
    std::vector<std::string> fwd_imap_, fwd_omap_;

    // Reverse mode directional derivatives
    std::vector<std::string> adj_imap_, adj_omap_;

    // Jacobian/gradient blocks
    std::vector<Block> jac_, grad_;

    // Hessian blocks
    std::vector<HBlock> hess_;

    // Read a Jacobian or gradient block
    Block block(const std::string& s) const;

    // Read a Hessian block
    HBlock hblock(const std::string& s) const;

    // Add an input expression
    void add_input(const std::string& s, const MatType& e, bool is_diff);

    // Add an output expression
    void add_output(const std::string& s, const MatType& e, bool is_diff);

    // Add the dual variables
    void add_dual(const Function::AuxOut& aux);

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

    // Calculate Hessian blocks, one expression
    void calculate_hess(const Dict& opts, size_t f);

    // Calculate Hessian blocks
    void calculate_hess(const Dict& opts);

    // Calculate requested outputs
    void calculate(const Dict& opts = Dict());

    // Get input index from string
    size_t imap(const std::string& s) const;

    // Get output index from string
    size_t omap(const std::string& s) const;

    // Retrieve an input
    MatType get_input(const std::string& s);

    // Retrieve an output
    MatType get_output(const std::string& s);

    // Helper function
    static bool has_prefix(const std::string& s);

    // Split prefix
    static std::pair<std::string, std::string> split_prefix(const std::string& s);

    // Check if input exists
    bool has_in(const std::string& s) const { return imap_.find(s)!=imap_.end();}

    // Check if out exists
    bool has_out(const std::string& s) const { return omap_.find(s) != omap_.end();}

    // Get input scheme
    std::vector<std::string> name_in() const;

    // Get output scheme
    std::vector<std::string> name_out() const;
  };

  template<typename MatType>
  void Factory<MatType>::
  add_input(const std::string& s, const MatType& e, bool is_diff) {
    size_t ind = in_.size();
    auto it = imap_.insert(make_pair(s, ind));
    casadi_assert(it.second, "Duplicate input expression \"" + s + "\"");
    is_diff_imap_.insert(make_pair(s, is_diff));
    in_.push_back(e);
    iname_.push_back(s);
  }

  template<typename MatType>
  void Factory<MatType>::
  add_output(const std::string& s, const MatType& e, bool is_diff) {
    size_t ind = out_.size();
    auto it = omap_.insert(make_pair(s, ind));
    casadi_assert(it.second, "Duplicate output expression \"" + s + "\"");
    is_diff_omap_.insert(make_pair(s, is_diff));
    out_.push_back(e);
    oname_.push_back(s);
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
      fwd_imap_.push_back(ss.second);
    } else if (ss.first=="adj") {
      // Reverse mode directional derivative
      casadi_assert(has_out(ss.second), "Cannot process \"" + ss.second + "\""
                                                " (from \"" + s + "\") as output."
                                                " Available: " + join(name_out()) + ".");
      adj_imap_.push_back(ss.second);
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
      fwd_omap_.push_back(ss.second);
    } else if (ss.first=="adj") {
      // Reverse mode directional derivative
      casadi_assert(has_in(ss.second),
        "Cannot process \"" + ss.second + "\" (from \"" + s + "\") as input. "
        "Available: " + join(name_in()) + ".");
      adj_omap_.push_back(ss.second);
    } else if (ss.first=="jac") {
      jac_.push_back(block(ss.second));
    } else if (ss.first=="grad") {
      grad_.push_back(block(ss.second));
    } else if (ss.first=="hess") {
      hess_.push_back(hblock(ss.second));
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
    if (fwd_omap_.empty()) return;
    casadi_assert_dev(!fwd_imap_.empty());

    std::vector<MatType> arg, res;
    std::vector<std::vector<MatType>> seed(1), sens(1);
    // Inputs and forward mode seeds
    for (const std::string& s : fwd_imap_) {
      arg.push_back(in_[imap_[s]]);
      Sparsity sp = is_diff_imap_[s] ? arg.back().sparsity() : Sparsity(arg.back().size());
      seed[0].push_back(MatType::sym("fwd_" + s, sp));
      add_input("fwd:" + s, seed[0].back(), true);
    }
    // Outputs
    for (const std::string& s : fwd_omap_) {
      res.push_back(out_.at(omap_[s]));
    }
    // Calculate directional derivatives
    Dict local_opts = opts;
    local_opts["always_inline"] = true;
    sens = forward(res, arg, seed, local_opts);

    // Get directional derivatives
    for (casadi_int i=0; i<fwd_omap_.size(); ++i) {
      std::string s = fwd_omap_[i];
      Sparsity sp = is_diff_omap_[s] ? res.at(i).sparsity() : Sparsity(res.at(i).size());
      add_output("fwd:" + s, project(sens[0].at(i), sp), is_diff_omap_[s]);
    }
  }

  template<typename MatType>
  void Factory<MatType>::calculate_adj(const Dict& opts) {
    if (adj_omap_.empty()) return;
    casadi_assert_dev(!adj_imap_.empty());
    std::vector<MatType> arg, res;
    std::vector<std::vector<MatType>> seed(1), sens(1);
    // Inputs
    for (const std::string& s : adj_omap_) {
      arg.push_back(in_[imap_[s]]);
    }
    // Outputs and reverse mode seeds
    for (const std::string& s : adj_imap_) {
      res.push_back(out_.at(omap_[s]));
      Sparsity sp = is_diff_omap_[s] ? res.back().sparsity() : Sparsity(res.back().size());
      seed[0].push_back(MatType::sym("adj_" + s, sp));
      add_input("adj:" + s, seed[0].back(), true);
    }
    // Calculate directional derivatives
    Dict local_opts;
    local_opts["always_inline"] = true;
    sens = reverse(res, arg, seed, local_opts);

    // Get directional derivatives
    for (casadi_int i=0; i<adj_omap_.size(); ++i) {
      std::string s = adj_omap_[i];
      Sparsity sp = is_diff_imap_[s] ? arg.at(i).sparsity() : Sparsity(arg.at(i).size());
      add_output("adj:" + s, project(sens[0].at(i), sp), is_diff_imap_[s]);
    }
  }

  template<typename MatType>
  void Factory<MatType>::calculate_jac(const Dict& opts) {
    // Calculate blocks for all non-differentiable inputs and outputs
    for (auto &&b : jac_) {
      std::string s = "jac:" + oname_.at(b.f) + ":" + iname_.at(b.x);
      if (is_diff_omap_.at(oname_.at(b.f)) && is_diff_imap_.at(iname_.at(b.x))) {
        is_diff_omap_[s] = true;
      } else {
        add_output(s, MatType(out_[b.f].numel(), in_[b.x].numel()), false);
      }
    }
    // Calculate regular blocks
    for (auto &&b : jac_) {
      // Get block name, skip if already calculated
      std::string s = "jac:" + oname_[b.f] + ":" + iname_[b.x];
      if (omap_.find(s) != omap_.end()) continue;
      // Find other blocks with the same input, but different (not yet calculated) outputs
      std::vector<MatType> ex;
      std::vector<std::string> all_ex;
      for (auto &&b1 : jac_) {
        // Check if same input
        if (b1.x != b.x) continue;
        // Check if already calculated
        std::string s1 = "jac:" + oname_.at(b1.f) + ":" + iname_[b1.x];
        if (omap_.find(s1) != omap_.end()) continue;
        // Collect expressions
        all_ex.push_back(oname_.at(b1.f));
        ex.push_back(out_.at(b1.f));
      }
      // Now find other blocks with *all* the same outputs, but different inputs
      std::vector<MatType> arg{in_[b.x]};
      std::vector<std::string> all_arg{iname_[b.x]};
      for (auto &&b1 : jac_) {
        // Candidate b1.arg: Check if already added
        bool skip = false;
        for (const std::string& a : all_arg) {
          if (a == iname_[b1.x]) {
            skip = true;
            break;
          }
        }
        if (skip) continue;
        // Check if all blocks corresponding to the same input are needed
        for (const std::string& e : all_ex) {
          std::string s1 = "jac:" + e + ":" + iname_[b1.x];
          if (is_diff_omap_.find(s1) == is_diff_omap_.end() || omap_.find(s1) != omap_.end()) {
            // Block is not requested or has already been calculated
            skip = true;
            break;
          }
        }
        if (skip) continue;
        // Keep candidate
        arg.push_back(in_[b1.x]);
        all_arg.push_back(iname_[b1.x]);
      }
      try {
        // Calculate Jacobian block(s)
        if (ex.size() == 1 && arg.size() == 1) {
          casadi_assert(omap_.find(s) == omap_.end(), "here");
          omap_[s] = out_.size();
          out_.push_back(MatType::jacobian(ex[0], arg[0], opts));
          oname_.push_back(s);
        } else {
          // Calculate Jacobian of all outputs with respect to all inputs
          MatType J = MatType::jacobian(vertcat(ex), vertcat(arg), opts);
          // Split Jacobian into blocks
          std::vector<std::vector<MatType>> J_all = blocksplit(J, offset(ex), offset(arg));
          // Save blocks
          for (size_t i = 0; i < all_ex.size(); ++i) {
            for (size_t j = 0; j < all_arg.size(); ++j) {
              std::string sJ = "jac:" + all_ex[i] + ":" + all_arg[j];
              casadi_assert(omap_.find(sJ) == omap_.end(), "here");
              omap_[sJ] = out_.size();
              out_.push_back(J_all.at(i).at(j));
              oname_.push_back(sJ);
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
      const MatType& ex = out_.at(b.f);
      const MatType& arg = in_[b.x];
      if (is_diff_omap_.at(oname_[b.f]) && is_diff_imap_.at(iname_[b.x])) {
        add_output("grad:" + oname_[b.f] + ":" + iname_[b.x],
          project(gradient(ex, arg, opts), arg.sparsity()), true);
      } else {
        casadi_assert(ex.is_scalar(), "Can only take gradient of scalar expression.");
        add_output("grad:" + oname_[b.f] + ":" + iname_[b.x], MatType(1, arg.numel()), false);
      }
    }
  }

  template<typename MatType>
  void Factory<MatType>::calculate_hess(const Dict& opts, size_t f) {
    // Handle all blocks for this expression
    for (auto &&b : hess_) {
      if (b.f != f) continue;
      // Skip if already calculated
      if (b.calculated) continue;
      // Calculate Hessian blocks
      const MatType& x1 = in_[b.x1];
      const MatType& x2 = in_[b.x2];
      MatType H = b.x1 == b.x2 ? hessian(out_.at(f), x1, opts)
        : jacobian(gradient(out_.at(f), x1), x2);
      add_output(b.s, H, true);
      b.calculated = true;
    }
  }

  template<typename MatType>
  void Factory<MatType>::calculate_hess(const Dict& opts) {
    // Calculate blocks for all non-differentiable inputs and outputs
    for (auto &&b : hess_) {
      if (is_diff_omap_.at(oname_[b.f]) && is_diff_imap_.at(iname_[b.x1])
          && is_diff_imap_.at(iname_[b.x2])) {
        b.calculated = false;
      } else {
        add_output(b.s, MatType(in_[b.x1].numel(), in_[b.x2].numel()), false);
        b.calculated = true;
      }
      // Consistency check
      casadi_assert(out_.at(b.f).is_scalar(),
        "Can only take Hessian of scalar expression.");
    }
    // Calculate regular blocks
    for (auto &&b : hess_) {
      // Skip if already calculated
      if (b.calculated) continue;
      // Calculate all Hessian blocks for b.f
      calculate_hess(opts, b.f);
    }
  }

  template<typename MatType>
  void Factory<MatType>::add_dual(const Function::AuxOut& aux) {
    // Dual variables
    for (auto&& e : omap_) {
      Sparsity sp = is_diff_omap_[e.first] ? out_.at(e.second).sparsity()
        : Sparsity(out_.at(e.second).size());
      add_input("lam:" + e.first, MatType::sym("lam_" + e.first, sp), true);
    }
    // Add linear combinations
    for (auto i : aux) {
      MatType lc = 0;
      for (auto j : i.second) {
        lc += dot(in_.at(imap_.at("lam:" + j)), out_.at(omap_.at(j)));
      }
      add_output(i.first, lc, true);
    }
  }

  template<typename MatType>
  void Factory<MatType>::calculate(const Dict& opts) {
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
    auto it = imap_.find(s);
    casadi_assert(it!=imap_.end(), "Cannot retrieve \"" + s + "\"");
    return in_.at(it->second);
  }

  template<typename MatType>
  MatType Factory<MatType>::get_output(const std::string& s) {
    using namespace std;

    // Quick return if output
    auto it = omap_.find(s);
    if (it!=omap_.end()) return out_.at(it->second);

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
      casadi_warning("Attribute 'sym' has been deprecated. Hessians are symmetric by default.");
      return r;
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
    for (auto i : imap_) {
      ret.push_back(i.first);
    }
    return ret;
  }

  template<typename MatType>
  std::vector<std::string> Factory<MatType>::name_out() const {
    std::vector<std::string> ret;
    for (auto i : omap_) {
      ret.push_back(i.first);
    }
    return ret;
  }

  template<typename MatType>
  size_t Factory<MatType>::imap(const std::string& s) const {
    auto iind = imap_.find(s);
    casadi_assert(iind != imap_.end(),
      "Cannot process \"" + s + "\" as input. Available: " + join(name_out()) + ".");
    return iind->second;
  }

  template<typename MatType>
  size_t Factory<MatType>::omap(const std::string& s) const {
    auto oind = omap_.find(s);
    casadi_assert(oind != omap_.end(),
      "Cannot process \"" + s + "\" as output. Available: " + join(name_out()) + ".");
    return oind->second;
  }

  template<typename MatType>
  Block Factory<MatType>::block(const std::string& s) const {
    Block b;
    size_t pos = s.find(':');
    if (pos<s.size()) {
      b.f = omap(s.substr(0, pos));
      b.x = imap(s.substr(pos+1, std::string::npos));
    }
    return b;
  }

  template<typename MatType>
  HBlock Factory<MatType>::hblock(const std::string& s) const {
    HBlock b;
    b.s = "hess:" + s;
    size_t pos1 = s.find(':');
    if (pos1 < s.size()) {
      size_t pos2 = s.find(':', pos1 + 1);
      if (pos2 < s.size()) {
        b.f = omap(s.substr(0, pos1));
        b.x1 = imap(s.substr(pos1 + 1, pos2 - pos1 - 1));
        b.x2 = imap(s.substr(pos2 + 1, std::string::npos));
      }
    }
    return b;
  }

} // namespace casadi
/// \endcond

#endif // CASADI_FACTORY_HPP
