/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
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

#include <algorithm>
#include "function.hpp"

/// \cond INTERNAL

namespace casadi {

  // A Jacobian or gradient block
  struct Block {
    size_t f, x;
    std::string s;
    bool calculated;
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
    std::vector<bool> is_diff_in_, is_diff_out_;

    // Forward mode directional derivatives
    std::vector<size_t> fwd_in_, fwd_out_;

    // Reverse mode directional derivatives
    std::vector<size_t> adj_in_, adj_out_;

    // Jacobian/gradient blocks
    std::vector<Block> jac_, grad_;

    // Hessian blocks
    std::vector<HBlock> hess_;

    // Read a Jacobian or gradient block
    Block block(const std::string& s1, const std::string& s) const;

    // Read a Hessian block
    HBlock hblock(const std::string& s1, const std::string& s) const;

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

    // Find a Jacobian block
    std::vector<Block>::iterator find_jac(size_t f, size_t x);

    // Find a Hessian block
    std::vector<HBlock>::iterator find_hess(size_t f, size_t x1, size_t x2);

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
    const std::vector<std::string>& iname() const {return iname_;}
    std::vector<std::string> iname(const std::vector<size_t>& ind) const;

    // Get output scheme
    const std::vector<std::string>& oname() const {return oname_;}
    std::vector<std::string> oname(const std::vector<size_t>& ind) const;
  };

  template<typename MatType>
  void Factory<MatType>::
  add_input(const std::string& s, const MatType& e, bool is_diff) {
    size_t ind = in_.size();
    auto it = imap_.insert(std::make_pair(s, ind));
    casadi_assert(it.second, "Duplicate input expression \"" + s + "\"");
    is_diff_in_.push_back(is_diff);
    in_.push_back(e);
    iname_.push_back(s);
  }

  template<typename MatType>
  void Factory<MatType>::
  add_output(const std::string& s, const MatType& e, bool is_diff) {
    size_t ind = out_.size();
    auto it = omap_.insert(std::make_pair(s, ind));
    casadi_assert(it.second, "Duplicate output expression \"" + s + "\"");
    is_diff_out_.push_back(is_diff);
    out_.push_back(e);
    oname_.push_back(s);
  }

  template<typename MatType>
  std::string Factory<MatType>::
  request_input(const std::string& s) {
    // Add input if not already available
    if (!has_in(s)) {
      // Get prefix
      casadi_assert(has_prefix(s), "Cannot process \"" + s + "\" as input."
        " Available: " + join(iname()) + ".");
      std::pair<std::string, std::string> ss = split_prefix(s);
      // Process specific prefixes
      if (ss.first=="fwd") {
        // Forward mode directional derivative
        fwd_in_.push_back(imap(ss.second));
      } else if (ss.first=="adj") {
        // Reverse mode directional derivative
        adj_in_.push_back(omap(ss.second));
      }
    }
    // Replace colons with underscore
    std::string ret = s;
    std::replace(ret.begin(), ret.end(), ':', '_');
    return ret;
  }

  template<typename MatType>
  std::string Factory<MatType>::
  request_output(const std::string& s) {
    // Quick return if already available
    if (has_out(s)) return s;

    // Get prefix
    casadi_assert(has_prefix(s), "Cannot process \"" + s + "\" as output."
      " Available: " + join(oname()) + ".");
    std::pair<std::string, std::string> ss = split_prefix(s);

    if (ss.first=="fwd") {
      fwd_out_.push_back(omap(ss.second));
    } else if (ss.first=="adj") {
      adj_out_.push_back(imap(ss.second));
    } else if (ss.first=="jac") {
      jac_.push_back(block(ss.second, s));
    } else if (ss.first=="grad") {
      grad_.push_back(block(ss.second, s));
    } else if (ss.first=="hess") {
      hess_.push_back(hblock(ss.second, s));
    } else {
      // Assume attribute
      request_output(ss.second);
    }

    // Replace colons with underscore
    std::string ret = s;
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
    for (size_t iind : fwd_in_) {
      arg.push_back(in_[iind]);
      Sparsity sp = is_diff_in_.at(iind) ? arg.back().sparsity() : Sparsity(arg.back().size());
      seed[0].push_back(MatType::sym("fwd_" + iname_[iind], sp));
      add_input("fwd:" + iname_[iind], seed[0].back(), true);
    }
    // Outputs
    for (size_t oind : fwd_out_) res.push_back(out_.at(oind));
    // Calculate directional derivatives
    Dict local_opts = opts;
    local_opts["always_inline"] = true;
    sens = forward(res, arg, seed, local_opts);

    // Get directional derivatives
    for (size_t i = 0; i < fwd_out_.size(); ++i) {
      std::string s = oname_.at(fwd_out_[i]);
      Sparsity sp = is_diff_out_.at(fwd_out_[i]) ? res.at(i).sparsity()
        : Sparsity(res.at(i).size());
      add_output("fwd:" + s, project(sens[0].at(i), sp), is_diff_out_.at(fwd_out_[i]));
    }
  }

  template<typename MatType>
  void Factory<MatType>::calculate_adj(const Dict& opts) {
    if (adj_out_.empty()) return;
    casadi_assert_dev(!adj_in_.empty());
    std::vector<MatType> arg, res;
    std::vector<std::vector<MatType>> seed(1), sens(1);
    // Inputs
    for (size_t ind : adj_out_) arg.push_back(in_[ind]);
    // Outputs and reverse mode seeds
    for (size_t ind : adj_in_) {
      res.push_back(out_.at(ind));
      Sparsity sp = is_diff_out_.at(ind) ? res.back().sparsity() : Sparsity(res.back().size());
      seed[0].push_back(MatType::sym("adj_" + oname_[ind], sp));
      add_input("adj:" + oname_[ind], seed[0].back(), true);
    }
    // Calculate directional derivatives
    Dict local_opts;
    local_opts["always_inline"] = true;
    sens = reverse(res, arg, seed, local_opts);

    // Get directional derivatives
    for (size_t i=0; i < adj_out_.size(); ++i) {
      std::string s = iname_[adj_out_[i]];
      Sparsity sp = is_diff_in_.at(adj_out_[i]) ? arg.at(i).sparsity() : Sparsity(arg.at(i).size());
      add_output("adj:" + s, project(sens[0].at(i), sp), is_diff_in_.at(adj_out_[i]));
    }
  }

  template<typename MatType>
  std::vector<Block>::iterator Factory<MatType>::find_jac(size_t f, size_t x) {
    for (std::vector<Block>::iterator it = jac_.begin(); it != jac_.end(); ++it) {
      if (it->f == f && it->x == x) return it;
    }
    // Not in list
    return jac_.end();
  }

  template<typename MatType>
  std::vector<HBlock>::iterator Factory<MatType>::find_hess(size_t f, size_t x1, size_t x2) {
    for (std::vector<HBlock>::iterator it = hess_.begin(); it != hess_.end(); ++it) {
      if (it->f == f && it->x1 == x1 && it->x2 == x2) return it;
    }
    // Not in list
    return hess_.end();
  }

  template<typename MatType>
  void Factory<MatType>::calculate_jac(const Dict& opts) {
    // Calculate blocks for all non-differentiable inputs and outputs
    for (auto &&b : jac_) {
      if (is_diff_out_.at(b.f) && is_diff_in_.at(b.x)) {
        b.calculated = false;
      } else {
        add_output(b.s, MatType(out_[b.f].numel(), in_[b.x].numel()), false);
        b.calculated = true;
      }
    }
    // Calculate regular blocks
    for (auto &&b : jac_) {
      // Skip if already calculated
      if (b.calculated) continue;
      // Find other blocks with the same input, but different (not yet calculated) outputs
      std::vector<size_t> all_f;
      for (auto &&b1 : jac_) {
        if (b1.x == b.x && !b1.calculated) all_f.push_back(b1.f);
      }
      // Now find other blocks with *all* the same outputs, but different inputs
      std::vector<size_t> all_x{b.x};
      for (auto &&b1 : jac_) {
        // Candidate b1.arg: Check if already added
        if (std::count(all_x.begin(), all_x.end(), b1.x)) continue;
        // Skip if all block are not requested or any block has already been calculated
        bool skip = false;
        for (size_t f1 : all_f) {
          auto it = find_jac(f1, b1.x);
          if (it == jac_.end() || it->calculated) {
            skip = true;
            break;
          }
        }
        if (skip) continue;
        // Keep candidate
        all_x.push_back(b1.x);
      }
      try {
        // Calculate Jacobian block(s)
        if (all_f.size() == 1 && all_x.size() == 1) {
          // Single block
          add_output(b.s, MatType::jacobian(out_[b.f], in_[b.x], opts), true);
          b.calculated = true;
        } else {
          // Sort blocks
          std::sort(all_x.begin(), all_x.end());
          std::sort(all_f.begin(), all_f.end());
          // Collect components
          std::vector<MatType> x(all_x.size()), f(all_f.size());
          for (size_t i = 0; i < x.size(); ++i) x[i] = in_.at(all_x[i]);
          for (size_t i = 0; i < f.size(); ++i) f[i] = out_.at(all_f[i]);
          // Calculate Jacobian of all outputs with respect to all inputs
          MatType J = MatType::jacobian(veccat(f), veccat(x), opts);
          // Split Jacobian into blocks
          std::vector<std::vector<MatType>> J_all = blocksplit(J, offset(f), offset(x));
          // Save blocks
          for (size_t i = 0; i < all_f.size(); ++i) {
            for (size_t j = 0; j < all_x.size(); ++j) {
              auto J_it = find_jac(all_f[i], all_x[j]);
              if (J_it != jac_.end()) {
                add_output(J_it->s, J_all.at(i).at(j), true);
                J_it->calculated = true;
              }
            }
          }
        }
      } catch (std::exception& e) {
        std::stringstream ss;
        ss << "Calculating Jacobian of " << oname(all_f) << " w.r.t. " << iname(all_x)
          << ": " << e.what();
        casadi_error(ss.str());
      }
    }
  }

  template<typename MatType>
  void Factory<MatType>::calculate_grad(const Dict& opts) {
    for (auto &&b : grad_) {
      const MatType& ex = out_.at(b.f);
      const MatType& arg = in_[b.x];
      if (is_diff_out_.at(b.f) && is_diff_in_.at(b.x)) {
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
      // Find other blocks with one of the arguments matching
      std::vector<size_t> all_x1;
      for (auto &&b1 : hess_) {
        if (b1.f == b.f && !b1.calculated) {
          if (b1.x1 == b.x1) {
            // Block found
            all_x1.push_back(b1.x2);
          } else if (b1.x2 == b.x1) {
            // Opposite block found
            all_x1.push_back(b1.x1);
          }
        }
      }
      // Find other blocks with both of the arguments matching
      std::vector<size_t> all_x2;
      for (auto &&b1 : hess_) {
        if (b1.f != f || b1.calculated) continue;
        // Can either b1.x1 or b1.x2 be added to all_x2?
        for (bool test_x1 : {false, true}) {
          size_t cand = test_x1 ? b1.x1 : b1.x2;
          size_t other = test_x1 ? b1.x2 : b1.x1;
          bool cand_ok = true;
          bool other_ok = false;
          // Skip if already in all_x2
          if (std::count(all_x2.begin(), all_x2.end(), cand)) continue;
          // Loop over existing entries in x1
          for (size_t a : all_x1) {
            // The other argument must already be in all_x1
            if (other == a) other_ok = true;
            // Is block not requested?
            auto it = find_hess(f, a, cand);
            if (it == hess_.end() || it->calculated) {
              // Also check mirror block, if there is one
              if (a != cand) {
                it = find_hess(f, cand, a);
                if (it != hess_.end() && !it->calculated) continue;
              }
              // Not a candidate
              cand_ok = false;
              break;
            }
          }
          // Keep candidate
          if (cand_ok && other_ok) all_x2.push_back(cand);
        }
      }
      // Calculate Hessian blocks
      try {
        if (all_x1.size() == 1 && all_x2.size() == 1) {
          // Single block
          MatType H = b.x1 == b.x2 ? hessian(out_.at(f), in_[b.x1], opts)
            : jacobian(gradient(out_.at(f), in_[b.x1]), in_[b.x2]);
          add_output(b.s, H, true);
          b.calculated = true;
        } else {
          // Sort blocks
          std::sort(all_x1.begin(), all_x1.end());
          std::sort(all_x2.begin(), all_x2.end());
          // Symmetric extended Hessian?
          bool symmetric = all_x1 == all_x2;
          // Collect components
          std::vector<MatType> x1(all_x1.size()), x2(all_x2.size());
          for (size_t i = 0; i < x1.size(); ++i) x1[i] = in_.at(all_x1[i]);
          for (size_t i = 0; i < x2.size(); ++i) x2[i] = in_.at(all_x2[i]);
          // Calculate extended Hessian
          MatType H;
          if (symmetric) {
            H = hessian(out_.at(f), vertcat(x1));
          } else {
            H = jacobian(gradient(out_.at(f), vertcat(x1)), vertcat(x2));
          }
          // Split into blocks
          std::vector<std::vector<MatType>> H_all = blocksplit(H, offset(x1), offset(x2));
          // Collect Hessian blocks
          for (auto &&b1 : hess_) {
            if (b1.f == f && !b1.calculated) {
              // Find arguments in all_x1 and all_x2
              auto it_x1 = std::find(all_x1.begin(), all_x1.end(), b1.x1);
              auto it_x2 = std::find(all_x2.begin(), all_x2.end(), b1.x2);
              if (it_x1 != all_x1.end() && it_x2 != all_x2.end()) {
                // Block located
                const MatType& Hb = H_all.at(it_x1 - all_x1.begin()).at(it_x2 - all_x2.begin());
                add_output(b1.s, Hb, true);
                b1.calculated = true;
              } else if (!symmetric) {
                // Check mirror block
                it_x1 = std::find(all_x1.begin(), all_x1.end(), b1.x2);
                it_x2 = std::find(all_x2.begin(), all_x2.end(), b1.x1);
                if (it_x1 != all_x1.end() && it_x2 != all_x2.end()) {
                  // Transpose of block located
                  const MatType& Hb = H_all.at(it_x1 - all_x1.begin()).at(it_x2 - all_x2.begin());
                  add_output(b1.s, Hb.T(), true);
                  b1.calculated = true;
                }
              }
            }
          }
        }
      } catch (std::exception& e) {
        std::stringstream ss;
        ss << "Calculating Hessian of " << oname_.at(f) << " w.r.t. " << iname(all_x1) << " and "
          << iname(all_x2) << ": " << e.what();
        casadi_error(ss.str());
      }
    }
  }

  template<typename MatType>
  void Factory<MatType>::calculate_hess(const Dict& opts) {
    // Calculate blocks for all non-differentiable inputs and outputs
    for (auto &&b : hess_) {
      if (is_diff_out_.at(b.f) && is_diff_in_.at(b.x1) && is_diff_in_.at(b.x2)) {
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
    for (size_t k = 0; k < out_.size(); ++k) {
      Sparsity sp = is_diff_out_[k] ? out_.at(k).sparsity() : Sparsity(out_.at(k).size());
      add_input("lam:" + oname_[k], MatType::sym("lam_" + oname_[k], sp), true);
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
    // Quick return if output
    auto it = omap_.find(s);
    if (it!=omap_.end()) return out_.at(it->second);

    // Assume attribute
    casadi_assert(has_prefix(s), "Cannot process \"" + s + "\"");
    std::pair<std::string, std::string> ss = split_prefix(s);
    std::string a = ss.first;
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
    return std::make_pair(s.substr(0, pos), s.substr(pos+1, std::string::npos));
  }

  template<typename MatType>
  std::vector<std::string> Factory<MatType>::iname(const std::vector<size_t>& ind) const {
    std::vector<std::string> ret;
    for (size_t i : ind) ret.push_back(iname_.at(i));
    return ret;
  }

  template<typename MatType>
  std::vector<std::string> Factory<MatType>::oname(const std::vector<size_t>& ind) const {
    std::vector<std::string> ret;
    for (size_t i : ind) ret.push_back(oname_.at(i));
    return ret;
  }

  template<typename MatType>
  size_t Factory<MatType>::imap(const std::string& s) const {
    auto iind = imap_.find(s);
    casadi_assert(iind != imap_.end(),
      "Cannot process \"" + s + "\" as input. Available: " + join(oname()) + ".");
    return iind->second;
  }

  template<typename MatType>
  size_t Factory<MatType>::omap(const std::string& s) const {
    auto oind = omap_.find(s);
    casadi_assert(oind != omap_.end(),
      "Cannot process \"" + s + "\" as output. Available: " + join(oname()) + ".");
    return oind->second;
  }

  template<typename MatType>
  Block Factory<MatType>::block(const std::string& s2, const std::string& s) const {
    Block b;
    b.s = s;
    size_t pos = s2.find(':');
    if (pos < s2.size()) {
      b.f = omap(s2.substr(0, pos));
      b.x = imap(s2.substr(pos+1, std::string::npos));
    }
    return b;
  }

  template<typename MatType>
  HBlock Factory<MatType>::hblock(const std::string& s2, const std::string& s) const {
    HBlock b;
    b.s = s;
    size_t pos1 = s2.find(':');
    if (pos1 < s2.size()) {
      size_t pos2 = s2.find(':', pos1 + 1);
      if (pos2 < s2.size()) {
        b.f = omap(s2.substr(0, pos1));
        b.x1 = imap(s2.substr(pos1 + 1, pos2 - pos1 - 1));
        b.x2 = imap(s2.substr(pos2 + 1, std::string::npos));
      }
    }
    return b;
  }

} // namespace casadi
/// \endcond

#endif // CASADI_FACTORY_HPP
