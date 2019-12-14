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


#include "casadi_low.hpp"

using namespace std;

namespace casadi {

  Low::Low(const MX& v, const MX& p, const Dict& opts) {
    casadi_assert_dev(v.is_vector() && v.is_dense());
    casadi_assert_dev(p.is_vector() && v.is_dense());
    set_dep(v, p);
    set_sparsity(p.sparsity());

    std::string lookup_mode = "auto";
    for (auto&& e : opts) {
      if (e.first=="lookup_mode") {
        lookup_mode = e.second.to_string();
      } else {
        casadi_error("Unrecongnized option: " + str(e.first));
      }
    }

    lookup_mode_ = interpret_lookup_mode(lookup_mode, v.numel());
  }

  std::string Low::lookup_mode_from_enum(casadi_int lookup_mode) {
    switch (lookup_mode) {
      case LOOKUP_LINEAR:
        return "linear";
      case LOOKUP_EXACT:
        return "exact";
      case LOOKUP_BINARY:
        return "binary";
      default:
        casadi_assert_dev(false);
    }
  }

  casadi_int Low::interpret_lookup_mode(const std::string& lookup_mode, casadi_int n) {
    if (lookup_mode=="auto") {
      if (n>100) return interpret_lookup_mode("binary", n);
      return interpret_lookup_mode("linear", n);
    } else if (lookup_mode=="binary") {
      return LOOKUP_BINARY;
    } else if (lookup_mode=="linear") {
      return LOOKUP_LINEAR;
    } else if (lookup_mode=="exact") {
      return LOOKUP_EXACT;
    } else {
      casadi_error("Invalid lookup mode '" + lookup_mode + "'. "
        "Available modes: linear|binary|exact|auto");
    }
  }

  std::string Low::disp(const std::vector<std::string>& arg) const {
    return "low(" + arg.at(0) + ", " + arg.at(1) + ")";
  }

  int Low::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    for (casadi_int i=0;i<dep(1).nnz();++i) {
      res[0][i] = casadi_low(arg[1][i], arg[0], dep(0).nnz(), lookup_mode_);
    }
    return 0;
  }

  void Low::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = low(arg[0], arg[1]);
  }

  void Low::ad_forward(const std::vector<std::vector<MX> >& fseed,
                     std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = 0;
    }
  }

  void Low::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                     std::vector<std::vector<MX> >& asens) const {
  }

  int Low::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    fill_n(res[0], nnz(), 0);
    return 0;
  }

  int Low::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    fill_n(res[0], nnz(), 0);
    return 0;
  }

  void Low::generate(CodeGenerator& g,
                      const std::vector<casadi_int>& arg,
                      const std::vector<casadi_int>& res) const {
    casadi_int n = dep(1).nnz();
    casadi_int ng = dep(0).nnz();
    g.local("cr", "const casadi_real", "*");
    g.local("rr", "casadi_real", "*");
    g << "for (cr=" << g.work(arg[1], n) << ", rr=" << g.work(res[0], n);
    g << ";cr!=" << g.work(arg[1], n) << "+" << n << ";++cr) ";
    g << "*rr++ = ";
    g << g.low("*cr", g.work(arg[0], ng), ng, lookup_mode_) << "\n";
  }

  void Low::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("Low::lookup_mode", static_cast<casadi_int>(lookup_mode_));
  }

  Low::Low(DeserializingStream& s) : MXNode(s) {
    casadi_int lookup_mode;
    s.unpack("Low::lookup_mode", lookup_mode);
    lookup_mode_ = static_cast<casadi_int>(lookup_mode);
  }

} // namespace casadi
