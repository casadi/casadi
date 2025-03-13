/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#include "monitor.hpp"

namespace casadi {

  Monitor::Monitor(const MX& x, const std::string& comment) : comment_(comment) {
    casadi_assert_dev(x.nnz()>0);
    set_dep(x);
    set_sparsity(x.sparsity());
  }

  std::string Monitor::disp(const std::vector<std::string>& arg) const {
    return "monitor(" + arg.at(0) + ", " + comment_ + ")";
  }

  void Monitor::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0].monitor(comment_);
  }

  void Monitor::propagate_interval(const std::vector<MX>& arg_L, const std::vector<MX>& arg_R,
      std::vector<MX>& res_L, std::vector<MX>& res_R) const {
    eval_mx(arg_L, res_L);
    eval_mx(arg_R, res_R);
  }

  void Monitor::ad_forward(const std::vector<std::vector<MX> >& fseed,
                        std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d=0; d<fsens.size(); ++d) {
      std::stringstream ss;
      ss << "fwd(" << d << ") of " << comment_;
      fsens[d][0] = fseed[d][0].monitor(ss.str());
    }
  }

  void Monitor::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                        std::vector<std::vector<MX> >& asens) const {
    for (casadi_int d=0; d<aseed.size(); ++d) {
      std::stringstream ss;
      ss << "adj(" << d << ") of " << comment_;
      asens[d][0] += aseed[d][0].monitor(ss.str());
    }
  }

  int Monitor::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    if (arg[0]!=res[0]) {
      std::copy(arg[0], arg[0]+nnz(), res[0]);
    }
    return 0;
  }

  int Monitor::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    // Print comment
    uout() << comment_ << ":" << std::endl;
    uout() << "[";
    casadi_int n = nnz();
    for (casadi_int i=0; i<n; ++i) {
      if (i!=0) uout() << ", ";
      uout() << arg[0][i];
    }
    uout() << "]" << std::endl;

    // Perform operation
    if (arg[0]!=res[0]) {
      std::copy(arg[0], arg[0]+n, res[0]);
    }
    return 0;
  }

  int Monitor::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    if (arg[0]!=res[0]) {
      std::copy(arg[0], arg[0]+nnz(), res[0]);
    }
    return 0;
  }

  int Monitor::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *a = arg[0];
    bvec_t *r = res[0];
    casadi_int n = nnz();
    if (a != r) {
      for (casadi_int i=0; i<n; ++i) {
        *a++ |= *r;
        *r++ = 0;
      }
    }
    return 0;
  }

  void Monitor::generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res,
                          const std::vector<bool>& arg_is_ref,
                          std::vector<bool>& res_is_ref) const {
    // Print comment
    g.local("cr", "const casadi_real", "*");
    g.local("i", "casadi_int");
    g << g.printf(comment_ + "\\n[") << "\n"
      << "  for (i=0, cr=" << g.work(arg[0], dep(0).nnz(), arg_is_ref[0])
      << "; i!=" << nnz() << "; ++i) {\n"
      << "    if (i!=0) " << g.printf(", ") << "\n"
      << "    " << g.printf("%g", "*cr++") << "\n"
      << "  }\n"
      << "  " << g.printf("]\\n") << "\n";

    generate_copy(g, arg, res, arg_is_ref, res_is_ref, 0);
  }

  void Monitor::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("Monitor::comment", comment_);
  }

  Monitor::Monitor(DeserializingStream& s) : MXNode(s) {
    s.unpack("Monitor::comment", comment_);
  }

} // namespace casadi
