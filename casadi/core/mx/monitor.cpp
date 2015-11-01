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


#include "monitor.hpp"

using namespace std;

namespace casadi {

  Monitor::Monitor(const MX& x, const std::string& comment) : comment_(comment) {
    casadi_assert(x.nnz()>0);
    setDependencies(x);
    setSparsity(x.sparsity());
  }

  std::string Monitor::print(const std::vector<std::string>& arg) const {
    return "monitor(" + arg.at(0) + ", " + comment_ + ")";
  }

  void Monitor::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = arg[0].monitor(comment_);
  }

  void Monitor::evalFwd(const std::vector<std::vector<MX> >& fseed,
                        std::vector<std::vector<MX> >& fsens) {
    for (int d=0; d<fsens.size(); ++d) {
      stringstream ss;
      ss << "fwd(" << d << ") of " << comment_;
      fsens[d][0] = fseed[d][0].monitor(ss.str());
    }
  }

  void Monitor::evalAdj(const std::vector<std::vector<MX> >& aseed,
                        std::vector<std::vector<MX> >& asens) {
    for (int d=0; d<aseed.size(); ++d) {
      stringstream ss;
      ss << "adj(" << d << ") of " << comment_;
      asens[d][0] += aseed[d][0].monitor(ss.str());
    }
  }

  void Monitor::evalSX(void* mem, const SXElem** arg, SXElem** res, int* iw, SXElem* w) {
    if (arg[0]!=res[0]) {
      copy(arg[0], arg[0]+nnz(), res[0]);
    }
  }

  void Monitor::evalD(void* mem, const double** arg, double** res, int* iw, double* w) {
    // Print comment
    userOut() << comment_ << ":" << endl;
    userOut() << "[";
    int n = nnz();
    for (int i=0; i<n; ++i) {
      if (i!=0) userOut() << ", ";
      userOut() << arg[0][i];
    }
    userOut() << "]" << endl;

    // Perform operation
    if (arg[0]!=res[0]) {
      copy(arg[0], arg[0]+n, res[0]);
    }
  }

  void Monitor::spFwd(void* mem, const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    if (arg[0]!=res[0]) {
      copy(arg[0], arg[0]+nnz(), res[0]);
    }
  }

  void Monitor::spAdj(void* mem, bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    bvec_t *a = arg[0];
    bvec_t *r = res[0];
    int n = nnz();
    if (a != r) {
      for (int i=0; i<n; ++i) {
        *a++ |= *r;
        *r++ = 0;
      }
    }
  }

  void Monitor::generate(CodeGenerator& g, const std::string& mem,
                         const std::vector<int>& arg, const std::vector<int>& res) const {
    // Print comment
    g.body << "  " << g.printf(comment_ + "\\n[") << endl
           << "  for (i=0, rr=" << g.work(arg[0], dep(0).nnz())
           << "; i!=" << nnz() << "; ++i) {" << endl
           << "    if (i!=0) " << g.printf(", ") << endl
           << "    " << g.printf("%g", "*rr++") << endl
           << "  }" << endl
           << "  " << g.printf("]\\n") << endl;

    // Copy if not inplace
    if (arg[0]!=res[0]) {
      if (nnz()==1) {
        g.body << "  " << g.workel(res[0]) << " = " << g.workel(arg[0]) << ";" << endl;
      } else {
        g.body << "  " << g.copy_n(g.work(arg[0], nnz()), nnz(), g.work(res[0], nnz())) << endl;
      }
    }
  }

} // namespace casadi
