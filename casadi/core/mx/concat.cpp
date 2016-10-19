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


#include "concat.hpp"
#include "../std_vector_tools.hpp"

using namespace std;

namespace casadi {

  Concat::Concat(const vector<MX>& x) {
    setDependencies(x);
  }

  Concat::~Concat() {
  }

  void Concat::eval(const double** arg, double** res, int* iw, double* w, int mem) const {
    evalGen<double>(arg, res, iw, w);
  }

  void Concat::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    evalGen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  void Concat::evalGen(const T* const* arg, T* const* res, int* iw, T* w) const {
    T* r = res[0];
    for (int i=0; i<ndep(); ++i) {
      int n = dep(i).nnz();
      copy(arg[i], arg[i]+n, r);
      r += n;
    }
  }

  void Concat::sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    bvec_t *res_ptr = res[0];
    for (int i=0; i<ndep(); ++i) {
      int n_i = dep(i).nnz();
      const bvec_t *arg_i_ptr = arg[i];
      copy(arg_i_ptr, arg_i_ptr+n_i, res_ptr);
      res_ptr += n_i;
    }
  }

  void Concat::sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    bvec_t *res_ptr = res[0];
    for (int i=0; i<ndep(); ++i) {
      int n_i = dep(i).nnz();
      bvec_t *arg_i_ptr = arg[i];
      for (int k=0; k<n_i; ++k) {
        *arg_i_ptr++ |= *res_ptr;
        *res_ptr++ = 0;
      }
    }
  }

  void Concat::generate(CodeGenerator& g, const std::string& mem,
                        const std::vector<int>& arg, const std::vector<int>& res) const {
    g.body << "  rr=" << g.work(res[0], nnz()) << ";" << endl;
    for (int i=0; i<arg.size(); ++i) {
      int nz = dep(i).nnz();
      if (nz==1) {
        g.body << "  *rr++ = " << g.workel(arg[i]) << ";" << endl;
      } else if (nz!=0) {
        g.body << "  for (i=0, " << "cs=" << g.work(arg[i], nz) << "; "
               << "i<" << nz << "; ++i) *rr++ = *cs++;" << endl;
      }
    }
  }

  MX Concat::getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const {
    // Get the first nonnegative nz
    int nz_test = -1;
    for (auto&& i : nz) {
      if (i>=0) {
        nz_test = i;
        break;
      }
    }

    // Quick return if none
    if (nz_test<0) return MX::zeros(sp);

    // Find out to which dependency it might depend
    int begin=0, end=0;
    int i;
    for (i=0; i<ndep(); ++i) {
      begin = end;
      end += dep(i).nnz();
      if (nz_test < end) break;
    }

    // Check if any nz refer to a different nonzero
    for (auto&& j : nz) {
      if (j>=0 && (j < begin || j >= end)) {

        // Fallback to the base class
        return MXNode::getGetNonzeros(sp, nz);
      }
    }

    // All nz refer to the same dependency, update the nonzero indices
    if (begin==0) {
      return dep(i)->getGetNonzeros(sp, nz);
    } else {
      vector<int> nz_new(nz);
      for (auto&& j : nz_new) if (j>=0) j -= begin;
      return dep(i)->getGetNonzeros(sp, nz_new);
    }
  }


  Diagcat::Diagcat(const std::vector<MX>& x) : Concat(x) {
    casadi_assert(x.size()>1);
    std::vector<Sparsity> sp(x.size());
    for (int i=0; i<x.size(); ++i) sp[i] = x[i].sparsity();
    setSparsity(diagcat(sp));
  }

  std::string Diagcat::print(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << "diagcat(" << arg.at(0);
    for (int i=1; i<ndep(); ++i) ss << ", " << arg.at(i);
    ss << ")";
    return ss.str();
  }

  void Diagcat::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = diagcat(arg);
  }

  void Diagcat::evalFwd(const std::vector<std::vector<MX> >& fseed,
                        std::vector<std::vector<MX> >& fsens) {
    int nfwd = fsens.size();
    for (int d = 0; d<nfwd; ++d) fsens[d][0] = diagcat(fseed[d]);
  }

  std::pair<std::vector<int>, std::vector<int> > Diagcat::offset() const {
    vector<int> offset1(ndep()+1, 0);
    vector<int> offset2(ndep()+1, 0);
    for (int i=0; i<ndep(); ++i) {
      int ncol = dep(i).sparsity().size2();
      int nrow = dep(i).sparsity().size1();
      offset2[i+1] = offset2[i] + ncol;
      offset1[i+1] = offset1[i] + nrow;
    }
    return make_pair(offset1, offset2);
  }

  void Diagcat::evalAdj(const std::vector<std::vector<MX> >& aseed,
                        std::vector<std::vector<MX> >& asens) {
    // Get offsets for each row and column
    auto off = offset();

    // Adjoint sensitivities
    int nadj = aseed.size();
    for (int d=0; d<nadj; ++d) {
      vector<MX> s = diagsplit(aseed[d][0], off.first, off.second);
      for (int i=0; i<ndep(); ++i) {
        asens[d][i] += s[i];
      }
    }
  }

  Horzcat::Horzcat(const std::vector<MX>& x) : Concat(x) {
    casadi_assert(x.size()>1);
    std::vector<Sparsity> sp(x.size());
    for (int i=0; i<x.size(); ++i)
      sp[i] = x[i].sparsity();
    setSparsity(horzcat(sp));
  }

  std::string Horzcat::print(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << "horzcat(" << arg.at(0);
    for (int i=1; i<ndep(); ++i) ss << ", " << arg.at(i);
    ss << ")";
    return ss.str();
  }

  void Horzcat::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = horzcat(arg);
  }

  void Horzcat::evalFwd(const std::vector<std::vector<MX> >& fseed,
                        std::vector<std::vector<MX> >& fsens) {
    int nfwd = fsens.size();
    for (int d = 0; d<nfwd; ++d) {
      fsens[d][0] = horzcat(fseed[d]);
    }
  }

  std::vector<int> Horzcat::offset() const {
    vector<int> col_offset(ndep()+1, 0);
    for (int i=0; i<ndep(); ++i) {
      int ncol = dep(i).sparsity().size2();
      col_offset[i+1] = col_offset[i] + ncol;
    }
    return col_offset;
  }

  void Horzcat::evalAdj(const std::vector<std::vector<MX> >& aseed,
                        std::vector<std::vector<MX> >& asens) {
    // Get offsets for each column
    vector<int> col_offset = offset();

    // Adjoint sensitivities
    int nadj = aseed.size();
    for (int d=0; d<nadj; ++d) {
      vector<MX> s = horzsplit(aseed[d][0], col_offset);
      for (int i=0; i<ndep(); ++i) {
        asens[d][i] += s[i];
      }
    }
  }

  Vertcat::Vertcat(const std::vector<MX>& x) : Concat(x) {
    casadi_assert(x.size()>1);
    std::vector<Sparsity> sp(x.size());
    for (int i=0; i<x.size(); ++i) sp[i] = x[i].sparsity();
    setSparsity(vertcat(sp));
  }

  std::string Vertcat::print(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << "vertcat(" << arg.at(0);
    for (int i=1; i<ndep(); ++i) ss << ", " << arg.at(i);
    ss << ")";
    return ss.str();
  }

  void Vertcat::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = vertcat(arg);
  }

  void Vertcat::evalFwd(const std::vector<std::vector<MX> >& fseed,
                        std::vector<std::vector<MX> >& fsens) {
    int nfwd = fsens.size();
    for (int d = 0; d<nfwd; ++d) {
      fsens[d][0] = vertcat(fseed[d]);
    }
  }

  std::vector<int> Vertcat::offset() const {
    vector<int> row_offset(ndep()+1, 0);
    for (int i=0; i<ndep(); ++i) {
      int nrow = dep(i).sparsity().size1();
      row_offset[i+1] = row_offset[i] + nrow;
    }
    return row_offset;
  }

  void Vertcat::evalAdj(const std::vector<std::vector<MX> >& aseed,
                        std::vector<std::vector<MX> >& asens) {
    // Get offsets for each row
    vector<int> row_offset = offset();

    // Adjoint sensitivities
    int nadj = aseed.size();
    for (int d=0; d<nadj; ++d) {
      vector<MX> s = vertsplit(aseed[d][0], row_offset);
      for (int i=0; i<ndep(); ++i) {
        asens[d][i] += s[i];
      }
    }
  }

  bool Concat::is_valid_input() const {
    for (int i=0; i<ndep(); ++i) {
      if (!dep(i)->is_valid_input()) return false;
    }
    return true;
  }

  int Concat::n_primitives() const {
    int nprim = 0;
    for (int i=0; i<ndep(); ++i) {
      nprim +=  dep(i)->n_primitives();
    }
    return nprim;
  }

  void Horzcat::split_primitives(const MX& x, std::vector<MX>::iterator& it) const {
    vector<MX> s = horzsplit(x, offset());
    for (int i=0; i<s.size(); ++i) {
      dep(i)->split_primitives(s[i], it);
    }
  }

  MX Horzcat::join_primitives(std::vector<MX>::const_iterator& it) const {
    vector<MX> s(ndep());
    for (int i=0; i<s.size(); ++i) {
      s[i] = dep(i)->join_primitives(it);
    }
    return horzcat(s);
  }

  void Vertcat::split_primitives(const MX& x, std::vector<MX>::iterator& it) const {
    vector<MX> s = vertsplit(x, offset());
    for (int i=0; i<s.size(); ++i) {
      dep(i)->split_primitives(s[i], it);
    }
  }

  MX Vertcat::join_primitives(std::vector<MX>::const_iterator& it) const {
    vector<MX> s(ndep());
    for (int i=0; i<s.size(); ++i) {
      s[i] = dep(i)->join_primitives(it);
    }
    return vertcat(s);
  }

  void Diagcat::split_primitives(const MX& x, std::vector<MX>::iterator& it) const {
    std::pair<std::vector<int>, std::vector<int> > off = offset();
    vector<MX> s = diagsplit(x, off.first, off.second);
    for (int i=0; i<s.size(); ++i) {
      dep(i)->split_primitives(s[i], it);
    }
  }

  MX Diagcat::join_primitives(std::vector<MX>::const_iterator& it) const {
    vector<MX> s(ndep());
    for (int i=0; i<s.size(); ++i) {
      s[i] = dep(i)->join_primitives(it);
    }
    return diagcat(s);
  }

  bool Concat::has_duplicates() {
    bool has_duplicates = false;
    for (int i=0; i<ndep(); ++i) {
      has_duplicates = dep(i)->has_duplicates() || has_duplicates;
    }
    return has_duplicates;
  }

  void Concat::resetInput() {
    for (int i=0; i<ndep(); ++i) {
      dep(i)->resetInput();
    }
  }

  void Concat::primitives(std::vector<MX>::iterator& it) const {
    for (int i=0; i<ndep(); ++i) {
      dep(i)->primitives(it);
    }
  }

} // namespace casadi
