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
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../function/sx_function.hpp"

using namespace std;

namespace casadi {

  Concat::Concat(const vector<MX>& x) {
    setDependencies(x);
  }

  Concat::~Concat() {
  }

  void Concat::evalD(cp_double* input, p_double* output,
                         int* itmp, double* rtmp) {
    evalGen<double>(input, output, itmp, rtmp);
  }

  void Concat::evalSX(cp_SXElement* input, p_SXElement* output,
                          int* itmp, SXElement* rtmp) {
    evalGen<SXElement>(input, output, itmp, rtmp);
  }

  template<typename T>
  void Concat::evalGen(const T* const* arg, T* const* res,
                       int* itmp, T* rtmp) {
    T* r = res[0];
    for (int i=0; i<ndep(); ++i) {
      int n = dep(i).nnz();
      copy(arg[i], arg[i]+n, r);
      r += n;
    }
  }

  void Concat::spFwd(cp_bvec_t* arg,
                     p_bvec_t* res, int* itmp, bvec_t* rtmp) {
    bvec_t *res_ptr = res[0];
    for (int i=0; i<ndep(); ++i) {
      int n_i = dep(i).nnz();
      const bvec_t *arg_i_ptr = arg[i];
      copy(arg_i_ptr, arg_i_ptr+n_i, res_ptr);
      res_ptr += n_i;
    }
  }

  void Concat::spAdj(p_bvec_t* arg,
                     p_bvec_t* res, int* itmp, bvec_t* rtmp) {
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

  void Concat::generate(std::ostream &stream, const std::vector<int>& arg,
                        const std::vector<int>& res, CodeGenerator& gen) const {
    for (int i=0; i<arg.size(); ++i) {
      int nz = dep(i).nnz();
      stream << "  for (i=0, ";
      if (i==0) stream << "rr=" << gen.work(res[0]) << ", ";
      stream << "cs=" << gen.work(arg[i]) << "; i<" << nz << "; ++i) *rr++ = *cs++;" << endl;
    }
  }

  MX Concat::getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const {
    // Get the first nonnegative nz
    int nz_test = -1;
    for (vector<int>::const_iterator i=nz.begin(); i!=nz.end(); ++i) {
      if (*i>=0) {
        nz_test = *i;
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
    for (vector<int>::const_iterator j=nz.begin(); j!=nz.end(); ++j) {
      if (*j>=0 && (*j < begin || *j >= end)) {

        // Fallback to the base class
        return MXNode::getGetNonzeros(sp, nz);
      }
    }

    // All nz refer to the same dependency, update the nonzero indices
    if (begin==0) {
      return dep(i)->getGetNonzeros(sp, nz);
    } else {
      vector<int> nz_new(nz);
      for (vector<int>::iterator j=nz_new.begin(); j!=nz_new.end(); ++j) {
        if (*j>=0) *j -= begin;
      }
      return dep(i)->getGetNonzeros(sp, nz_new);
    }
  }


  Diagcat::Diagcat(const std::vector<MX>& x) : Concat(x) {
    // Construct the sparsity
    casadi_assert(!x.empty());
    std::vector<Sparsity> sp;
    for (int i=0;i<x.size(); ++i) {
      sp.push_back(x[i].sparsity());
    }

    setSparsity(diagcat(sp));
  }

  void Diagcat::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "diagcat(";
    } else if (part==ndep()) {
      stream << ")";
    } else {
      stream << ", ";
    }
  }

  void Diagcat::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = diagcat(arg);
  }

  void Diagcat::evalFwd(const std::vector<std::vector<MX> >& fseed,
                        std::vector<std::vector<MX> >& fsens) {
    int nfwd = fsens.size();
    for (int d = 0; d<nfwd; ++d) {
      fsens[d][0] = diagcat(fseed[d]);
    }
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
    std::pair<std::vector<int>, std::vector<int> > off = offset();

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
    // Construct the sparsity
    casadi_assert(!x.empty());
    Sparsity sp = x.front().sparsity();
    for (vector<MX>::const_iterator i=x.begin()+1; i!=x.end(); ++i) {
      sp.appendColumns(i->sparsity());
    }

    setSparsity(sp);
  }

  void Horzcat::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "horzcat(";
    } else if (part==ndep()) {
      stream << ")";
    } else {
      stream << ", ";
    }
  }

  void Horzcat::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
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
    // Construct the sparsity
    casadi_assert(!x.empty());
    Sparsity sp = x.front().sparsity();
    for (vector<MX>::const_iterator i=x.begin()+1; i!=x.end(); ++i) {
      sp.append(i->sparsity());
    }

    setSparsity(sp);
  }

  void Vertcat::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "vertcat(";
    } else if (part==ndep()) {
      stream << ")";
    } else {
      stream << ", ";
    }
  }

  void Vertcat::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
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

  bool Concat::isValidInput() const {
    for (int i=0; i<ndep(); ++i) {
      if (!dep(i)->isValidInput()) return false;
    }
    return true;
  }

  int Concat::numPrimitives() const {
    int nprim = 0;
    for (int i=0; i<ndep(); ++i) {
      nprim +=  dep(i)->numPrimitives();
    }
    return nprim;
  }

  void Horzcat::splitPrimitives(const MX& x, std::vector<MX>::iterator& it) const {
    vector<MX> s = horzsplit(x, offset());
    for (int i=0; i<s.size(); ++i) {
      dep(i)->splitPrimitives(s[i], it);
    }
  }

  MX Horzcat::joinPrimitives(std::vector<MX>::const_iterator& it) const {
    vector<MX> s(ndep());
    for (int i=0; i<s.size(); ++i) {
      s[i] = dep(i)->joinPrimitives(it);
    }
    return horzcat(s);
  }

  void Vertcat::splitPrimitives(const MX& x, std::vector<MX>::iterator& it) const {
    vector<MX> s = vertsplit(x, offset());
    for (int i=0; i<s.size(); ++i) {
      dep(i)->splitPrimitives(s[i], it);
    }
  }

  MX Vertcat::joinPrimitives(std::vector<MX>::const_iterator& it) const {
    vector<MX> s(ndep());
    for (int i=0; i<s.size(); ++i) {
      s[i] = dep(i)->joinPrimitives(it);
    }
    return vertcat(s);
  }

  void Diagcat::splitPrimitives(const MX& x, std::vector<MX>::iterator& it) const {
    std::pair<std::vector<int>, std::vector<int> > off = offset();
    vector<MX> s = diagsplit(x, off.first, off.second);
    for (int i=0; i<s.size(); ++i) {
      dep(i)->splitPrimitives(s[i], it);
    }
  }

  MX Diagcat::joinPrimitives(std::vector<MX>::const_iterator& it) const {
    vector<MX> s(ndep());
    for (int i=0; i<s.size(); ++i) {
      s[i] = dep(i)->joinPrimitives(it);
    }
    return diagcat(s);
  }

  bool Concat::hasDuplicates() {
    bool has_duplicates = false;
    for (int i=0; i<ndep(); ++i) {
      has_duplicates = dep(i)->hasDuplicates() || has_duplicates;
    }
    return has_duplicates;
  }

  void Concat::resetInput() {
    for (int i=0; i<ndep(); ++i) {
      dep(i)->resetInput();
    }
  }

  void Concat::getPrimitives(std::vector<MX>::iterator& it) const {
    for (int i=0; i<ndep(); ++i) {
      dep(i)->getPrimitives(it);
    }
  }

} // namespace casadi
