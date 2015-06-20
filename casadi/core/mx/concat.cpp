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

  void Concat::evalD(const cpv_double& input, const pv_double& output,
                         int* itmp, double* rtmp) {
    evalGen<double>(input, output, itmp, rtmp);
  }

  void Concat::evalSX(const cpv_SXElement& input, const pv_SXElement& output,
                          int* itmp, SXElement* rtmp) {
    evalGen<SXElement>(input, output, itmp, rtmp);
  }

  template<typename T>
  void Concat::evalGen(const std::vector<const T*>& input,
                       const std::vector<T*>& output, int* itmp, T* rtmp) {
    T* res = output[0];
    for (int i=0; i<ndep(); ++i) {
      const T* arg_i = input[i];
      copy(arg_i, arg_i+dep(i).nnz(), res);
      res += dep(i).nnz();
    }
  }

  void Concat::spFwd(const cpv_bvec_t& arg,
                     const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    bvec_t *res_ptr = res[0];
    for (int i=0; i<ndep(); ++i) {
      int n_i = dep(i).nnz();
      const bvec_t *arg_i_ptr = arg[i];
      copy(arg_i_ptr, arg_i_ptr+n_i, res_ptr);
      res_ptr += n_i;
    }
  }

  void Concat::spAdj(const pv_bvec_t& arg,
                     const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
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
    casadi_assert(x.size()>1);
    std::vector<Sparsity> sp(x.size());
    for (int i=0; i<x.size(); ++i)
      sp[i] = x[i].sparsity();
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

  void Diagcat::eval(const cpv_MX& arg, const pv_MX& res) {
    *res[0] = diagcat(getVector(arg, ndep()));
  }

  void Diagcat::evalFwd(const std::vector<cpv_MX>& fwdSeed, const std::vector<pv_MX>& fwdSens) {
    int nfwd = fwdSens.size();
    for (int d = 0; d<nfwd; ++d) {
      *fwdSens[d][0] = diagcat(getVector(fwdSeed[d], ndep()));
    }
  }

  void Diagcat::evalAdj(const std::vector<pv_MX>& adjSeed, const std::vector<pv_MX>& adjSens) {
    // Get offsets for each row and column
    vector<int> offset1(ndep()+1, 0);
    vector<int> offset2(ndep()+1, 0);
    for (int i=0; i<ndep(); ++i) {
      int ncol = dep(i).sparsity().size2();
      int nrow = dep(i).sparsity().size1();
      offset2[i+1] = offset2[i] + ncol;
      offset1[i+1] = offset1[i] + nrow;
    }

    // Adjoint sensitivities
    int nadj = adjSeed.size();
    for (int d=0; d<nadj; ++d) {
      MX& aseed = *adjSeed[d][0];
      vector<MX> s = diagsplit(aseed, offset1, offset2);
      aseed = MX();
      for (int i=0; i<ndep(); ++i) {
        adjSens[d][i]->addToSum(s[i]);
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

  void Horzcat::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "horzcat(";
    } else if (part==ndep()) {
      stream << ")";
    } else {
      stream << ", ";
    }
  }

  void Horzcat::eval(const cpv_MX& arg, const pv_MX& res) {
    *res[0] = horzcat(getVector(arg, ndep()));
  }

  void Horzcat::evalFwd(const std::vector<cpv_MX>& fwdSeed, const std::vector<pv_MX>& fwdSens) {
    int nfwd = fwdSens.size();
    for (int d = 0; d<nfwd; ++d) {
      *fwdSens[d][0] = horzcat(getVector(fwdSeed[d], ndep()));
    }
  }

  void Horzcat::evalAdj(const std::vector<pv_MX>& adjSeed, const std::vector<pv_MX>& adjSens) {
    // Get offsets for each column
    vector<int> col_offset(ndep()+1, 0);
    for (int i=0; i<ndep(); ++i) {
      int ncol = dep(i).sparsity().size2();
      col_offset[i+1] = col_offset[i] + ncol;
    }

    // Adjoint sensitivities
    int nadj = adjSeed.size();
    for (int d=0; d<nadj; ++d) {
      MX& aseed = *adjSeed[d][0];
      vector<MX> s = horzsplit(aseed, col_offset);
      aseed = MX();
      for (int i=0; i<ndep(); ++i) {
        adjSens[d][i]->addToSum(s[i]);
      }
    }
  }

  Vertcat::Vertcat(const std::vector<MX>& x) : Concat(x) {
    casadi_assert(x.size()>1);
    std::vector<Sparsity> sp(x.size());
    for (int i=0; i<x.size(); ++i)
      sp[i] = x[i].sparsity();
    setSparsity(vertcat(sp));
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

  void Vertcat::eval(const cpv_MX& arg, const pv_MX& res) {
    *res[0] = vertcat(getVector(arg, ndep()));
  }

  void Vertcat::evalFwd(const std::vector<cpv_MX>& fwdSeed, const std::vector<pv_MX>& fwdSens) {
    int nfwd = fwdSens.size();
    for (int d = 0; d<nfwd; ++d) {
      *fwdSens[d][0] = vertcat(getVector(fwdSeed[d], ndep()));
    }
  }

  void Vertcat::evalAdj(const std::vector<pv_MX>& adjSeed, const std::vector<pv_MX>& adjSens) {
    // Get offsets for each row
    vector<int> row_offset(ndep()+1, 0);
    for (int i=0; i<ndep(); ++i) {
      int nrow = dep(i).sparsity().size1();
      row_offset[i+1] = row_offset[i] + nrow;
    }

    // Adjoint sensitivities
    int nadj = adjSeed.size();
    for (int d=0; d<nadj; ++d) {
      MX& aseed = *adjSeed[d][0];
      vector<MX> s = vertsplit(aseed, row_offset);
      aseed = MX();
      for (int i=0; i<ndep(); ++i) {
        adjSens[d][i]->addToSum(s[i]);
      }
    }
  }


} // namespace casadi
