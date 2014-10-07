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
#include "../matrix/sparsity_tools.hpp"

using namespace std;

namespace casadi {

  Concat::Concat(const vector<MX>& x) {
    setDependencies(x);
  }

  Concat::~Concat() {
  }

  void Concat::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                         std::vector<double>& rtmp) {
    evaluateGen<double, DMatrixPtrV, DMatrixPtrVV>(input, output, itmp, rtmp);
  }

  void Concat::evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                          std::vector<SXElement>& rtmp) {
    evaluateGen<SXElement, SXPtrV, SXPtrVV>(input, output, itmp, rtmp);
  }

  template<typename T, typename MatV, typename MatVV>
  void Concat::evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp,
                           std::vector<T>& rtmp) {
    typename vector<T>::iterator res_it = output[0]->data().begin();
    for (int i=0; i<input.size(); ++i) {
      const vector<T>& arg_i = input[i]->data();
      copy(arg_i.begin(), arg_i.end(), res_it);
      res_it += arg_i.size();
    }
  }

  void Concat::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd) {
    bvec_t *res_ptr = get_bvec_t(output[0]->data());
    for (int i=0; i<input.size(); ++i) {
      vector<double>& arg_i = input[i]->data();
      bvec_t *arg_i_ptr = get_bvec_t(arg_i);
      if (fwd) {
        copy(arg_i_ptr, arg_i_ptr+arg_i.size(), res_ptr);
        res_ptr += arg_i.size();
      } else {
        for (int k=0; k<arg_i.size(); ++k) {
          *arg_i_ptr++ |= *res_ptr;
          *res_ptr++ = 0;
        }
      }
    }
  }

  void Concat::generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                 const std::vector<std::string>& res, CodeGenerator& gen) const {
    int nz_offset = 0;
    for (int i=0; i<arg.size(); ++i) {
      int nz = dep(i).size();
      stream << "  for (i=0; i<" << nz << "; ++i) " << res.front() << "[i+" << nz_offset
             << "] = " << arg.at(i) << "[i];" << endl;
      nz_offset += nz;
    }
    casadi_assert(nz_offset == size());
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
      end += dep(i).size();
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

    setSparsity(blkdiag(sp));
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

  void Diagcat::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                           MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                           bool output_given) {
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Non-differentiated output
    if (!output_given) {
      *output[0] = diagcat(getVector(input));
    }

    // Forward sensitivities
    for (int d = 0; d<nfwd; ++d) {
      *fwdSens[d][0] = diagcat(getVector(fwdSeed[d]));
    }

    // Quick return?
    if (nadj==0) return;

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

  void Horzcat::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                           MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                           bool output_given) {
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Non-differentiated output
    if (!output_given) {
      *output[0] = horzcat(getVector(input));
    }

    // Forward sensitivities
    for (int d = 0; d<nfwd; ++d) {
      *fwdSens[d][0] = horzcat(getVector(fwdSeed[d]));
    }

    // Quick return?
    if (nadj==0) return;

    // Get offsets for each column
    vector<int> col_offset(ndep()+1, 0);
    for (int i=0; i<ndep(); ++i) {
      int ncol = dep(i).sparsity().size2();
      col_offset[i+1] = col_offset[i] + ncol;
    }

    // Adjoint sensitivities
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

  void Vertcat::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                           MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                           bool output_given) {
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Non-differentiated output
    if (!output_given) {
      *output[0] = vertcat(getVector(input));
    }

    // Forward sensitivities
    for (int d = 0; d<nfwd; ++d) {
      *fwdSens[d][0] = vertcat(getVector(fwdSeed[d]));
    }

    // Quick return?
    if (nadj==0) return;

    // Get offsets for each row
    vector<int> row_offset(ndep()+1, 0);
    for (int i=0; i<ndep(); ++i) {
      int nrow = dep(i).sparsity().size1();
      row_offset[i+1] = row_offset[i] + nrow;
    }

    // Adjoint sensitivities
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
