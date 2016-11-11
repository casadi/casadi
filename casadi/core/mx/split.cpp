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


#include "split.hpp"
#include "../std_vector_tools.hpp"
#include "../global_options.hpp"

using namespace std;

namespace casadi {

  Split::Split(const MX& x, const std::vector<int>& offset) : offset_(offset) {
    setDependencies(x);
    setSparsity(Sparsity::scalar());
  }

  Split::~Split() {
  }

  void Split::eval(const double** arg, double** res, int* iw, double* w, int mem) const {
    evalGen<double>(arg, res, iw, w, mem);
  }

  void Split::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    evalGen<SXElem>(arg, res, iw, w, mem);
  }

  template<typename T>
  void Split::evalGen(const T** arg, T** res, int* iw, T* w, int mem) const {
    // Number of derivatives
    int nx = offset_.size()-1;

    for (int i=0; i<nx; ++i) {
      int nz_first = offset_[i];
      int nz_last = offset_[i+1];
      if (res[i]!=0) {
        copy(arg[0]+nz_first, arg[0]+nz_last, res[i]);
      }
    }
  }

  void Split::sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    int nx = offset_.size()-1;
    for (int i=0; i<nx; ++i) {
      if (res[i]!=0) {
        const bvec_t *arg_ptr = arg[0] + offset_[i];
        int n_i = sparsity(i).nnz();
        bvec_t *res_i_ptr = res[i];
        for (int k=0; k<n_i; ++k) {
          *res_i_ptr++ = *arg_ptr++;
        }
      }
    }
  }

  void Split::sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    int nx = offset_.size()-1;
    for (int i=0; i<nx; ++i) {
      if (res[i]!=0) {
        bvec_t *arg_ptr = arg[0] + offset_[i];
        int n_i = sparsity(i).nnz();
        bvec_t *res_i_ptr = res[i];
        for (int k=0; k<n_i; ++k) {
          *arg_ptr++ |= *res_i_ptr;
          *res_i_ptr++ = 0;
        }
      }
    }
  }

  void Split::generate(CodeGenerator& g, const std::string& mem,
                       const std::vector<int>& arg, const std::vector<int>& res) const {
    int nx = nout();
    for (int i=0; i<nx; ++i) {
      int nz_first = offset_[i];
      int nz_last = offset_[i+1];
      int nz = nz_last-nz_first;
      if (res[i]>=0 && nz>0) { // if anything to assign
        if (nz==1) { // assign scalar
          g.body << "  " << g.workel(res[i]) << " = ";
          if (dep(0).nnz()==1) {
            // rhs is also scalar
            casadi_assert(nz_first==0);
            g.body << g.workel(arg[0]) << ";" << endl;
          } else {
            // rhs is an element in a vector
            g.body << g.work(arg[0], dep(0).nnz()) << "[" << nz_first << "];" << endl;
          }
        } else {
          // assign vector
          std::string r = g.work(arg[0], dep(0).nnz());
          if (nz_first!=0) r = r + "+" + g.to_string(nz_first);
          g.body << "  " << g.copy(r, nz, g.work(res[i], nnz(i))) << endl;
        }
      }
    }
  }

  Horzsplit::Horzsplit(const MX& x, const std::vector<int>& offset) : Split(x, offset) {

    // Split up the sparsity pattern
    output_sparsity_ = horzsplit(x.sparsity(), offset_);

    // Have offset_ refer to the nonzero offsets instead of column offsets
    offset_.resize(1);
    for (auto&& s : output_sparsity_) {
      offset_.push_back(offset_.back() + s.nnz());
    }
  }

  std::string Horzsplit::print(const std::vector<std::string>& arg) const {
    return "horzsplit(" + arg.at(0) + ")";
  }

  void Horzsplit::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) {
    // Get column offsets
    vector<int> col_offset;
    col_offset.reserve(offset_.size());
    col_offset.push_back(0);
    for (auto&& s : output_sparsity_) {
      col_offset.push_back(col_offset.back() + s.size2());
    }

    res = horzsplit(arg[0], col_offset);
  }

  void Horzsplit::evalFwd(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) {
    int nfwd = fsens.size();

    // Get column offsets
    vector<int> col_offset;
    col_offset.reserve(offset_.size());
    col_offset.push_back(0);
    for (auto&& s : output_sparsity_) {
      col_offset.push_back(col_offset.back() + s.size2());
    }

    // Non-differentiated output and forward sensitivities
    for (int d=0; d<nfwd; ++d) {
      fsens[d] = horzsplit(fseed[d][0], col_offset);
    }
  }

  void Horzsplit::evalAdj(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) {
    int nadj = aseed.size();

    // Get column offsets
    vector<int> col_offset;
    col_offset.reserve(offset_.size());
    col_offset.push_back(0);
    for (auto&& s : output_sparsity_) {
      col_offset.push_back(col_offset.back() + s.size2());
    }

    for (int d=0; d<nadj; ++d) {
      asens[d][0] += horzcat(aseed[d]);
    }
  }

  Diagsplit::Diagsplit(const MX& x,
    const std::vector<int>& offset1,
    const std::vector<int>& offset2) : Split(x, offset1) {

    // Split up the sparsity pattern
    output_sparsity_ = diagsplit(x.sparsity(), offset1, offset2);

    // Have offset_ refer to the nonzero offsets instead of column offsets
    offset_.resize(1);
    for (auto&& s : output_sparsity_) {
      offset_.push_back(offset_.back() + s.nnz());
    }

    casadi_assert_message(offset_.back()==x.nnz(),
      "DiagSplit:: the presence of nonzeros outside the diagonal blocks in unsupported.");
  }

  std::string Diagsplit::print(const std::vector<std::string>& arg) const {
    return "diagsplit(" + arg.at(0) + ")";
  }

  void Diagsplit::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) {
    // Get offsets
    vector<int> offset1;
    offset1.reserve(offset_.size());
    offset1.push_back(0);
    vector<int> offset2;
    offset2.reserve(offset_.size());
    offset2.push_back(0);
    for (auto&& s : output_sparsity_) {
      offset1.push_back(offset1.back() + s.size1());
      offset2.push_back(offset2.back() + s.size2());
    }

    res = diagsplit(arg[0], offset1, offset2);
  }

  void Diagsplit::evalFwd(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) {
    int nfwd = fsens.size();
    // Get offsets
    vector<int> offset1;
    offset1.reserve(offset_.size());
    offset1.push_back(0);
    vector<int> offset2;
    offset2.reserve(offset_.size());
    offset2.push_back(0);
    for (auto&& s : output_sparsity_) {
      offset1.push_back(offset1.back() + s.size1());
      offset2.push_back(offset2.back() + s.size2());
    }

    // Non-differentiated output and forward sensitivities
    for (int d=0; d<nfwd; ++d) {
      fsens[d] = diagsplit(fseed[d][0], offset1, offset2);
    }
  }

  void Diagsplit::evalAdj(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) {
    int nadj = asens.size();

    // Get offsets
    vector<int> offset1;
    offset1.reserve(offset_.size());
    offset1.push_back(0);
    vector<int> offset2;
    offset2.reserve(offset_.size());
    offset2.push_back(0);
    for (auto&& s : output_sparsity_) {
      offset1.push_back(offset1.back() + s.size1());
      offset2.push_back(offset2.back() + s.size2());
    }

    for (int d=0; d<nadj; ++d) {
      asens[d][0] += diagcat(aseed[d]);
    }
  }

  Vertsplit::Vertsplit(const MX& x, const std::vector<int>& offset) : Split(x, offset) {

    // Split up the sparsity pattern
    output_sparsity_ = vertsplit(x.sparsity(), offset_);

    // Have offset_ refer to the nonzero offsets instead of column offsets
    offset_.resize(1);
    for (auto&& s : output_sparsity_) {
      offset_.push_back(offset_.back() + s.nnz());
    }
  }

  std::string Vertsplit::print(const std::vector<std::string>& arg) const {
    return "vertsplit(" + arg.at(0) + ")";
  }

  void Vertsplit::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) {
    // Get row offsets
    vector<int> row_offset;
    row_offset.reserve(offset_.size());
    row_offset.push_back(0);
    for (auto&& s : output_sparsity_) {
      row_offset.push_back(row_offset.back() + s.size1());
    }

    res = vertsplit(arg[0], row_offset);
  }

  void Vertsplit::evalFwd(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) {
    int nfwd = fsens.size();

    // Get row offsets
    vector<int> row_offset;
    row_offset.reserve(offset_.size());
    row_offset.push_back(0);
    for (auto&& s : output_sparsity_) {
      row_offset.push_back(row_offset.back() + s.size1());
    }

    for (int d=0; d<nfwd; ++d) {
      fsens[d] = vertsplit(fseed[d][0], row_offset);
    }
  }

  void Vertsplit::evalAdj(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) {
    int nadj = aseed.size();

    // Get row offsets
    vector<int> row_offset;
    row_offset.reserve(offset_.size());
    row_offset.push_back(0);
    for (auto&& s : output_sparsity_) {
      row_offset.push_back(row_offset.back() + s.size1());
    }

    for (int d=0; d<nadj; ++d) {
      asens[d][0] += vertcat(aseed[d]);
    }
  }

  MX Horzsplit::getHorzcat(const std::vector<MX>& x) const {
    // Check x length
    if (x.size()!=nout()) {
      return MXNode::getHorzcat(x);
    }

    // Check x content
    for (int i=0; i<x.size(); ++i) {
      if (!(x[i]->isOutputNode() && x[i]->getFunctionOutput()==i && x[i]->dep().get()==this)) {
        return MXNode::getHorzcat(x);
      }
    }

    // OK if reached this point
    return dep();
  }

  MX Vertsplit::getVertcat(const std::vector<MX>& x) const {
    // Check x length
    if (x.size()!=nout()) {
      return MXNode::getVertcat(x);
    }

    // Check x content
    for (int i=0; i<x.size(); ++i) {
      if (!(x[i]->isOutputNode() && x[i]->getFunctionOutput()==i && x[i]->dep().get()==this)) {
        return MXNode::getVertcat(x);
      }
    }

    // OK if reached this point
    return dep();
  }

  MX Diagsplit::get_diagcat(const std::vector<MX>& x) const {
    // Check x length
    if (x.size()!=nout()) {
      return MXNode::get_diagcat(x);
    }

    // Check x content
    for (int i=0; i<x.size(); ++i) {
      if (!(x[i]->isOutputNode() && x[i]->getFunctionOutput()==i && x[i]->dep().get()==this)) {
        return MXNode::get_diagcat(x);
      }
    }

    // OK if reached this point
    return dep();
  }

} // namespace casadi
