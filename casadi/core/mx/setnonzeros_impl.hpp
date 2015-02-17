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


#ifndef CASADI_SETNONZEROS_IMPL_HPP
#define CASADI_SETNONZEROS_IMPL_HPP

#include "setnonzeros.hpp"
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../function/sx_function.hpp"

/// \cond INTERNAL

using namespace std;

namespace casadi {

  template<bool Add>
  SetNonzeros<Add>::SetNonzeros(const MX& y, const MX& x) {
    this->setSparsity(y.sparsity());
    this->setDependencies(y, x);
  }

  template<bool Add>
  SetNonzerosVector<Add>::SetNonzerosVector(const MX& y, const MX& x, const std::vector<int>& nz) :
      SetNonzeros<Add>(y, x), nz_(nz) {
    // Ignore duplicate assignments
    if (!Add) {
      vector<bool> already_set(this->nnz(), false);
      for (vector<int>::reverse_iterator i=nz_.rbegin(); i!=nz_.rend(); ++i) {
        if (*i>=0) {
          if (already_set[*i]) {
            *i = -1;
          } else {
            already_set[*i] = true;
          }
        }
      }
    }
  }

  template<bool Add>
  SetNonzeros<Add>:: ~SetNonzeros() {
  }

  template<bool Add>
  void SetNonzeros<Add>::eval(const cpv_MX& input, const pv_MX& output) {
    // Get all the nonzeros
    vector<int> nz = getAll();

    // Output sparsity
    const Sparsity &osp = sparsity();
    const int* orow = osp.row();
    vector<int> ocol = osp.getCol();

    // Input sparsity (first input same as output)
    const Sparsity &isp = dep(1).sparsity();
    const int* irow = isp.row();
    vector<int> icol = isp.getCol();

    // We next need to resort the assignment vector by outputs instead of inputs
    // Start by counting the number of output nonzeros corresponding to each input nonzero
    vector<int> onz_count(osp.nnz()+2, 0);
    for (vector<int>::const_iterator it=nz.begin(); it!=nz.end(); ++it) {
      onz_count[*it+2]++;
    }

    // Cumsum to get index offset for output nonzero
    for (int i=0; i<onz_count.size()-1; ++i) {
      onz_count[i+1] += onz_count[i];
    }

    // Get the order of assignments
    vector<int> nz_order(nz.size());
    for (int k=0; k<nz.size(); ++k) {
      // Save the new index
      nz_order[onz_count[1+nz[k]]++] = k;
    }

    // Find out which elements are being set
    vector<int>& with_duplicates = onz_count; // Reuse memory
    onz_count.resize(nz.size());
    for (int k=0; k<nz.size(); ++k) {
      // Get output nonzero
      int onz_k = nz[nz_order[k]];

      // Get element (note: may contain duplicates)
      if (onz_k>=0) {
        with_duplicates[k] = ocol[onz_k]*osp.size1() + orow[onz_k];
      } else {
        with_duplicates[k] = -1;
      }
    }

    // Get all output elements (this time without duplicates)
    vector<int> el_output;
    osp.find(el_output);

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<int> r_colind, r_row, r_nz, r_ind;

    // Get references to arguments and results
    const MX& arg = *input[1];
    const MX& arg0 = *input[0];
    MX& res = *output[0];
    if (&arg0 != &res) {
      res = arg0;
    }

    // Entries in res with elements zero'ed out
    if (!Add) {

      // Get the nz locations in res corresponding to the output sparsity pattern
      r_nz.resize(with_duplicates.size());
      copy(with_duplicates.begin(), with_duplicates.end(), r_nz.begin());
      res.sparsity().getNZ(r_nz);

      // Zero out the corresponding entries
      res = MX::zeros(isp)->getSetNonzeros(res, r_nz);
    }

    // Get the nz locations of the elements in arg corresponding to the argument sparsity pattern
    arg.sparsity().find(r_nz);
    isp.getNZ(r_nz);

    // Filter out ignored entries and check if there is anything to add at all
    bool elements_to_add = false;
    for (vector<int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
      if (*k>=0) {
        if (nz[*k]>=0) {
          elements_to_add = true;
        } else {
          *k = -1;
        }
      }
    }

    // Quick continue of no elements to set/add
    if (!elements_to_add) return;

    // Get the nz locations in the argument corresponding to the inputs
    r_ind.resize(el_output.size());
    copy(el_output.begin(), el_output.end(), r_ind.begin());
    res.sparsity().getNZ(r_ind);

    // Enlarge the sparsity pattern of the arguments if not all assignments fit
    for (vector<int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
      if (*k>=0 && nz[*k]>=0 && r_ind[nz[*k]]<0) {

        // Create a new pattern which includes both the the previous seed
        // and the addition/assignment
        Sparsity sp = res.sparsity().patternUnion(osp);
        res = res->getSetSparse(sp);

        // Recalculate the nz locations in the arguments corresponding to the inputs
        copy(el_output.begin(), el_output.end(), r_ind.begin());
        res.sparsity().getNZ(r_ind);

        break;
      }
    }

    // Have r_nz point to locations in the result instead of the output
    for (vector<int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
      if (*k>=0) {
        *k = r_ind[nz[*k]];
      }
    }

    // Add to the element to the sensitivity, if any
    res = arg->getAddNonzeros(res, r_nz);
  }

  template<bool Add>
  void SetNonzeros<Add>::evalFwd(const std::vector<cpv_MX>& fwdSeed,
                                 const std::vector<pv_MX>& fwdSens) {
    // Get all the nonzeros
    vector<int> nz = getAll();

    // Number of derivative directions
    int nfwd = fwdSens.size();

    // Output sparsity
    const Sparsity &osp = sparsity();
    const int* orow = osp.row();
    vector<int> ocol = osp.getCol();

    // Input sparsity (first input same as output)
    const Sparsity &isp = dep(1).sparsity();
    const int* irow = isp.row();
    vector<int> icol = isp.getCol();

    // We next need to resort the assignment vector by outputs instead of inputs
    // Start by counting the number of output nonzeros corresponding to each input nonzero
    vector<int> onz_count(osp.nnz()+2, 0);
    for (vector<int>::const_iterator it=nz.begin(); it!=nz.end(); ++it) {
      onz_count[*it+2]++;
    }

    // Cumsum to get index offset for output nonzero
    for (int i=0; i<onz_count.size()-1; ++i) {
      onz_count[i+1] += onz_count[i];
    }

    // Get the order of assignments
    vector<int> nz_order(nz.size());
    for (int k=0; k<nz.size(); ++k) {
      // Save the new index
      nz_order[onz_count[1+nz[k]]++] = k;
    }

    // Find out which elements are being set
    vector<int>& with_duplicates = onz_count; // Reuse memory
    onz_count.resize(nz.size());
    for (int k=0; k<nz.size(); ++k) {
      // Get output nonzero
      int onz_k = nz[nz_order[k]];

      // Get element (note: may contain duplicates)
      if (onz_k>=0) {
        with_duplicates[k] = ocol[onz_k]*osp.size1() + orow[onz_k];
      } else {
        with_duplicates[k] = -1;
      }
    }

    // Get all output elements (this time without duplicates)
    vector<int> el_output;
    osp.find(el_output);

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<int> r_colind, r_row, r_nz, r_ind;

    // Nondifferentiated function and forward sensitivities
    for (int d=0; d<nfwd; ++d) {

      // Get references to arguments and results
      const MX& arg = *fwdSeed[d][1];
      const MX& arg0 = *fwdSeed[d][0];
      MX& res = *fwdSens[d][0];
      if (&arg0 != &res) {
        res = arg0;
      }

      // Entries in res with elements zero'ed out
      if (!Add) {

        // Get the nz locations in res corresponding to the output sparsity pattern
        r_nz.resize(with_duplicates.size());
        copy(with_duplicates.begin(), with_duplicates.end(), r_nz.begin());
        res.sparsity().getNZ(r_nz);

        // Zero out the corresponding entries
        res = MX::zeros(isp)->getSetNonzeros(res, r_nz);
      }

      // Get the nz locations of the elements in arg corresponding to the argument sparsity pattern
      arg.sparsity().find(r_nz);
      isp.getNZ(r_nz);

      // Filter out ignored entries and check if there is anything to add at all
      bool elements_to_add = false;
      for (vector<int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
        if (*k>=0) {
          if (nz[*k]>=0) {
            elements_to_add = true;
          } else {
            *k = -1;
          }
        }
      }

      // Quick continue of no elements to set/add
      if (!elements_to_add) continue;

      // Get the nz locations in the argument corresponding to the inputs
      r_ind.resize(el_output.size());
      copy(el_output.begin(), el_output.end(), r_ind.begin());
      res.sparsity().getNZ(r_ind);

      // Enlarge the sparsity pattern of the arguments if not all assignments fit
      for (vector<int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
        if (*k>=0 && nz[*k]>=0 && r_ind[nz[*k]]<0) {

          // Create a new pattern which includes both the the previous seed
          // and the addition/assignment
          Sparsity sp = res.sparsity().patternUnion(osp);
          res = res->getSetSparse(sp);

          // Recalculate the nz locations in the arguments corresponding to the inputs
          copy(el_output.begin(), el_output.end(), r_ind.begin());
          res.sparsity().getNZ(r_ind);

          break;
        }
      }

      // Have r_nz point to locations in the result instead of the output
      for (vector<int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
        if (*k>=0) {
          *k = r_ind[nz[*k]];
        }
      }

      // Add to the element to the sensitivity, if any
      res = arg->getAddNonzeros(res, r_nz);
    }
  }

  template<bool Add>
  void SetNonzeros<Add>::evalAdj(const std::vector<pv_MX>& adjSeed,
                                 const std::vector<pv_MX>& adjSens) {
    // Get all the nonzeros
    vector<int> nz = getAll();

    // Number of derivative directions
    int nadj = adjSeed.size();

    // Output sparsity
    const Sparsity &osp = sparsity();
    const int* orow = osp.row();
    vector<int> ocol = osp.getCol();

    // Input sparsity (first input same as output)
    const Sparsity &isp = dep(1).sparsity();
    const int* irow = isp.row();
    vector<int> icol = isp.getCol();

    // We next need to resort the assignment vector by outputs instead of inputs
    // Start by counting the number of output nonzeros corresponding to each input nonzero
    vector<int> onz_count(osp.nnz()+2, 0);
    for (vector<int>::const_iterator it=nz.begin(); it!=nz.end(); ++it) {
      onz_count[*it+2]++;
    }

    // Cumsum to get index offset for output nonzero
    for (int i=0; i<onz_count.size()-1; ++i) {
      onz_count[i+1] += onz_count[i];
    }

    // Get the order of assignments
    vector<int> nz_order(nz.size());
    for (int k=0; k<nz.size(); ++k) {
      // Save the new index
      nz_order[onz_count[1+nz[k]]++] = k;
    }

    // Find out which elements are being set
    vector<int>& with_duplicates = onz_count; // Reuse memory
    onz_count.resize(nz.size());
    for (int k=0; k<nz.size(); ++k) {
      // Get output nonzero
      int onz_k = nz[nz_order[k]];

      // Get element (note: may contain duplicates)
      if (onz_k>=0) {
        with_duplicates[k] = ocol[onz_k]*osp.size1() + orow[onz_k];
      } else {
        with_duplicates[k] = -1;
      }
    }

    // Get all output elements (this time without duplicates)
    vector<int> el_output;
    osp.find(el_output);

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<int> r_colind, r_row, r_nz, r_ind;

    for (int d=0; d<nadj; ++d) {

      // Get an owning references to the seeds and sensitivities
      // and clear the seeds for the next run
      MX& aseed = *adjSeed[d][0];
      MX& asens0 = *adjSens[d][0];
      MX& asens = *adjSens[d][1];

      // Get the matching nonzeros
      r_ind.resize(el_output.size());
      copy(el_output.begin(), el_output.end(), r_ind.begin());
      aseed.sparsity().getNZ(r_ind);

      // Sparsity pattern for the result
      r_colind.resize(isp.size2()+1); // Col count
      fill(r_colind.begin(), r_colind.end(), 0);
      r_row.clear();

      // Perform the assignments
      r_nz.clear();
      for (int k=0; k<nz.size(); ++k) {

        // Get the corresponding nonzero for the input
        int el = nz[k];

        // Skip if zero assignment
        if (el==-1) continue;

        // Get the corresponding nonzero in the argument
        int el_arg = r_ind[el];

        // Skip if no argument
        if (el_arg==-1) continue;

        // Save the assignment
        r_nz.push_back(el_arg);

        // Get the corresponding element
        int i=icol[k], j=irow[k];

        // Add to sparsity pattern
        r_row.push_back(j);
        r_colind[1+i]++;
      }

      // col count -> col offset
      for (int i=1; i<r_colind.size(); ++i) r_colind[i] += r_colind[i-1];

      // If anything to set/add
      if (!r_nz.empty()) {
        // Create a sparsity pattern from vectors
        Sparsity f_sp(isp.size1(), isp.size2(), r_colind, r_row);
        asens.addToSum(aseed->getGetNonzeros(f_sp, r_nz));
        if (!Add) {
          aseed = MX::zeros(f_sp)->getSetNonzeros(aseed, r_nz);
        }
      }

      if (&aseed != &asens0) {
        asens0.addToSum(aseed);
        aseed = MX();
      }
    }
  }

  template<bool Add>
  void SetNonzerosVector<Add>::evalD(const cpv_double& input, const pv_double& output,
                                         int* itmp, double* rtmp) {
    evalGen<double>(input, output, itmp, rtmp);
  }

  template<bool Add>
  void SetNonzerosVector<Add>::evalSX(const cpv_SXElement& input, const pv_SXElement& output,
                                          int* itmp, SXElement* rtmp) {
    evalGen<SXElement>(input, output, itmp, rtmp);
  }

  template<bool Add>
  template<typename T>
  void SetNonzerosVector<Add>::evalGen(const std::vector<const T*>& input,
                                       const std::vector<T*>& output, int* itmp, T* rtmp) {
    const T* idata0 = input[0];
    const T* idata = input[1];
    T* odata = output[0];
    if (idata0 != odata) {
      copy(idata0, idata0+this->dep(0).nnz(), odata);
    }
    for (vector<int>::const_iterator k=this->nz_.begin(); k!=this->nz_.end(); ++k, ++idata) {
      if (Add) {
        if (*k>=0) odata[*k] += *idata;
      } else {
        if (*k>=0) odata[*k] = *idata;
      }
    }
  }

  template<bool Add>
  void SetNonzerosSlice<Add>::evalD(const cpv_double& input, const pv_double& output,
                                    int* itmp, double* rtmp) {
    evalGen<double>(input, output, itmp, rtmp);
  }

  template<bool Add>
  void SetNonzerosSlice<Add>::evalSX(const cpv_SXElement& input, const pv_SXElement& output,
                                         int* itmp, SXElement* rtmp) {
    evalGen<SXElement>(input, output, itmp, rtmp);
  }

  template<bool Add>
  template<typename T>
  void SetNonzerosSlice<Add>::evalGen(const std::vector<const T*>& input,
                                      const std::vector<T*>& output, int* itmp, T* rtmp) {
    const T* idata0 = input[0];
    const T* idata = input[1];
    T* odata = output[0];
    if (idata0 != odata) {
      copy(idata0, idata0+this->dep(0).nnz(), odata);
    }
    T* odata_stop = odata + s_.stop_;
    for (odata += s_.start_; odata != odata_stop; odata += s_.step_) {
      if (Add) {
        *odata += *idata++;
      } else {
        *odata = *idata++;
      }
    }
  }

  template<bool Add>
  void SetNonzerosSlice2<Add>::evalD(const cpv_double& input, const pv_double& output,
                                         int* itmp, double* rtmp) {
    evalGen<double>(input, output, itmp, rtmp);
  }

  template<bool Add>
  void SetNonzerosSlice2<Add>::evalSX(const cpv_SXElement& input, const pv_SXElement& output,
                                          int* itmp, SXElement* rtmp) {
    evalGen<SXElement>(input, output, itmp, rtmp);
  }

  template<bool Add>
  template<typename T>
  void SetNonzerosSlice2<Add>::evalGen(const std::vector<const T*>& input,
                                       const std::vector<T*>& output, int* itmp, T* rtmp) {
    const T* idata0 = input[0];
    const T* idata = input[1];
    T* odata = output[0];
    if (idata0 != odata) {
      copy(idata0, idata0 + this->dep(0).nnz(), odata);
    }
    T* outer_stop = odata + outer_.stop_;
    T* outer = odata + outer_.start_;
    for (; outer != outer_stop; outer += outer_.step_) {
      for (T* inner = outer+inner_.start_;
          inner != outer+inner_.stop_;
          inner += inner_.step_) {
        if (Add) {
          *inner += *idata++;
        } else {
          *inner = *idata++;
        }
      }
    }
  }

  template<bool Add>
  void SetNonzerosVector<Add>::
  spFwd(const cpv_bvec_t& arg,
        const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    const bvec_t *a0 = arg[0];
    const bvec_t *a = arg[1];
    bvec_t *r = res[0];
    int n = this->nnz();

    // Propagate sparsity
    if (r != a0) copy(a0, a0+n, r);
    for (vector<int>::const_iterator k=this->nz_.begin(); k!=this->nz_.end(); ++k, ++a) {
      if (Add) {
        if (*k>=0) r[*k] |= *a;
      } else {
        if (*k>=0) r[*k] = *a;
      }
    }
  }

  template<bool Add>
  void SetNonzerosVector<Add>::
  spAdj(const pv_bvec_t& arg,
        const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    bvec_t *a = arg[1];
    bvec_t *r = res[0];
    for (vector<int>::const_iterator k=this->nz_.begin(); k!=this->nz_.end(); ++k, ++a) {
      if (*k>=0) {
        *a |= r[*k];
        if (!Add) {
          r[*k] = 0;
        }
      }
    }
    MXNode::copyAdj(arg[0], r, this->nnz());
  }

  template<bool Add>
  void SetNonzerosSlice<Add>::
  spFwd(const cpv_bvec_t& arg,
        const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    const bvec_t *a0 = arg[0];
    const bvec_t *a = arg[1];
    bvec_t *r = res[0];
    int n = this->nnz();

    // Propagate sparsity
    if (r != a0) copy(a0, a0+n, r);
    for (int k=s_.start_; k!=s_.stop_; k+=s_.step_) {
      if (Add) {
        r[k] |= *a++;
      } else {
        r[k] = *a++;
      }
    }
  }

  template<bool Add>
  void SetNonzerosSlice<Add>::
  spAdj(const pv_bvec_t& arg,
        const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    bvec_t *a = arg[1];
    bvec_t *r = res[0];
    for (int k=s_.start_; k!=s_.stop_; k+=s_.step_) {
      *a++ |= r[k];
      if (!Add) {
        r[k] = 0;
      }
    }
    MXNode::copyAdj(arg[0], r, this->nnz());
  }

  template<bool Add>
  void SetNonzerosSlice2<Add>::
  spFwd(const cpv_bvec_t& arg,
        const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    const bvec_t *a0 = arg[0];
    const bvec_t *a = arg[1];
    bvec_t *r = res[0];
    int n = this->nnz();

    // Propagate sparsity
    if (r != a0) copy(a0, a0+n, r);
    for (int k1=outer_.start_; k1!=outer_.stop_; k1+=outer_.step_) {
      for (int k2=k1+inner_.start_; k2!=k1+inner_.stop_; k2+=inner_.step_) {
        if (Add) {
          r[k2] |= *a++;
        } else {
          r[k2] = *a++;
        }
      }
    }
  }

  template<bool Add>
  void SetNonzerosSlice2<Add>::
  spAdj(const pv_bvec_t& arg,
        const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    bvec_t *a = arg[1];
    bvec_t *r = res[0];
    for (int k1=outer_.start_; k1!=outer_.stop_; k1+=outer_.step_) {
      for (int k2=k1+inner_.start_; k2!=k1+inner_.stop_; k2+=inner_.step_) {
        *a++ |= r[k2];
        if (!Add) {
          r[k2] = 0;
        }
      }
    }
    MXNode::copyAdj(arg[0], r, this->nnz());
  }

  template<bool Add>
  void SetNonzerosVector<Add>::printPart(std::ostream &stream, int part) const {
    switch (part) {
    case 0: stream << "(";           break;
    case 1: stream << this->nz_ << (Add ? " += " : " = ") ; break;
    case 2: stream << ")";           break;
    }
  }

  template<bool Add>
  void SetNonzerosSlice<Add>::printPart(std::ostream &stream, int part) const {
    switch (part) {
    case 0: stream << "(";           break;
    case 1: stream << "[" << s_ << "]" << (Add ? " += " : " = "); break;
    case 2: stream << ")";           break;
    }
  }

  template<bool Add>
  void SetNonzerosSlice2<Add>::printPart(std::ostream &stream, int part) const {
    switch (part) {
    case 0: stream << "(";           break;
    case 1: stream << "[" << outer_ << ";" << inner_ << "]" << (Add ? " += " : " = "); break;
    case 2: stream << ")";           break;
    }
  }

  template<bool Add>
  Matrix<int> SetNonzeros<Add>::mapping() const {
    vector<int> nz = getAll();
    return Matrix<int>(this->dep(1).sparsity(), nz, false);
  }

  template<bool Add>
  bool SetNonzerosVector<Add>::zz_isEqual(const MXNode* node, int depth) const {
    // Check dependencies
    if (!this->sameOpAndDeps(node, depth)) return false;

    // Check if same node
    const SetNonzerosVector<Add>* n = dynamic_cast<const SetNonzerosVector<Add>*>(node);
    if (n==0) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (this->nz_.size()!=n->nz_.size()) return false;
    if (!std::equal(this->nz_.begin(), this->nz_.end(), n->nz_.begin())) return false;

    return true;
  }

  template<bool Add>
  bool SetNonzerosSlice<Add>::zz_isEqual(const MXNode* node, int depth) const {
    // Check dependencies
    if (!this->sameOpAndDeps(node, depth)) return false;

    // Check if same node
    const SetNonzerosSlice<Add>* n = dynamic_cast<const SetNonzerosSlice<Add>*>(node);
    if (n==0) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (this->s_ != n->s_) return false;

    return true;
  }

  template<bool Add>
  bool SetNonzerosSlice2<Add>::zz_isEqual(const MXNode* node, int depth) const {
    // Check dependencies
    if (!this->sameOpAndDeps(node, depth)) return false;

    // Check if same node
    const SetNonzerosSlice2<Add>* n = dynamic_cast<const SetNonzerosSlice2<Add>*>(node);
    if (n==0) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (this->inner_ != n->inner_ || this->outer_!=n->outer_) return false;

    return true;
  }

  template<bool Add>
  bool SetNonzerosSlice<Add>::isAssignment() const {
    // Check sparsity
    if (!(this->sparsity() == this->dep(1).sparsity()))
      return false;

    // Check if the nonzeros follow in increasing order
    if (s_.start_ != 0) return false;
    if (s_.step_ != 1) return false;
    if (s_.stop_ != this->nnz()) return false;

    // True if reached this point
    return true;
  }

  template<bool Add>
  void SetNonzerosVector<Add>::generate(std::ostream &stream,
                                                 const std::vector<int>& arg,
                                                 const std::vector<int>& res,
                                                 CodeGenerator& gen) const {
    // Copy first argument if not inplace
    if (arg[0]!=res[0]) {
      gen.copyVector(stream, gen.work(arg[0]), this->nnz(), gen.work(res[0]));
    }

    // Condegen the indices
    int ind = gen.getConstant(this->nz_, true);

    // Perform the operation inplace
    stream << "  for (ii=s" << ind << ", rr=" << gen.work(res[0]) << ", "
           << "ss=" << gen.work(arg[1]) << "; ii!=s" << ind
           << "+" << this->nz_.size() << "; ++ii, ++ss)";
    stream << " if (*ii>=0) rr[*ii] " << (Add?"+=":"=") << " *ss;" << endl;
  }

  template<bool Add>
  void SetNonzerosSlice<Add>::generate(std::ostream &stream,
                                                const std::vector<int>& arg,
                                                const std::vector<int>& res,
                                                CodeGenerator& gen) const {
    // Copy first argument if not inplace
    if (arg[0]!=res[0]) {
      gen.copyVector(stream, gen.work(arg[0]), this->nnz(), gen.work(res[0]));
    }

    // Perform the operation inplace
    stream << "  for (rr=" << gen.work(res[0]+s_.start_) << ", ss="
           << gen.work(arg[1]) << "; rr!=" << gen.work(res[0]+s_.stop_)
           << "; rr+=" << s_.step_ << ")";
    stream << " *rr " << (Add?"+=":"=") << " *ss++;" << endl;
  }

  template<bool Add>
  void SetNonzerosSlice2<Add>::generate(std::ostream &stream,
                                                 const std::vector<int>& arg,
                                                 const std::vector<int>& res,
                                                 CodeGenerator& gen) const {
    // Copy first argument if not inplace
    if (arg[0]!=res[0]) {
      gen.copyVector(stream, gen.work(arg[0]), this->nnz(), gen.work(res[0]));
    }

    // Perform the operation inplace
    stream << "  for (rr=" << gen.work(res[0]+outer_.start_)
           << ", ss=" << gen.work(arg[1]) << "; rr!=" << gen.work(res[0]+outer_.stop_)
           << "; rr+=" << outer_.step_ << ")";
    stream << " for (tt=rr+" << inner_.start_ << "; tt!=rr+" << inner_.stop_
           << "; tt+=" << inner_.step_ << ")";
    stream << " *tt " << (Add?"+=":"=") << " *ss++;" << endl;
  }

  template<bool Add>
  void SetNonzerosSlice<Add>::simplifyMe(MX& ex) {
    // Simplify if addition
    if (isAssignment()) {
      MX t = this->dep(1);
      if (Add) {
        ex += t;
      } else {
        ex = t;
      }
    }
  }

} // namespace casadi

/// \endcond

#endif // CASADI_SETNONZEROS_IMPL_HPP
