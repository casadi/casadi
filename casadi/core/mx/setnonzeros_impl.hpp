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
#include "../matrix/sparsity_tools.hpp"

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
      vector<bool> already_set(this->size(), false);
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
  void SetNonzeros<Add>::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                                    MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                                    bool output_given) {
    // Get all the nonzeros
    vector<int> nz = getAll();

    // Number of derivative directions
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Output sparsity
    const Sparsity &osp = sparsity();
    const vector<int>& orow = osp.row();
    vector<int> ocol = osp.getCol();

    // Input sparsity (first input same as output)
    const Sparsity &isp = dep(1).sparsity();
    const vector<int>& irow = isp.row();
    vector<int> icol = isp.getCol();

    // We next need to resort the assignment vector by outputs instead of inputs
    // Start by counting the number of output nonzeros corresponding to each input nonzero
    vector<int> onz_count(orow.size()+2, 0);
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
        with_duplicates[k] = ocol[onz_k] + orow[onz_k]*osp.size2();
      } else {
        with_duplicates[k] = -1;
      }
    }

    // Get all output elements (this time without duplicates)
    vector<int> el_output;
    osp.getElements(el_output, false);

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<int> r_colind, r_row, r_nz, r_ind;

    // Nondifferentiated function and forward sensitivities
    int first_d = output_given ? 0 : -1;
    for (int d=first_d; d<nfwd; ++d) {

      // Get references to arguments and results
      const MX& arg = d<0 ? *input[1] : *fwdSeed[d][1];
      const MX& arg0 = d<0 ? *input[0] : *fwdSeed[d][0];
      MX& res = d<0 ? *output[0] : *fwdSens[d][0];
      if (&arg0 != &res) {
        res = arg0;
      }

      // Entries in res with elements zero'ed out
      if (!Add) {

        // Get the nz locations in res corresponding to the output sparsity pattern
        r_nz.resize(with_duplicates.size());
        copy(with_duplicates.begin(), with_duplicates.end(), r_nz.begin());
        res.sparsity().getNZInplace(r_nz);

        // Zero out the corresponding entries
        res = MX::zeros(isp)->getSetNonzeros(res, r_nz);
      }

      // Get the nz locations of the elements in arg corresponding to the argument sparsity pattern
      arg.sparsity().getElements(r_nz, false);
      isp.getNZInplace(r_nz);

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
      res.sparsity().getNZInplace(r_ind);

      // Enlarge the sparsity pattern of the arguments if not all assignments fit
      for (vector<int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
        if (*k>=0 && nz[*k]>=0 && r_ind[nz[*k]]<0) {

          // Create a new pattern which includes both the the previous seed
          // and the addition/assignment
          Sparsity sp = res.sparsity().patternUnion(osp);
          res = res->getSetSparse(sp);

          // Recalculate the nz locations in the arguments corresponding to the inputs
          copy(el_output.begin(), el_output.end(), r_ind.begin());
          res.sparsity().getNZInplace(r_ind);

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

    // Adjoint sensitivities
    for (int d=0; d<nadj; ++d) {

      // Get an owning references to the seeds and sensitivities
      // and clear the seeds for the next run
      MX& aseed = *adjSeed[d][0];
      MX& asens0 = *adjSens[d][0];
      MX& asens = *adjSens[d][1];

      // Get the matching nonzeros
      r_ind.resize(el_output.size());
      copy(el_output.begin(), el_output.end(), r_ind.begin());
      aseed.sparsity().getNZInplace(r_ind);

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
  void SetNonzerosVector<Add>::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output,
                                         std::vector<int>& itmp, std::vector<double>& rtmp) {
    evaluateGen<double, DMatrixPtrV, DMatrixPtrVV>(input, output, itmp, rtmp);
  }

  template<bool Add>
  void SetNonzerosVector<Add>::evaluateSX(const SXPtrV& input, SXPtrV& output,
                                          std::vector<int>& itmp, std::vector<SXElement>& rtmp) {
    evaluateGen<SXElement, SXPtrV, SXPtrVV>(input, output, itmp, rtmp);
  }

  template<bool Add>
  template<typename T, typename MatV, typename MatVV>
  void SetNonzerosVector<Add>::evaluateGen(const MatV& input, MatV& output,
                                           std::vector<int>& itmp, std::vector<T>& rtmp) {

    const vector<T>& idata0 = input[0]->data();
    typename vector<T>::const_iterator idata_it = input[1]->begin();
    vector<T>& odata = output[0]->data();
    if (&idata0 != &odata) {
      copy(idata0.begin(), idata0.end(), odata.begin());
    }
    for (vector<int>::const_iterator k=this->nz_.begin(); k!=this->nz_.end(); ++k, ++idata_it) {
      if (Add) {
        if (*k>=0) odata[*k] += *idata_it;
      } else {
        if (*k>=0) odata[*k] = *idata_it;
      }
    }
  }

  template<bool Add>
  void SetNonzerosSlice<Add>::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output,
                                        std::vector<int>& itmp, std::vector<double>& rtmp) {
    evaluateGen<double, DMatrixPtrV, DMatrixPtrVV>(input, output, itmp, rtmp);
  }

  template<bool Add>
  void SetNonzerosSlice<Add>::evaluateSX(const SXPtrV& input, SXPtrV& output,
                                         std::vector<int>& itmp, std::vector<SXElement>& rtmp) {
    evaluateGen<SXElement, SXPtrV, SXPtrVV>(input, output, itmp, rtmp);
  }

  template<bool Add>
  template<typename T, typename MatV, typename MatVV>
  void SetNonzerosSlice<Add>::evaluateGen(const MatV& input, MatV& output,
                                          std::vector<int>& itmp, std::vector<T>& rtmp) {

    const vector<T>& idata0 = input[0]->data();
    vector<T>& odata = output[0]->data();
    if (&idata0 != &odata) {
      copy(idata0.begin(), idata0.end(), odata.begin());
    }
    const vector<T>& idata = input[1]->data();
    const T* idata_ptr = getPtr(idata);
    T* odata_ptr = getPtr(odata) + s_.start_;
    T* odata_stop = getPtr(odata) + s_.stop_;
    for (; odata_ptr != odata_stop; odata_ptr += s_.step_) {
      if (Add) {
        *odata_ptr += *idata_ptr++;
      } else {
        *odata_ptr = *idata_ptr++;
      }
    }
  }

  template<bool Add>
  void SetNonzerosSlice2<Add>::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output,
                                         std::vector<int>& itmp, std::vector<double>& rtmp) {
    evaluateGen<double, DMatrixPtrV, DMatrixPtrVV>(input, output, itmp, rtmp);
  }

  template<bool Add>
  void SetNonzerosSlice2<Add>::evaluateSX(const SXPtrV& input, SXPtrV& output,
                                          std::vector<int>& itmp, std::vector<SXElement>& rtmp) {
    evaluateGen<SXElement, SXPtrV, SXPtrVV>(input, output, itmp, rtmp);
  }

  template<bool Add>
  template<typename T, typename MatV, typename MatVV>
  void SetNonzerosSlice2<Add>::evaluateGen(const MatV& input, MatV& output,
                                           std::vector<int>& itmp, std::vector<T>& rtmp) {

    const vector<T>& idata0 = input[0]->data();
    vector<T>& odata = output[0]->data();
    if (&idata0 != &odata) {
      copy(idata0.begin(), idata0.end(), odata.begin());
    }
    const vector<T>& idata = input[1]->data();
    const T* idata_ptr = getPtr(idata);
    T* outer_ptr = getPtr(odata) + outer_.start_;
    T* outer_stop = getPtr(odata) + outer_.stop_;
    for (; outer_ptr != outer_stop; outer_ptr += outer_.step_) {
      for (T* inner_ptr = outer_ptr+inner_.start_;
          inner_ptr != outer_ptr+inner_.stop_;
          inner_ptr += inner_.step_) {
        if (Add) {
          *inner_ptr += *idata_ptr++;
        } else {
          *inner_ptr = *idata_ptr++;
        }
      }
    }
  }

  template<bool Add>
  void SetNonzerosVector<Add>::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output,
                                                 bool fwd) {
    // Get references to the assignment operations and data
    bvec_t *outputd = get_bvec_t(output[0]->data());
    bvec_t *inputd0 = get_bvec_t(input[0]->data());
    bvec_t *inputd = get_bvec_t(input[1]->data());

    // Propagate sparsity
    if (fwd) {
      if (outputd != inputd0) {
        copy(inputd0, inputd0+input[0]->size(), outputd);
      }
      for (vector<int>::const_iterator k=this->nz_.begin(); k!=this->nz_.end(); ++k, ++inputd) {
        if (Add) {
          if (*k>=0) outputd[*k] |= *inputd;
        } else {
          if (*k>=0) outputd[*k] = *inputd;
        }
      }
    } else {
      for (vector<int>::const_iterator k=this->nz_.begin(); k!=this->nz_.end(); ++k, ++inputd) {
        if (*k>=0) {
          *inputd |= outputd[*k];
          if (!Add) {
            outputd[*k] = 0;
          }
        }
      }
      if (outputd != inputd0) {
        int n = input[0]->size();
        for (int k=0; k<n; ++k) {
          inputd0[k] |= outputd[k];
          outputd[k] = 0;
        }
      }
    }
  }

  template<bool Add>
  void SetNonzerosSlice<Add>::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd) {
    // Get references to the assignment operations and data
    bvec_t *outputd = get_bvec_t(output[0]->data());
    bvec_t *inputd0 = get_bvec_t(input[0]->data());
    bvec_t *inputd = get_bvec_t(input[1]->data());

    // Propagate sparsity
    if (fwd) {
      if (outputd != inputd0) {
        copy(inputd0, inputd0+input[0]->size(), outputd);
      }
      for (int k=s_.start_; k!=s_.stop_; k+=s_.step_) {
        if (Add) {
          outputd[k] |= *inputd++;
        } else {
          outputd[k] = *inputd++;
        }
      }
    } else {
      for (int k=s_.start_; k!=s_.stop_; k+=s_.step_) {
        *inputd++ |= outputd[k];
        if (!Add) {
          outputd[k] = 0;
        }
      }
      if (outputd != inputd0) {
        int n = input[0]->size();
        for (int k=0; k<n; ++k) {
          inputd0[k] |= outputd[k];
          outputd[k] = 0;
        }
      }
    }
  }

  template<bool Add>
  void SetNonzerosSlice2<Add>::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output,
                                                 bool fwd) {
    // Get references to the assignment operations and data
    bvec_t *outputd = get_bvec_t(output[0]->data());
    bvec_t *inputd0 = get_bvec_t(input[0]->data());
    bvec_t *inputd = get_bvec_t(input[1]->data());

    // Propagate sparsity
    if (fwd) {
      if (outputd != inputd0) {
        copy(inputd0, inputd0+input[0]->size(), outputd);
      }
      for (int k1=outer_.start_; k1!=outer_.stop_; k1+=outer_.step_) {
        for (int k2=k1+inner_.start_; k2!=k1+inner_.stop_; k2+=inner_.step_) {
          if (Add) {
            outputd[k2] |= *inputd++;
          } else {
            outputd[k2] = *inputd++;
          }
        }
      }
    } else {
      for (int k1=outer_.start_; k1!=outer_.stop_; k1+=outer_.step_) {
        for (int k2=k1+inner_.start_; k2!=k1+inner_.stop_; k2+=inner_.step_) {
          *inputd++ |= outputd[k2];
          if (!Add) {
            outputd[k2] = 0;
          }
        }
      }
      if (outputd != inputd0) {
        int n = input[0]->size();
        for (int k=0; k<n; ++k) {
          inputd0[k] |= outputd[k];
          outputd[k] = 0;
        }
      }
    }
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
    return Matrix<int>(this->dep(1).sparsity(), nz);
  }

  template<bool Add>
  bool SetNonzerosVector<Add>::isEqual(const MXNode* node, int depth) const {
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
  bool SetNonzerosSlice<Add>::isEqual(const MXNode* node, int depth) const {
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
  bool SetNonzerosSlice2<Add>::isEqual(const MXNode* node, int depth) const {
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
    if (s_.stop_ != this->size()) return false;

    // True if reached this point
    return true;
  }

  template<bool Add>
  void SetNonzerosVector<Add>::generateOperation(std::ostream &stream,
                                                 const std::vector<std::string>& arg,
                                                 const std::vector<std::string>& res,
                                                 CodeGenerator& gen) const {
    // Check if inplace
    bool inplace = arg.at(0).compare(res.front())==0;

    // Copy first argument if not inplace
    if (!inplace) {
      stream << "  for (i=0; i<" << this->size() << "; ++i) " << res.front() << "[i]="
             << arg.at(0) << "[i];" << endl;
    }

    // Condegen the indices
    int ind = gen.getConstant(this->nz_, true);

    // Perform the operation inplace
    stream << "  for (ii=s" << ind << ", rr=" << res.front() << ", ss=" << arg.at(1) << "; ii!=s"
           << ind << "+" << this->nz_.size() << "; ++ii, ++ss)";
    stream << " if (*ii>=0) rr[*ii] " << (Add?"+=":"=") << " *ss;" << endl;
  }

  template<bool Add>
  void SetNonzerosSlice<Add>::generateOperation(std::ostream &stream,
                                                const std::vector<std::string>& arg,
                                                const std::vector<std::string>& res,
                                                CodeGenerator& gen) const {
    // Check if inplace
    bool inplace = arg.at(0).compare(res.front())==0;

    // Copy first argument if not inplace
    if (!inplace) {
      stream << "  for (i=0; i<" << this->size() << "; ++i) " << res.front()
             << "[i]=" << arg.at(0) << "[i];" << endl;
    }

    // Perform the operation inplace
    stream << "  for (rr=" << res.front() << "+" << s_.start_ << ", ss="
           << arg.at(1) << "; rr!=" << res.front() << "+" << s_.stop_
           << "; rr+=" << s_.step_ << ")";
    stream << " *rr " << (Add?"+=":"=") << " *ss++;" << endl;
  }

  template<bool Add>
  void SetNonzerosSlice2<Add>::generateOperation(std::ostream &stream,
                                                 const std::vector<std::string>& arg,
                                                 const std::vector<std::string>& res,
                                                 CodeGenerator& gen) const {
    // Check if inplace
    bool inplace = arg.at(0).compare(res.front())==0;

    // Copy first argument if not inplace
    if (!inplace) {
      stream << "  for (i=0; i<" << this->size() << "; ++i) "
             << res.front() << "[i]=" << arg.at(0) << "[i];" << endl;
    }

    // Perform the operation inplace
    stream << "  for (rr=" << res.front() << "+" << outer_.start_
           << ", ss=" << arg.at(1) << "; rr!=" << res.front()
           << "+" << outer_.stop_ << "; rr+=" << outer_.step_ << ")";
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

#endif // CASADI_SETNONZEROS_IMPL_HPP
