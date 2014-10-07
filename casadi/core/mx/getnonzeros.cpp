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


#include "getnonzeros.hpp"
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../function/sx_function.hpp"
#include "../matrix/sparsity_tools.hpp"

using namespace std;

namespace casadi {

  GetNonzeros::GetNonzeros(const Sparsity& sp, const MX& y) {
    setSparsity(sp);
    setDependencies(y);
  }

  void GetNonzerosVector::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output,
                                    std::vector<int>& itmp, std::vector<double>& rtmp) {
    evaluateGen<double, DMatrixPtrV, DMatrixPtrVV>(input, output, itmp, rtmp);
  }

  void GetNonzerosVector::evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                                     std::vector<SXElement>& rtmp) {
    evaluateGen<SXElement, SXPtrV, SXPtrVV>(input, output, itmp, rtmp);
  }

  template<typename T, typename MatV, typename MatVV>
  void GetNonzerosVector::evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp,
                                      std::vector<T>& rtmp) {
    const vector<T>& idata = input[0]->data();
    typename vector<T>::iterator odata_it = output[0]->begin();
    for (vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k) {
      *odata_it++ = *k>=0 ? idata[*k] : 0;
    }
  }

  void GetNonzerosSlice::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output,
                                   std::vector<int>& itmp, std::vector<double>& rtmp) {
    evaluateGen<double, DMatrixPtrV, DMatrixPtrVV>(input, output, itmp, rtmp);
  }

  void GetNonzerosSlice::evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                                    std::vector<SXElement>& rtmp) {
    evaluateGen<SXElement, SXPtrV, SXPtrVV>(input, output, itmp, rtmp);
  }

  template<typename T, typename MatV, typename MatVV>
  void GetNonzerosSlice::evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp,
                                     std::vector<T>& rtmp) {

    const vector<T>& idata = input[0]->data();
    const T* idata_ptr = getPtr(idata) + s_.start_;
    const T* idata_stop = getPtr(idata) + s_.stop_;
    T* odata_ptr = getPtr(output[0]->data());
    for (; idata_ptr != idata_stop; idata_ptr += s_.step_) {
      *odata_ptr++ = *idata_ptr;
    }
  }

  void GetNonzerosSlice2::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output,
                                    std::vector<int>& itmp, std::vector<double>& rtmp) {
    evaluateGen<double, DMatrixPtrV, DMatrixPtrVV>(input, output, itmp, rtmp);
  }

  void GetNonzerosSlice2::evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                                     std::vector<SXElement>& rtmp) {
    evaluateGen<SXElement, SXPtrV, SXPtrVV>(input, output, itmp, rtmp);
  }

  template<typename T, typename MatV, typename MatVV>
  void GetNonzerosSlice2::evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp,
                                      std::vector<T>& rtmp) {

    const vector<T>& idata = input[0]->data();
    const T* outer_ptr = getPtr(idata) + outer_.start_;
    const T* outer_stop = getPtr(idata) + outer_.stop_;
    T* odata_ptr = getPtr(output[0]->data());
    for (; outer_ptr != outer_stop; outer_ptr += outer_.step_) {
      for (const T* inner_ptr = outer_ptr+inner_.start_;
          inner_ptr != outer_ptr+inner_.stop_;
          inner_ptr += inner_.step_) {
        *odata_ptr++ = *inner_ptr;
      }
    }
  }

  void GetNonzerosVector::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd) {
    // Get references to the assignment operations and data
    bvec_t *outputd = get_bvec_t(output[0]->data());
    bvec_t *inputd = get_bvec_t(input[0]->data());

    // Propagate sparsity
    if (fwd) {
      for (vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k) {
        *outputd++ = *k>=0 ? inputd[*k] : 0;
      }
    } else {
      for (vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k) {
        if (*k>=0) inputd[*k] |= *outputd;
        *outputd++ = 0;
      }
    }
  }

  void GetNonzerosSlice::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd) {
    // Get references to the assignment operations and data
    bvec_t *outputd = get_bvec_t(output[0]->data());
    bvec_t *inputd = get_bvec_t(input[0]->data());

    // Propagate sparsity
    if (fwd) {
      for (int k=s_.start_; k!=s_.stop_; k+=s_.step_) {
        *outputd++ = inputd[k];
      }
    } else {
      for (int k=s_.start_; k!=s_.stop_; k+=s_.step_) {
        inputd[k] |= *outputd;
        *outputd++ = 0;
      }
    }
  }

  void GetNonzerosSlice2::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd) {
    // Get references to the assignment operations and data
    bvec_t *outputd = get_bvec_t(output[0]->data());
    bvec_t *inputd = get_bvec_t(input[0]->data());

    // Propagate sparsity
    if (fwd) {
      for (int k1=outer_.start_; k1!=outer_.stop_; k1+=outer_.step_) {
        for (int k2=k1+inner_.start_; k2!=k1+inner_.stop_; k2+=inner_.step_) {
          *outputd++ = inputd[k2];
        }
      }
    } else {
      for (int k1=outer_.start_; k1!=outer_.stop_; k1+=outer_.step_) {
        for (int k2=k1+inner_.start_; k2!=k1+inner_.stop_; k2+=inner_.step_) {
          inputd[k2] |= *outputd;
          *outputd++ = 0;
        }
      }
    }
  }

  void GetNonzerosVector::printPart(std::ostream &stream, int part) const {
    switch (part) {
    case 1: stream << nz_; break;
    }
  }

  void GetNonzerosSlice::printPart(std::ostream &stream, int part) const {
    switch (part) {
    case 1: stream << "[" << s_ << "]"; break;
    }
  }

  void GetNonzerosSlice2::printPart(std::ostream &stream, int part) const {
    switch (part) {
    case 1: stream << "[" << outer_ << ";" << inner_ << "]"; break;
    }
  }

  void GetNonzeros::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                               MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                               bool output_given) {
    // Get all the nonzeros
    vector<int> nz = getAll();

    // Number of derivative directions
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Output sparsity
    const Sparsity& osp = sparsity();
    const vector<int>& orow = osp.row();
    vector<int> ocol = osp.getCol();

    // Input sparsity
    const Sparsity& isp = dep().sparsity();
    //const vector<int>& irow = isp.row();
    vector<int> icol = isp.getCol();

    // Get all input elements
    vector<int> el_input;
    isp.getElements(el_input, false);

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<int> r_colind, r_row, r_nz, r_ind;

    // Nondifferentiated function and forward sensitivities
    int first_d = output_given ? 0 : -1;
    for (int d=first_d; d<nfwd; ++d) {

      // Get references to arguments and results
      const MX& arg = d<0 ? *input[0] : *fwdSeed[d][0];
      MX& res = d<0 ? *output[0] : *fwdSens[d][0];

      // Get the matching nonzeros
      r_ind.resize(el_input.size());
      copy(el_input.begin(), el_input.end(), r_ind.begin());
      arg.sparsity().getNZInplace(r_ind);

      // Sparsity pattern for the result
      r_colind.resize(osp.size2()+1); // Col count
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
        int i=ocol[k], j=orow[k];

        // Add to sparsity pattern
        r_row.push_back(j);
        r_colind[1+i]++;
      }

      // col count -> col offset
      for (int i=1; i<r_colind.size(); ++i) r_colind[i] += r_colind[i-1];

      // Create a sparsity pattern from vectors
      if (r_nz.size()==0) {
        res = MX::sparse(osp.shape());
      } else {
        Sparsity f_sp(osp.size1(), osp.size2(), r_colind, r_row);
        res = arg->getGetNonzeros(f_sp, r_nz);
      }
    }

    // Adjoint sensitivities
    for (int d=0; d<nadj; ++d) {

      // Get an owning references to the seeds and sensitivities
      // and clear the seeds for the next run
      MX aseed = *adjSeed[d][0];
      *adjSeed[d][0] = MX();
      MX& asens = *adjSens[d][0]; // Sensitivity after addition
      MX asens0 = asens; // Sensitivity before addition

      // Get the corresponding nz locations in the output sparsity pattern
      aseed.sparsity().getElements(r_nz, false);
      osp.getNZInplace(r_nz);

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

      // Quick continue of no elements to add
      if (!elements_to_add) continue;

      // Get the nz locations in the adjoint sensitivity corresponding to the inputs
      r_ind.resize(el_input.size());
      copy(el_input.begin(), el_input.end(), r_ind.begin());
      asens0.sparsity().getNZInplace(r_ind);

      // Enlarge the sparsity pattern of the sensitivity if not all additions fit
      for (vector<int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
        if (*k>=0 && r_ind[nz[*k]]<0) {

          // Create a new pattern which includes both the the previous seed and the addition
          Sparsity sp = asens0.sparsity().patternUnion(dep().sparsity());
          asens0 = asens0->getSetSparse(sp);

          // Recalculate the nz locations in the adjoint sensitivity corresponding to the inputs
          copy(el_input.begin(), el_input.end(), r_ind.begin());
          asens0.sparsity().getNZInplace(r_ind);

          break;
        }
      }

      // Have r_nz point to locations in the sensitivity instead of the output
      for (vector<int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
        if (*k>=0) {
          *k = r_ind[nz[*k]];
        }
      }

      // Add to the element to the sensitivity
      asens = aseed->getAddNonzeros(asens0, r_nz);
    }
  }

  Matrix<int> GetNonzeros::mapping() const {
    vector<int> nz = getAll();
    return Matrix<int>(sparsity(), nz);
  }

  bool GetNonzerosSlice::isIdentity() const {
    // Check sparsity
    if (!(sparsity() == dep().sparsity()))
      return false;

    // Check if the nonzeros follow in increasing order
    if (s_.start_ != 0) return false;
    if (s_.step_ != 1) return false;
    if (s_.stop_ != size()) return false;

    // True if reached this point
    return true;
  }

  void GetNonzerosVector::generateOperation(std::ostream &stream,
                                            const std::vector<std::string>& arg,
                                            const std::vector<std::string>& res,
                                            CodeGenerator& gen) const {
    // Condegen the indices
    int ind = gen.getConstant(nz_, true);

    // Codegen the assignments
    stream << "  for (ii=s" << ind << ", rr=" << res.front() << ", ss=" << arg.front()
           << "; ii!=s" << ind << "+" << nz_.size()
           << "; ++ii) *rr++ = *ii>=0 ? ss[*ii] : 0;" << endl;
  }

  void GetNonzerosSlice::simplifyMe(MX& ex) {
    // Simplify if identity
    if (isIdentity()) {
      MX t = dep(0);
      ex = t;
    }
  }

  MX GetNonzeros::getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const {
    // Get all the nonzeros
    vector<int> nz_all = getAll();

    // Eliminate recursive calls
    vector<int> nz_new(nz);
    for (vector<int>::iterator i=nz_new.begin(); i!=nz_new.end(); ++i) {
      if (*i>=0) *i = nz_all[*i];
    }
    return dep()->getGetNonzeros(sp, nz_new);
  }

  void GetNonzerosSlice::generateOperation(std::ostream &stream,
                                           const std::vector<std::string>& arg,
                                           const std::vector<std::string>& res,
                                           CodeGenerator& gen) const {
    stream << "  for (rr=" << res.front() << ", ss=" << arg.front() << "+" << s_.start_
           << "; ss!=" << arg.front() << "+" << s_.stop_ << "; ss+=" << s_.step_ << ") ";
    stream << "*rr++ = *ss;" << endl;
  }

  void GetNonzerosSlice2::generateOperation(std::ostream &stream,
                                            const std::vector<std::string>& arg,
                                            const std::vector<std::string>& res,
                                            CodeGenerator& gen) const {
    stream << "  for (rr=" << res.front() << ", ss=" << arg.front() << "+" << outer_.start_
           << "; ss!=" << arg.front() << "+" << outer_.stop_ << "; ss+="
           << outer_.step_ << ") ";
    stream << "for (tt=ss+" << inner_.start_ << "; tt!=ss+" << inner_.stop_
           << "; tt+=" << inner_.step_ << ") ";
    stream << "*rr++ = *tt;" << endl;
  }

  bool GetNonzerosVector::isEqual(const MXNode* node, int depth) const {
    // Check dependencies
    if (!sameOpAndDeps(node, depth)) return false;

    // Check if same node
    const GetNonzerosVector* n = dynamic_cast<const GetNonzerosVector*>(node);
    if (n==0) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (this->nz_.size()!=n->nz_.size()) return false;
    if (!std::equal(this->nz_.begin(), this->nz_.end(), n->nz_.begin())) return false;

    return true;
  }

  bool GetNonzerosSlice::isEqual(const MXNode* node, int depth) const {
    // Check dependencies
    if (!sameOpAndDeps(node, depth)) return false;

    // Check if same node
    const GetNonzerosSlice* n = dynamic_cast<const GetNonzerosSlice*>(node);
    if (n==0) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (this->s_ != n->s_) return false;

    return true;
  }

  bool GetNonzerosSlice2::isEqual(const MXNode* node, int depth) const {
    // Check dependencies
    if (!sameOpAndDeps(node, depth)) return false;

    // Check if same node
    const GetNonzerosSlice2* n = dynamic_cast<const GetNonzerosSlice2*>(node);
    if (n==0) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (this->inner_ != n->inner_ || this->outer_!=n->outer_) return false;

    return true;
  }


} // namespace casadi
