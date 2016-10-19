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

using namespace std;

namespace casadi {

  GetNonzeros::GetNonzeros(const Sparsity& sp, const MX& y) {
    setSparsity(sp);
    setDependencies(y);
  }

  void GetNonzerosVector::
  eval(const double** arg, double** res, int* iw, double* w, int mem) const {
    evalGen<double>(arg, res, iw, w);
  }

  void GetNonzerosVector::
  eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    evalGen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  void GetNonzerosVector::
  evalGen(const T* const* arg, T* const* res, int* iw, T* w) const {
    const T* idata = arg[0];
    T* odata = res[0];
    for (auto&& k : nz_) {
      *odata++ = k>=0 ? idata[k] : 0;
    }
  }

  void GetNonzerosSlice::
  eval(const double** arg, double** res, int* iw, double* w, int mem) const {
    evalGen<double>(arg, res, iw, w);
  }

  void GetNonzerosSlice::
  eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    evalGen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  void GetNonzerosSlice::evalGen(const T* const* arg, T* const* res,
                                 int* iw, T* w) const {
    const T* idata = arg[0] + s_.start;
    const T* idata_stop = arg[0] + s_.stop;
    T* odata = res[0];
    for (; idata != idata_stop; idata += s_.step) {
      *odata++ = *idata;
    }
  }

  void GetNonzerosSlice2::
  eval(const double** arg, double** res, int* iw, double* w, int mem) const {
    evalGen<double>(arg, res, iw, w);
  }

  void GetNonzerosSlice2::
  eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    evalGen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  void GetNonzerosSlice2::
  evalGen(const T* const* arg, T* const* res, int* iw, T* w) const {
    const T* outer = arg[0] + outer_.start;
    const T* outer_stop = arg[0] + outer_.stop;
    T* odata = res[0];
    for (; outer != outer_stop; outer += outer_.step) {
      for (const T* inner = outer+inner_.start;
          inner != outer+inner_.stop;
          inner += inner_.step) {
        *odata++ = *inner;
      }
    }
  }

  void GetNonzerosVector::
  sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    const bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (auto&& k : nz_) *r++ = k>=0 ? a[k] : 0;
  }

  void GetNonzerosVector::
  sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (auto&& k : nz_) {
      if (k>=0) a[k] |= *r;
      *r++ = 0;
    }
  }

  void GetNonzerosSlice::
  sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    const bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (int k=s_.start; k!=s_.stop; k+=s_.step) {
      *r++ = a[k];
    }
  }

  void GetNonzerosSlice::
  sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (int k=s_.start; k!=s_.stop; k+=s_.step) {
      a[k] |= *r;
      *r++ = 0;
    }
  }

  void GetNonzerosSlice2::
  sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    const bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (int k1=outer_.start; k1!=outer_.stop; k1+=outer_.step) {
      for (int k2=k1+inner_.start; k2!=k1+inner_.stop; k2+=inner_.step) {
        *r++ = a[k2];
      }
    }
  }

  void GetNonzerosSlice2::
  sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (int k1=outer_.start; k1!=outer_.stop; k1+=outer_.step) {
      for (int k2=k1+inner_.start; k2!=k1+inner_.stop; k2+=inner_.step) {
        a[k2] |= *r;
        *r++ = 0;
      }
    }
  }

  std::string GetNonzerosVector::print(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << nz_;
    return ss.str();
  }

  std::string GetNonzerosSlice::print(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << "[" << s_ << "]";
    return ss.str();
  }

  std::string GetNonzerosSlice2::print(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << "[" << outer_ << ";" << inner_ << "]";
    return ss.str();
  }

  void GetNonzeros::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) {
    // Get all the nonzeros
    vector<int> nz = all();

    // Output sparsity
    const Sparsity& osp = sparsity();
    const int* orow = osp.row();
    vector<int> ocol = osp.get_col();

    // Input sparsity
    const Sparsity& isp = dep().sparsity();
    //const vector<int>& irow = isp.row();
    vector<int> icol = isp.get_col();

    // Get all input elements
    vector<int> el_input;
    isp.find(el_input);

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<int> r_colind, r_row, r_nz, r_ind;

    // Get the matching nonzeros
    r_ind.resize(el_input.size());
    copy(el_input.begin(), el_input.end(), r_ind.begin());
    arg[0].sparsity().get_nz(r_ind);

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
      res[0] = MX(osp.size());
    } else {
      Sparsity f_sp(osp.size1(), osp.size2(), r_colind, r_row);
      res[0] = arg[0]->getGetNonzeros(f_sp, r_nz);
    }
  }

  void GetNonzeros::evalFwd(const std::vector<std::vector<MX> >& fseed,
                            std::vector<std::vector<MX> >& fsens) {
    // Get all the nonzeros
    vector<int> nz = all();

    // Number of derivative directions
    int nfwd = fsens.size();

    // Output sparsity
    const Sparsity& osp = sparsity();
    const int* orow = osp.row();
    vector<int> ocol = osp.get_col();

    // Input sparsity
    const Sparsity& isp = dep().sparsity();
    //const vector<int>& irow = isp.row();
    vector<int> icol = isp.get_col();

    // Get all input elements
    vector<int> el_input;
    isp.find(el_input);

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<int> r_colind, r_row, r_nz, r_ind;

    // Nondifferentiated function and forward sensitivities
    for (int d=0; d<nfwd; ++d) {

      // Get references to arguments and results
      const MX& arg = fseed[d][0];
      MX& res = fsens[d][0];

      // Get the matching nonzeros
      r_ind.resize(el_input.size());
      copy(el_input.begin(), el_input.end(), r_ind.begin());
      arg.sparsity().get_nz(r_ind);

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
        res = MX(osp.size());
      } else {
        Sparsity f_sp(osp.size1(), osp.size2(), r_colind, r_row);
        res = arg->getGetNonzeros(f_sp, r_nz);
      }
    }
  }

  void GetNonzeros::evalAdj(const std::vector<std::vector<MX> >& aseed,
                            std::vector<std::vector<MX> >& asens) {
    // Get all the nonzeros
    vector<int> nz = all();

    // Number of derivative directions
    int nadj = aseed.size();

    // Output sparsity
    const Sparsity& osp = sparsity();
    vector<int> ocol = osp.get_col();

    // Input sparsity
    const Sparsity& isp = dep().sparsity();
    //const vector<int>& irow = isp.row();
    vector<int> icol = isp.get_col();

    // Get all input elements
    vector<int> el_input;
    isp.find(el_input);

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<int> r_colind, r_row, r_nz, r_ind;

    // Adjoint sensitivities
    for (int d=0; d<nadj; ++d) {

      // Get an owning references to the seeds and sensitivities
      // and clear the seeds for the next run
      MX aseed0 = aseed[d][0];
      MX asens0 = asens[d][0]; // Sensitivity before addition

      // Get the corresponding nz locations in the output sparsity pattern
      aseed0.sparsity().find(r_nz);
      osp.get_nz(r_nz);

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
      asens0.sparsity().get_nz(r_ind);

      // Enlarge the sparsity pattern of the sensitivity if not all additions fit
      for (vector<int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
        if (*k>=0 && r_ind[nz[*k]]<0) {

          // Create a new pattern which includes both the the previous seed and the addition
          Sparsity sp = asens0.sparsity().unite(dep().sparsity());
          asens0 = asens0->getProject(sp);

          // Recalculate the nz locations in the adjoint sensitivity corresponding to the inputs
          copy(el_input.begin(), el_input.end(), r_ind.begin());
          asens0.sparsity().get_nz(r_ind);

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
      asens[d][0] = aseed0->getAddNonzeros(asens0, r_nz);
    }
  }

  Matrix<int> GetNonzeros::mapping() const {
    vector<int> nz = all();
    return Matrix<int>(sparsity(), nz, false);
  }

  bool GetNonzerosSlice::is_identity() const {
    // Check sparsity
    if (!(sparsity() == dep().sparsity()))
      return false;

    // Check if the nonzeros follow in increasing order
    if (s_.start != 0) return false;
    if (s_.step != 1) return false;
    if (s_.stop != nnz()) return false;

    // True if reached this point
    return true;
  }

  void GetNonzerosVector::generate(CodeGenerator& g, const std::string& mem,
                                   const std::vector<int>& arg, const std::vector<int>& res) const {
    // Codegen the indices
    int ind = g.getConstant(nz_, true);

    // Codegen the assignments
    g.body << "  for (cii=s" << ind << ", rr=" << g.work(res[0], nnz())
           << ", ss=" << g.work(arg[0], dep(0).nnz())
           << "; cii!=s" << ind << "+" << nz_.size()
           << "; ++cii) *rr++ = *cii>=0 ? ss[*cii] : 0;" << endl;
  }

  void GetNonzerosSlice::simplifyMe(MX& ex) {
    // Simplify if identity
    if (is_identity()) {
      MX t = dep(0);
      ex = t;
    }
  }

  MX GetNonzeros::getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const {
    // Get all the nonzeros
    vector<int> nz_all = all();

    // Eliminate recursive calls
    vector<int> nz_new(nz);
    for (vector<int>::iterator i=nz_new.begin(); i!=nz_new.end(); ++i) {
      if (*i>=0) *i = nz_all[*i];
    }
    return dep()->getGetNonzeros(sp, nz_new);
  }

  void GetNonzerosSlice::generate(CodeGenerator& g, const std::string& mem,
                                  const std::vector<int>& arg, const std::vector<int>& res) const {
    g.body << "  for (rr=" << g.work(res[0], nnz()) << ", ss="
           << g.work(arg[0], dep(0).nnz()) << "+" << s_.start
           << "; ss!=" << g.work(arg[0], dep(0).nnz()) << "+" << s_.stop
           << "; ss+=" << s_.step << ") "
           << "*rr++ = *ss;" << endl;
  }

  void GetNonzerosSlice2::generate(CodeGenerator& g, const std::string& mem,
                                   const std::vector<int>& arg, const std::vector<int>& res) const {
    g.body << "  for (rr=" << g.work(res[0], nnz()) << ", ss="
           << g.work(arg[0], dep(0).nnz()) << "+" << outer_.start
           << "; ss!=" << g.work(arg[0], dep(0).nnz()) << "+" << outer_.stop << "; ss+="
           << outer_.step << ") "
           << "for (tt=ss+" << inner_.start << "; tt!=ss+" << inner_.stop
           << "; tt+=" << inner_.step << ") "
           << "*rr++ = *tt;" << endl;
  }

  bool GetNonzerosVector::is_equal(const MXNode* node, int depth) const {
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

  bool GetNonzerosSlice::is_equal(const MXNode* node, int depth) const {
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

  bool GetNonzerosSlice2::is_equal(const MXNode* node, int depth) const {
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
