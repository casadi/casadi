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


#ifndef CASADI_SETNONZEROS_PARAM_IMPL_HPP
#define CASADI_SETNONZEROS_PARAM_IMPL_HPP

#include "setnonzeros_param.hpp"
#include "casadi_misc.hpp"
#include "serializing_stream.hpp"

/// \cond INTERNAL

using namespace std;

namespace casadi {

  template<bool Add>
  MX SetNonzerosParam<Add>::create(const MX& y, const MX& x, const MX& nz) {
    return MX::create(new SetNonzerosParamVector<Add>(y, x, nz));
  }

  template<bool Add>
  MX SetNonzerosParam<Add>::create(const MX& y, const MX& x, const Slice& s, const MX& nz_offset) {
     casadi_assert(nz_offset.is_scalar(), "nz_offset must be scalar");
    return MX::create(new SetNonzerosParamSlice<Add>(y, x, s, nz_offset));
  }

  template<bool Add>
  MX SetNonzerosParam<Add>::create(const MX& y, const MX& x, const Slice& inner, const Slice& outer, const MX& nz_offset) {
    casadi_assert(nz_offset.numel()==2, "nz_offset must be 2-vector");
    return MX::create(new SetNonzerosParamSlice2<Add>(y, x, inner, outer, nz_offset));
  }

  template<bool Add>
  SetNonzerosParam<Add>::SetNonzerosParam(const MX& y, const MX& x, const MX& nz) {
    this->set_sparsity(y.sparsity());
    this->set_dep(y, x, nz);
  }

  template<bool Add>
  SetNonzerosParamVector<Add>::SetNonzerosParamVector(const MX& y, const MX& x,
      const MX& nz) : SetNonzerosParam<Add>(y, x, nz) {
  }

  template<bool Add>
  SetNonzerosParam<Add>:: ~SetNonzerosParam() {
  }

  template<bool Add>
  void SetNonzerosParam<Add>::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    // Add to the element to the sensitivity, if any
    MX arg0 = project(arg[0], dep(0).sparsity());
    MX arg1 = project(arg[1], dep(1).sparsity());
    MX nz = project(arg[2], dep(2).sparsity());
    if (Add) {
      res[0] = arg1->get_nzadd(arg0, nz);
    } else {
      res[0] = arg1->get_nzassign(arg0, nz);
    }
  }

  template<bool Add>
  void SetNonzerosParam<Add>::ad_forward(const std::vector<std::vector<MX> >& fseed,
                                 std::vector<std::vector<MX> >& fsens) const {
    const MX& nz = dep(2);
    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<fsens.size(); ++d) {

      // Get references to arguments and results
      MX arg0 = project(fseed[d][0], dep(0).sparsity());
      MX arg1 = project(fseed[d][1], dep(1).sparsity());

      /*
        dep(0) <-> y
        dep(1) <-> x

        y[nz]+=x

        dot(y)[nz]+=dot(x)

        dot(x)->get_nzadd(dot(y), nz)

      */
      
      MX& res = fsens[d][0];
      res = arg0;
      
      if (Add) {
        res = arg1->get_nzadd(res, nz);
      } else {
        res = arg1->get_nzassign(res, nz);
      }
    }
  }

  template<bool Add>
  void SetNonzerosParam<Add>::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                                 std::vector<std::vector<MX> >& asens) const {
    const MX& nz = dep(2);
    for (casadi_int d=0; d<aseed.size(); ++d) {
      MX seed = project(aseed[d][0], sparsity());

      /*
        dep(0) <-> y
        dep(1) <-> x

        z: y[nz]+=x

        bar(x) += bar(z)[nz]
        bar(y) += bar(z)
      */
      asens[d][1] += seed->get_nz_ref(nz);
      if (!Add) {
        asens[d][0] += MX::zeros(dep(1).sparsity())->get_nzassign(seed, nz);
      } else {
        asens[d][0] += seed;
      }
    }
  }

  template<bool Add>
  int SetNonzerosParamVector<Add>::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    const double* idata0 = arg[0];
    const double* idata = arg[1];
    const double* nz = arg[2];
    double* odata = res[0];
    // Dimensions
    casadi_int nnz = this->dep(1).nnz();
    casadi_int max_ind = this->dep(0).nnz();
    if (idata0 != odata) {
      copy(idata0, idata0+this->dep(0).nnz(), odata);
    }
    for (casadi_int k=0; k<nnz; ++k) {
      // Get index
      casadi_int index = static_cast<casadi_int>(*nz++);
      if (Add) {
        if (index>=0 && index<max_ind) odata[index] += *idata++;
      } else {
        if (index>=0 && index<max_ind) odata[index] = *idata++;
      }
    }
    return 0;
  }

  template<bool Add>
  int SetNonzerosParamSlice<Add>::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    /**const double* idata0 = arg[0];
    const double* idata = arg[1];
    double* odata = res[0];
    if (idata0 != odata) {
      copy(idata0, idata0+this->dep(0).nnz(), odata);
    }
    double* odata_stop = odata + s_.stop;
    for (odata += s_.start; odata != odata_stop; odata += s_.step) {
      if (Add) {
        *odata += *idata++;
      } else {
        *odata = *idata++;
      }
    }*/
    return 0;
  }

  template<bool Add>
  int SetNonzerosParamSlice2<Add>::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    /**const double* idata0 = arg[0];
    const double* idata = arg[1];
    double* odata = res[0];
    if (idata0 != odata) {
      copy(idata0, idata0 + this->dep(0).nnz(), odata);
    }
    double* outer_stop = odata + outer_.stop;
    double* outer = odata + outer_.start;
    for (; outer != outer_stop; outer += outer_.step) {
      for (double* inner = outer+inner_.start;
          inner != outer+inner_.stop;
          inner += inner_.step) {
        if (Add) {
          *inner += *idata++;
        } else {
          *inner = *idata++;
        }
      }
    }*/
    return 0;
  }

  template<bool Add>
  int SetNonzerosParamVector<Add>::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // Parametric index -> any input disturbance propagates to any output
    bvec_t arg0 = bvec_or(arg[0], this->dep(0).nnz());
    bvec_t arg1 = bvec_or(arg[1], this->dep(1).nnz());

    bvec_t *r = res[0];
    std::fill(r, r+this->nnz(), arg0 | arg1);
    return 0;
  }

  template<bool Add>
  int SetNonzerosParamVector<Add>::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *arg0 = arg[0];
    bvec_t *arg1 = arg[1];
    bvec_t r = bvec_or(res[0], this->nnz());
    std::fill(res[0], res[0]+this->nnz(), bvec_t(0));

    for (casadi_int i=0;i<this->dep(0).nnz();++i) {
      *arg0++ |= r;
    }
    for (casadi_int i=0;i<this->dep(1).nnz();++i) {
      *arg1++ |= r;
    }
    return 0;
  }

  template<bool Add>
  int SetNonzerosParamSlice<Add>::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    return 0;
  }

  template<bool Add>
  int SetNonzerosParamSlice<Add>::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    return 0;
  }

  template<bool Add>
  int SetNonzerosParamSlice2<Add>::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    return 0;
  }

  template<bool Add>
  int SetNonzerosParamSlice2<Add>::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    return 0;
  }

  template<bool Add>
  std::string SetNonzerosParamVector<Add>::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << "(" << arg.at(0) << "[" << arg.at(2) << "]" << (Add ? " += " : " = ") << arg.at(1) << ")";
    return ss.str();
  }

  template<bool Add>
  std::string SetNonzerosParamSlice<Add>::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << "(" << arg.at(0) << "[" << arg.at(2) << "+" << s_ << "]" << (Add ? " += " : " = ") << arg.at(1) << ")";
    return ss.str();
  }

  template<bool Add>
  std::string SetNonzerosParamSlice2<Add>::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << "(" << arg.at(0) << "[" << arg.at(2) << "+("<< outer_ << ";" << inner_ << ")]" << (Add ? " += " : " = ")
       << arg.at(1) << ")";
    return ss.str();
  }

  template<bool Add>
  void SetNonzerosParamVector<Add>::
  generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {
    // Copy first argument if not inplace
    if (arg[0]!=res[0]) {
      g << g.copy(g.work(arg[0], this->dep(0).nnz()), this->nnz(),
                          g.work(res[0], this->nnz())) << '\n';
    }

    casadi_int n = this->dep(1).nnz();

    g.local("i", "casadi_int");
    g.local("cr", "const casadi_real", "*");
    g.local("cs", "const casadi_real", "*");
    g << "for (cs=" << g.work(arg[1], n) << ", cr=" << g.work(arg[2], n)
      << "; cs!=" << g.work(arg[1], n) << "+" << n
      << "; ++cs) { i=(int) *cr++; if (i>=0 && i<" << this->dep(0).nnz() << ") "
      << g.work(res[0], this->nnz()) << "[i] " << (Add?"+= ":"= ")
      << "*cs; }\n";
  }

  template<bool Add>
  void SetNonzerosParamSlice<Add>::
  generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {

  }

  template<bool Add>
  void SetNonzerosParamSlice2<Add>::
  generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {

  }

  template<bool Add>
  void SetNonzerosParamVector<Add>::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
  }

  template<bool Add>
  SetNonzerosParamVector<Add>::SetNonzerosParamVector(DeserializingStream& s) : SetNonzerosParam<Add>(s) {
  }

  template<bool Add>
  void SetNonzerosParamVector<Add>::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("SetNonzerosParam::type", 'a');
  }

  template<bool Add>
  void SetNonzerosParamSlice<Add>::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("SetNonzerosParamSlice::slice", s_);
  }

  template<bool Add>
  SetNonzerosParamSlice<Add>::SetNonzerosParamSlice(DeserializingStream& s) : SetNonzerosParam<Add>(s) {
    s.unpack("SetNonzerosParamSlice::slice", s_);
  }

  template<bool Add>
  void SetNonzerosParamSlice<Add>::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("SetNonzerosParam::type", 'b');
  }

  template<bool Add>
  void SetNonzerosParamSlice2<Add>::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("SetNonzerosParamSlice2::inner", inner_);
    s.pack("SetNonzerosParamSlice2::outer", outer_);
  }

  template<bool Add>
  SetNonzerosParamSlice2<Add>::SetNonzerosParamSlice2(DeserializingStream& s) : SetNonzerosParam<Add>(s) {
    s.unpack("SetNonzerosParamSlice2::inner", inner_);
    s.unpack("SetNonzerosParamSlice2::outer", outer_);
  }

  template<bool Add>
  void SetNonzerosParamSlice2<Add>::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("SetNonzerosParam::type", 'c');
  }

  template<bool Add>
  MXNode* SetNonzerosParam<Add>::deserialize(DeserializingStream& s) {
    char t;
    s.unpack("SetNonzerosParam::type", t);
    switch (t) {
      case 'a': return new SetNonzerosParamVector<Add>(s);
      case 'b': return new SetNonzerosParamSlice<Add>(s);
      case 'c': return new SetNonzerosParamSlice2<Add>(s);
      default: casadi_assert_dev(false);
    }
  }

} // namespace casadi

/// \endcond

#endif // CASADI_SETNONZEROS_PARAM_IMPL_HPP
