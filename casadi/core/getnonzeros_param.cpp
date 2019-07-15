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


#include "getnonzeros_param.hpp"
#include "casadi_misc.hpp"
#include "serializing_stream.hpp"

using namespace std;

namespace casadi {

  MX GetNonzerosParam::create(const MX& x, const MX& nz) {
    // No elements at all
    if (nz.nnz()==0) return MX::zeros(nz.sparsity());
    return MX::create(new GetNonzerosParamVector(x, nz));
  }

  MX GetNonzerosParam::create(const Sparsity& sp, const MX& x, const Slice& s, const MX& nz_offset) {
    casadi_assert(nz_offset.is_scalar(), "nz_offset must be scalar");
    return MX::create(new GetNonzerosParamSlice(sp, x, s, nz_offset));
  }

  MX GetNonzerosParam::create(const Sparsity& sp, const MX& x,
                         const Slice& inner, const Slice& outer,
                         const MX& nz_offset) {
    casadi_assert(nz_offset.numel()==2, "nz_offset must be 2-vector");
    return MX::create(new GetNonzerosParamSlice2(sp, x, inner, outer, nz_offset));
  }

  GetNonzerosParam::GetNonzerosParam(const Sparsity& sp, const MX& y, const MX& nz) {
    set_sparsity(sp);
    set_dep(y, nz);
  }

  int GetNonzerosParamVector::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    const double* idata = arg[0];
    const double* nz = arg[1];
    double* odata = res[0];
    // Dimensions
    casadi_int nnz = dep(1).nnz();
    casadi_int max_ind = dep(0).nnz();
    // Get elements
    for (casadi_int i=0; i<nnz; ++i) {
      // Get index
      casadi_int index = static_cast<casadi_int>(*nz++);
      // Make assignment if in bounds, else NaN
      *odata++ = index>=0 && index<max_ind ? idata[index] : nan;
    }
    return 0;
  }

  int GetNonzerosParamSlice::eval(const double** arg, double** res,
                                 casadi_int* iw, double* w) const {
    casadi_int offset = static_cast<casadi_int>(arg[1][0]);
    const double* idata = arg[0] + offset + s_.start;
    const double* idata_stop = arg[0] + offset + s_.stop;
    double* odata = res[0];
    for (; idata != idata_stop; idata += s_.step) {
      *odata++ = *idata;
    }
    return 0;
  }

  int GetNonzerosParamSlice2::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    casadi_int offset_inner = static_cast<casadi_int>(arg[1][0]);
    casadi_int offset_outer = static_cast<casadi_int>(arg[1][1]);

    const double* outer = arg[0] + offset_outer + outer_.start;
    const double* outer_stop = arg[0] + offset_outer + outer_.stop;
    double* odata = res[0];
    for (; outer != outer_stop; outer += outer_.step) {
      for (const double* inner = outer+offset_inner+inner_.start;
          inner != outer+offset_inner+inner_.stop;
          inner += inner_.step) {
        *odata++ = *inner;
      }
    }
    return 0;
  }

  int GetNonzerosParamVector::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // Parametric index -> any input disturbance propagates to any output
    bvec_t a = bvec_or(arg[0], dep(0).nnz());
    bvec_t *r = res[0];
    std::fill(r, r+nnz(), a);
    return 0;
  }

  int GetNonzerosParamVector::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *a = arg[0];
    bvec_t r = bvec_or(res[0], nnz());
    std::fill(res[0], res[0]+nnz(), bvec_t(0));

    for (casadi_int i=0;i<dep(0).nnz();++i) {
      *a++ |= r;
    }
    return 0;
  }

  int GetNonzerosParamSlice::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    /**const bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (casadi_int k=s_.start; k!=s_.stop; k+=s_.step) {
      *r++ = a[k];
    }*/
    return 0;
  }

  int GetNonzerosParamSlice::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    /*bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (casadi_int k=s_.start; k!=s_.stop; k+=s_.step) {
      a[k] |= *r;
      *r++ = 0;
    }*/
    return 0;
  }

  int GetNonzerosParamSlice2::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    /*const bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (casadi_int k1=outer_.start; k1!=outer_.stop; k1+=outer_.step) {
      for (casadi_int k2=k1+inner_.start; k2!=k1+inner_.stop; k2+=inner_.step) {
        *r++ = a[k2];
      }
    }*/
    return 0;
  }

  int GetNonzerosParamSlice2::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    /*bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (casadi_int k1=outer_.start; k1!=outer_.stop; k1+=outer_.step) {
      for (casadi_int k2=k1+inner_.start; k2!=k1+inner_.stop; k2+=inner_.step) {
        a[k2] |= *r;
        *r++ = 0;
      }
    }*/
    return 0;
  }

  std::string GetNonzerosParamVector::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << "[" << arg.at(1) << "]";
    return ss.str();
  }

  std::string GetNonzerosParamSlice::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << "[" << arg.at(1) << "+" << s_ << "]";
    return ss.str();
  }

  std::string GetNonzerosParamSlice2::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << "[" << arg.at(1) << "+(" << outer_ << ";" << inner_ << ")]";
    return ss.str();
  }

  void GetNonzerosParam::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = project(arg[0], dep(0).sparsity())->get_nz_ref(arg[1]);
  }

  void GetNonzerosParam::ad_forward(const std::vector<std::vector<MX> >& fseed,
                            std::vector<std::vector<MX> >& fsens) const {
    const MX& nz = dep(1);
    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<fsens.size(); ++d) {
      // Get references to arguments and results
      MX arg = project(fseed[d][0], dep(0).sparsity());
      fsens[d][0] = arg->get_nz_ref(nz);
    }
  }

  void GetNonzerosParam::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                            std::vector<std::vector<MX> >& asens) const {
    const MX& nz = dep(1);
    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<asens.size(); ++d) {
      // Get references to arguments and results
      MX arg = project(aseed[d][0], sparsity());
      asens[d][0] += arg->get_nzadd(DM::zeros(dep(0).sparsity()), nz);
    }

  }

  void GetNonzerosParamVector::generate(CodeGenerator& g,
                                    const std::vector<casadi_int>& arg,
                                    const std::vector<casadi_int>& res) const {
    g.local("i", "casadi_int");
    g.local("rr", "casadi_real", "*");
    g.local("cr", "const casadi_real", "*");
    g << "for (rr=" << g.work(res[0], nnz()) << ", cr=" << g.work(arg[1], dep(1).nnz())
      << "; rr!=" << g.work(res[0], nnz()) << "+" << nnz()
      << "; ++rr) { i=(int) *cr++; "
      << "*rr = i>=0 && i<" << dep(0).nnz() << " ? "
      << g.work(arg[0], dep(0).nnz()) <<  "[i] : " << g.constant(nan) << "; }\n";
  }

  void GetNonzerosParamSlice::generate(CodeGenerator& g,
                                  const std::vector<casadi_int>& arg,
                                  const std::vector<casadi_int>& res) const {

  }

  void GetNonzerosParamSlice2::generate(CodeGenerator& g,
                                    const std::vector<casadi_int>& arg,
                                    const std::vector<casadi_int>& res) const {
  }

  void GetNonzerosParamVector::serialize_body(SerializingStream& s) const {
    GetNonzerosParam::serialize_body(s);
  }

  void GetNonzerosParamVector::serialize_type(SerializingStream& s) const {
    GetNonzerosParam::serialize_type(s);
    s.pack("GetNonzerosParam::type", 'a');
  }

  GetNonzerosParamVector::GetNonzerosParamVector(DeserializingStream& s) : GetNonzerosParam(s) {

  }

  void GetNonzerosParamSlice::serialize_body(SerializingStream& s) const {
    GetNonzerosParam::serialize_body(s);
    s.pack("GetNonzerosParamSlice::slice", s_);
  }

  void GetNonzerosParamSlice::serialize_type(SerializingStream& s) const {
    GetNonzerosParam::serialize_type(s);
    s.pack("GetNonzerosParam::type", 'b');
  }

  GetNonzerosParamSlice::GetNonzerosParamSlice(DeserializingStream& s) : GetNonzerosParam(s) {
    s.unpack("GetNonzerosParamSlice::slice", s_);
  }

  void GetNonzerosParamSlice2::serialize_body(SerializingStream& s) const {
    GetNonzerosParam::serialize_body(s);
    s.pack("GetNonzerosParamSlice2::inner", inner_);
    s.pack("GetNonzerosParamSlice2::outer", outer_);
  }

  void GetNonzerosParamSlice2::serialize_type(SerializingStream& s) const {
    GetNonzerosParam::serialize_type(s);
    s.pack("GetNonzerosParam::type", 'c');
  }

  GetNonzerosParamSlice2::GetNonzerosParamSlice2(DeserializingStream& s) : GetNonzerosParam(s) {
    s.unpack("GetNonzerosParamVector2::inner", inner_);
    s.unpack("GetNonzerosParamVector2::outer", outer_);
  }

  MXNode* GetNonzerosParam::deserialize(DeserializingStream& s) {
    char t;
    s.unpack("GetNonzerosParam::type", t);
    switch (t) {
      case 'a': return new GetNonzerosParamVector(s);
      case 'b': return new GetNonzerosParamSlice(s);
      case 'c': return new GetNonzerosParamSlice2(s);
      default: casadi_assert_dev(false);
    }
  }

} // namespace casadi
