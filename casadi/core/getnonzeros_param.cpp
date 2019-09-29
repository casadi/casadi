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

  MX GetNonzerosParam::create(const MX& x, const MX& inner, const Slice& outer) {
    casadi_assert(inner.is_vector() && inner.is_dense(), "inner must be dense vector");
    return MX::create(new GetNonzerosParamSlice(
      Sparsity::dense(inner.numel(), outer.size()), x, inner, outer));
  }

  MX GetNonzerosParam::create(const MX& x, const Slice& inner, const MX& outer) {
    casadi_assert(outer.is_vector() && outer.is_dense(), "outer must be dense vector");
    return MX::create(new GetNonzerosSliceParam(
      Sparsity::dense(inner.size(), outer.numel()), x, inner, outer));
  }

  MX GetNonzerosParam::create(const MX& x, const MX& inner, const MX& outer) {
    casadi_assert(outer.is_vector() && outer.is_dense(), "outer must be dense vector");
    casadi_assert(inner.is_vector() && inner.is_dense(), "inner must be dense vector");
    return MX::create(new GetNonzerosParamParam(
      Sparsity::dense(inner.numel(), outer.numel()), x, inner, outer));
  }

  GetNonzerosParam::GetNonzerosParam(const Sparsity& sp, const MX& y, const MX& nz) {
    set_sparsity(sp);
    set_dep(y, nz);
  }

  GetNonzerosParam::GetNonzerosParam(const Sparsity& sp, const MX& y,
      const MX& nz, const MX& nz_extra) {
    set_sparsity(sp);
    set_dep(y, nz, nz_extra);
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

  size_t GetNonzerosParamSlice::sz_iw() const {
    return dep(1).nnz();
  }

  int GetNonzerosParamSlice::eval(const double** arg, double** res,
                                 casadi_int* iw, double* w) const {
    const double* idata = arg[0];
    const double* nz = arg[1];
    double* odata = res[0];

    // Dimensions
    casadi_int nnz = dep(1).nnz();
    casadi_int max_ind = dep(0).nnz();

    casadi_int* inner = iw; iw += nnz;
    for (casadi_int i=0; i<nnz; ++i) {
      // Get index
      inner[i] = static_cast<casadi_int>(*nz++);
    }
    // Get elements
    for (casadi_int i=outer_.start;i<outer_.stop;i+= outer_.step) {
      // Get index
      for (casadi_int* inner_it=inner; inner_it!=inner+nnz; ++inner_it) {
        casadi_int index = i+*inner_it;
        // Make assignment if in bounds, else NaN
        *odata++ = index>=0 && index<max_ind ? idata[index] : nan;
      }
    }
    return 0;
  }

  int GetNonzerosSliceParam::
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
      casadi_int ind = static_cast<casadi_int>(*nz++);
      for (casadi_int j=inner_.start;j<inner_.stop;j+= inner_.step) {
        casadi_int index = ind+j;
        // Make assignment if in bounds, else NaN
        *odata++ = index>=0 && index<max_ind ? idata[index] : nan;
      }
    }
    return 0;
  }

  size_t GetNonzerosParamParam::sz_iw() const {
    return dep(1).nnz();
  }

  int GetNonzerosParamParam::eval(const double** arg, double** res,
                                 casadi_int* iw, double* w) const {
    const double* idata = arg[0];
    const double* nz = arg[1];
    const double* nz2 = arg[2];
    double* odata = res[0];


    // Dimensions
    casadi_int nnz = dep(1).nnz();
    casadi_int nnz2 = dep(2).nnz();
    casadi_int max_ind = dep(0).nnz();

    casadi_int* inner = iw; iw += nnz;
    for (casadi_int i=0; i<nnz; ++i) {
      // Get index
      inner[i] = static_cast<casadi_int>(*nz++);
    }
    for (casadi_int i=0; i<nnz2; ++i) {
      // Get index
      casadi_int ind = static_cast<casadi_int>(*nz2++);
      // Get index
      for (casadi_int* inner_it=inner; inner_it!=inner+nnz; ++inner_it) {
        casadi_int index = ind+*inner_it;
        // Make assignment if in bounds, else NaN
        *odata++ = index>=0 && index<max_ind ? idata[index] : nan;
      }
    }
    return 0;
  }

  int GetNonzerosParam::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // Parametric index -> any input disturbance propagates to any output
    bvec_t a = bvec_or(arg[0], dep(0).nnz());
    bvec_t *r = res[0];
    std::fill(r, r+nnz(), a);
    return 0;
  }

  int GetNonzerosParam::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *a = arg[0];
    bvec_t r = bvec_or(res[0], nnz());
    std::fill(res[0], res[0]+nnz(), bvec_t(0));

    for (casadi_int i=0;i<dep(0).nnz();++i) {
      *a++ |= r;
    }
    return 0;
  }

  std::string GetNonzerosParamVector::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << "[" << arg.at(1) << "]";
    return ss.str();
  }

  std::string GetNonzerosParamSlice::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << "[(" << arg.at(1) << ";" << outer_ << ")]";
    return ss.str();
  }

  std::string GetNonzerosSliceParam::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << "[(" << inner_ << ";" << arg.at(1) << ")]";
    return ss.str();
  }

  std::string GetNonzerosParamParam::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << "[(" << arg.at(1) << ";" << arg.at(2) << ")]";
    return ss.str();
  }

  void GetNonzerosParamVector::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = project(arg[0], dep(0).sparsity())->get_nz_ref(arg[1]);
  }

  void GetNonzerosParamSlice::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = project(arg[0], dep(0).sparsity())->get_nz_ref(arg[1], outer_);
  }

  void GetNonzerosSliceParam::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = project(arg[0], dep(0).sparsity())->get_nz_ref(inner_, arg[1]);
  }

  void GetNonzerosParamParam::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = project(arg[0], dep(0).sparsity())->get_nz_ref(arg[1], arg[2]);
  }

  void GetNonzerosParamVector::ad_forward(const std::vector<std::vector<MX> >& fseed,
                            std::vector<std::vector<MX> >& fsens) const {
    const MX& nz = dep(1);
    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<fsens.size(); ++d) {
      // Get references to arguments and results
      MX arg = project(fseed[d][0], dep(0).sparsity());
      fsens[d][0] = arg->get_nz_ref(nz);
    }
  }

  void GetNonzerosParamVector::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                            std::vector<std::vector<MX> >& asens) const {
    const MX& nz = dep(1);
    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<asens.size(); ++d) {
      // Get references to arguments and results
      MX arg = project(aseed[d][0], sparsity());
      asens[d][0] += arg->get_nzadd(DM::zeros(dep(0).sparsity()), nz);
    }
  }

  void GetNonzerosParamSlice::ad_forward(const std::vector<std::vector<MX> >& fseed,
                            std::vector<std::vector<MX> >& fsens) const {
    const MX& inner = dep(1);
    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<fsens.size(); ++d) {
      // Get references to arguments and results
      MX arg = project(fseed[d][0], dep(0).sparsity());
      fsens[d][0] = arg->get_nz_ref(inner, outer_);
    }
  }

  void GetNonzerosParamSlice::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                            std::vector<std::vector<MX> >& asens) const {
    const MX& inner = dep(1);
    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<asens.size(); ++d) {
      // Get references to arguments and results
      MX arg = project(aseed[d][0], sparsity());
      asens[d][0] += arg->get_nzadd(DM::zeros(dep(0).sparsity()), inner, outer_);
    }
  }


  void GetNonzerosSliceParam::ad_forward(const std::vector<std::vector<MX> >& fseed,
                            std::vector<std::vector<MX> >& fsens) const {
    const MX& outer = dep(1);
    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<fsens.size(); ++d) {
      // Get references to arguments and results
      MX arg = project(fseed[d][0], dep(0).sparsity());
      fsens[d][0] = arg->get_nz_ref(inner_, outer);
    }
  }

  void GetNonzerosSliceParam::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                            std::vector<std::vector<MX> >& asens) const {
    const MX& outer = dep(1);
    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<asens.size(); ++d) {
      // Get references to arguments and results
      MX arg = project(aseed[d][0], sparsity());
      asens[d][0] += arg->get_nzadd(DM::zeros(dep(0).sparsity()), inner_, outer);
    }

  }

  void GetNonzerosParamParam::ad_forward(const std::vector<std::vector<MX> >& fseed,
                            std::vector<std::vector<MX> >& fsens) const {
    const MX& inner = dep(1);
    const MX& outer = dep(2);
    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<fsens.size(); ++d) {
      // Get references to arguments and results
      MX arg = project(fseed[d][0], dep(0).sparsity());
      fsens[d][0] = arg->get_nz_ref(inner, outer);
    }
  }

  void GetNonzerosParamParam::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                            std::vector<std::vector<MX> >& asens) const {
    const MX& inner = dep(1);
    const MX& outer = dep(2);
    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<asens.size(); ++d) {
      // Get references to arguments and results
      MX arg = project(aseed[d][0], sparsity());
      asens[d][0] += arg->get_nzadd(DM::zeros(dep(0).sparsity()), inner, outer);
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
    g.local("cii", "const casadi_int", "*");
    g.local("i", "casadi_int");
    g << "for (i=0;i<" << dep(1).nnz() << ";++i) iw[i] = (int) "
      << g.work(arg[1], dep(1).nnz()) << "[i];\n";

    g.local("rr", "casadi_real", "*");
    g.local("k", "casadi_int");
    g << "for (rr=" << g.work(res[0], nnz()) << ", "
      << "k=" << outer_.start << ";k<" << outer_.stop << ";k+=" << outer_.step << ") ";
    g << "for (cii=iw; cii!=iw" << "+" << dep(1).nnz() << "; ++cii) { i=k+*cii; "
      << "*rr++ = i>=0 && i<" << dep(0).nnz() << " ? "
      << g.work(arg[0], dep(0).nnz()) <<  "[i] : " << g.constant(nan) << "; }\n";
  }

  void GetNonzerosSliceParam::generate(CodeGenerator& g,
                                    const std::vector<casadi_int>& arg,
                                    const std::vector<casadi_int>& res) const {
    g.local("i", "casadi_int");
    g.local("j", "casadi_int");
    g.local("rr", "casadi_real", "*");
    g.local("k", "casadi_int");
    g.local("cr", "const casadi_real", "*");
    g << "for (cr=" << g.work(arg[1], dep(1).nnz())
      << ", rr=" << g.work(res[0], nnz())
      << "; cr!=" << g.work(arg[1], dep(1).nnz()) << "+" << dep(1).nnz()
      << "; ++cr) ";
    g << "for (j=(int) *cr, "
      << "k=" << inner_.start << ";k<" << inner_.stop << ";k+=" << inner_.step << ") ";
    g << "{ i=k+j; "
      << "*rr++ = i>=0 && i<" << dep(0).nnz() << " ? "
      << g.work(arg[0], dep(0).nnz()) <<  "[i] : " << g.constant(nan) << "; }\n";

  }

  void GetNonzerosParamParam::generate(CodeGenerator& g,
                                    const std::vector<casadi_int>& arg,
                                    const std::vector<casadi_int>& res) const {
    g.local("cii", "const casadi_int", "*");
    g.local("i", "casadi_int");
    g << "for (i=0;i<" << dep(1).nnz() << ";++i) iw[i] = (int) "
      << g.work(arg[1], dep(1).nnz()) << "[i];\n";

    g.local("j", "casadi_int");
    g.local("cr", "const casadi_real", "*");
    g.local("rr", "casadi_real", "*");
    g << "for (cr=" << g.work(arg[2], dep(2).nnz())
      << ", rr=" << g.work(res[0], nnz())
      << "; cr!=" << g.work(arg[2], dep(2).nnz()) << "+" << dep(2).nnz()
      << "; ++cr) ";
    g << "for (j=(int) *cr, cii=iw; cii!=iw" << "+" << dep(1).nnz()
      << "; ++cii) ";
    g << "{ i=j+*cii;"
      << "*rr++ = i>=0 && i<" << dep(0).nnz() << " ? "
      << g.work(arg[0], dep(0).nnz()) <<  "[i] : " << g.constant(nan) << "; }\n";
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
    s.pack("GetNonzerosParamSlice::outer", outer_);
  }

  void GetNonzerosParamSlice::serialize_type(SerializingStream& s) const {
    GetNonzerosParam::serialize_type(s);
    s.pack("GetNonzerosParam::type", 'b');
  }

  GetNonzerosParamSlice::GetNonzerosParamSlice(DeserializingStream& s) : GetNonzerosParam(s) {
    s.unpack("GetNonzerosParamSlice::outer", outer_);
  }

  void GetNonzerosSliceParam::serialize_body(SerializingStream& s) const {
    GetNonzerosParam::serialize_body(s);
    s.pack("GetNonzerosSliceParam::inner", inner_);
  }

  void GetNonzerosSliceParam::serialize_type(SerializingStream& s) const {
    GetNonzerosParam::serialize_type(s);
    s.pack("GetNonzerosParam::type", 'c');
  }

  GetNonzerosSliceParam::GetNonzerosSliceParam(DeserializingStream& s) : GetNonzerosParam(s) {
    s.unpack("GetNonzerosSliceParam::inner", inner_);
  }

  void GetNonzerosParamParam::serialize_type(SerializingStream& s) const {
    GetNonzerosParam::serialize_type(s);
    s.pack("GetNonzerosParam::type", 'd');
  }

  GetNonzerosParamParam::GetNonzerosParamParam(DeserializingStream& s) : GetNonzerosParam(s) {
  }

  MXNode* GetNonzerosParam::deserialize(DeserializingStream& s) {
    char t;
    s.unpack("GetNonzerosParam::type", t);
    switch (t) {
      case 'a': return new GetNonzerosParamVector(s);
      case 'b': return new GetNonzerosParamSlice(s);
      case 'c': return new GetNonzerosSliceParam(s);
      case 'd': return new GetNonzerosParamParam(s);
      default: casadi_assert_dev(false);
    }
  }

} // namespace casadi
