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
  MX SetNonzerosParam<Add>::create(const MX& y, const MX& x, const Slice& inner, const MX& outer) {
    casadi_assert(outer.is_vector() && outer.is_dense(), "outer must be dense vector");
    return MX::create(new SetNonzerosSliceParam<Add>(y, x, inner, outer));
  }

  template<bool Add>
  MX SetNonzerosParam<Add>::create(const MX& y, const MX& x, const MX& inner, const Slice& outer) {
    casadi_assert(inner.is_vector() && inner.is_dense(), "inner must be dense vector");
    return MX::create(new SetNonzerosParamSlice<Add>(y, x, inner, outer));
  }

  template<bool Add>
  MX SetNonzerosParam<Add>::create(const MX& y, const MX& x, const MX& inner, const MX& outer) {
    casadi_assert(inner.is_vector() && inner.is_dense(), "inner must be dense vector");
    casadi_assert(outer.is_vector() && outer.is_dense(), "outer must be dense vector");
    return MX::create(new SetNonzerosParamParam<Add>(y, x, inner, outer));
  }

  template<bool Add>
  SetNonzerosParam<Add>::SetNonzerosParam(const MX& y, const MX& x, const MX& nz) {
    this->set_sparsity(y.sparsity());
    this->set_dep(y, x, nz);
  }

  template<bool Add>
  SetNonzerosParam<Add>::SetNonzerosParam(const MX& y, const MX& x, const MX& nz, const MX& nz2) {
    this->set_sparsity(y.sparsity());
    this->set_dep({y, x, nz, nz2});
  }

  template<bool Add>
  SetNonzerosParamVector<Add>::SetNonzerosParamVector(const MX& y, const MX& x,
      const MX& nz) : SetNonzerosParam<Add>(y, x, nz) {
  }

  template<bool Add>
  SetNonzerosParam<Add>:: ~SetNonzerosParam() {
  }

  template<bool Add>
  void SetNonzerosParamVector<Add>::
  eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    // Add to the element to the sensitivity, if any
    MX arg0 = project(arg[0], this->dep(0).sparsity());
    MX arg1 = project(arg[1], this->dep(1).sparsity());
    MX nz = arg[2];
    if (Add) {
      res[0] = arg1->get_nzadd(arg0, nz);
    } else {
      res[0] = arg1->get_nzassign(arg0, nz);
    }
  }

  template<bool Add>
  void SetNonzerosParamSlice<Add>::
  eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    // Add to the element to the sensitivity, if any
    MX arg0 = project(arg[0], this->dep(0).sparsity());
    MX arg1 = project(arg[1], this->dep(1).sparsity());
    MX inner = arg[2];
    if (Add) {
      res[0] = arg1->get_nzadd(arg0, inner, outer_);
    } else {
      res[0] = arg1->get_nzassign(arg0,  inner, outer_);
    }
  }

  template<bool Add>
  void SetNonzerosSliceParam<Add>::
  eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    // Add to the element to the sensitivity, if any
    MX arg0 = project(arg[0], this->dep(0).sparsity());
    MX arg1 = project(arg[1], this->dep(1).sparsity());
    MX outer = arg[2];
    if (Add) {
      res[0] = arg1->get_nzadd(arg0, inner_, outer);
    } else {
      res[0] = arg1->get_nzassign(arg0, inner_, outer);
    }
  }

  template<bool Add>
  void SetNonzerosParamParam<Add>::
  eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    // Add to the element to the sensitivity, if any
    MX arg0 = project(arg[0], this->dep(0).sparsity());
    MX arg1 = project(arg[1], this->dep(1).sparsity());
    MX inner = arg[2];
    MX outer = arg[2];
    if (Add) {
      res[0] = arg1->get_nzadd(arg0, inner, outer);
    } else {
      res[0] = arg1->get_nzassign(arg0, inner, outer);
    }
  }

  template<bool Add>
  void SetNonzerosParamVector<Add>::ad_forward(const std::vector<std::vector<MX> >& fseed,
                                 std::vector<std::vector<MX> >& fsens) const {
    const MX& nz = this->dep(2);
    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<fsens.size(); ++d) {

      // Get references to arguments and results
      MX arg0 = project(fseed[d][0], this->dep(0).sparsity());
      MX arg1 = project(fseed[d][1], this->dep(1).sparsity());

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
  void SetNonzerosParamVector<Add>::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                                 std::vector<std::vector<MX> >& asens) const {
    const MX& nz = this->dep(2);
    for (casadi_int d=0; d<aseed.size(); ++d) {
      MX seed = project(aseed[d][0], this->sparsity());

      /*
        dep(0) <-> y
        dep(1) <-> x

        z: y[nz]+=x

        bar(x) += bar(z)[nz]
        bar(y) += bar(z)
      */
      asens[d][1] += seed->get_nz_ref(nz);
      if (!Add) {
        asens[d][0] += MX::zeros(this->dep(1).sparsity())->get_nzassign(seed, nz);
      } else {
        asens[d][0] += seed;
      }
    }
  }


  template<bool Add>
  void SetNonzerosParamSlice<Add>::ad_forward(const std::vector<std::vector<MX> >& fseed,
                                 std::vector<std::vector<MX> >& fsens) const {
    const MX& inner = this->dep(2);
    for (casadi_int d=0; d<fsens.size(); ++d) {
      MX arg0 = project(fseed[d][0], this->dep(0).sparsity());
      MX arg1 = project(fseed[d][1], this->dep(1).sparsity());

      MX& res = fsens[d][0];
      res = arg0;

      if (Add) {
        res = arg1->get_nzadd(res, inner, outer_);
      } else {
        res = arg1->get_nzassign(res, inner, outer_);
      }
    }
  }

  template<bool Add>
  void SetNonzerosParamSlice<Add>::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                                 std::vector<std::vector<MX> >& asens) const {
    const MX& inner = this->dep(2);
    for (casadi_int d=0; d<aseed.size(); ++d) {
      MX seed = project(aseed[d][0], this->sparsity());
      asens[d][1] += seed->get_nz_ref(inner, outer_);
      if (!Add) {
        asens[d][0] += MX::zeros(this->dep(1).sparsity())->get_nzassign(seed, inner, outer_);
      } else {
        asens[d][0] += seed;
      }
    }
  }

  template<bool Add>
  void SetNonzerosSliceParam<Add>::ad_forward(const std::vector<std::vector<MX> >& fseed,
                                 std::vector<std::vector<MX> >& fsens) const {
    const MX& outer = this->dep(2);
    for (casadi_int d=0; d<fsens.size(); ++d) {
      MX arg0 = project(fseed[d][0], this->dep(0).sparsity());
      MX arg1 = project(fseed[d][1], this->dep(1).sparsity());

      MX& res = fsens[d][0];
      res = arg0;

      if (Add) {
        res = arg1->get_nzadd(res, inner_, outer);
      } else {
        res = arg1->get_nzassign(res, inner_, outer);
      }
    }
  }

  template<bool Add>
  void SetNonzerosSliceParam<Add>::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                                 std::vector<std::vector<MX> >& asens) const {
    const MX& outer = this->dep(2);
    for (casadi_int d=0; d<aseed.size(); ++d) {
      MX seed = project(aseed[d][0], this->sparsity());
      asens[d][1] += seed->get_nz_ref(inner_, outer);
      if (!Add) {
        asens[d][0] += MX::zeros(this->dep(1).sparsity())->get_nzassign(seed, inner_, outer);
      } else {
        asens[d][0] += seed;
      }
    }
  }

  template<bool Add>
  void SetNonzerosParamParam<Add>::ad_forward(const std::vector<std::vector<MX> >& fseed,
                                 std::vector<std::vector<MX> >& fsens) const {
    const MX& inner = this->dep(2);
    const MX& outer = this->dep(3);
    for (casadi_int d=0; d<fsens.size(); ++d) {
      MX arg0 = project(fseed[d][0], this->dep(0).sparsity());
      MX arg1 = project(fseed[d][1], this->dep(1).sparsity());

      MX& res = fsens[d][0];
      res = arg0;

      if (Add) {
        res = arg1->get_nzadd(res, inner, outer);
      } else {
        res = arg1->get_nzassign(res, inner, outer);
      }
    }
  }

  template<bool Add>
  void SetNonzerosParamParam<Add>::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                                 std::vector<std::vector<MX> >& asens) const {
    const MX& inner = this->dep(2);
    const MX& outer = this->dep(3);
    for (casadi_int d=0; d<aseed.size(); ++d) {
      MX seed = project(aseed[d][0], this->sparsity());
      asens[d][1] += seed->get_nz_ref(inner, outer);
      if (!Add) {
        asens[d][0] += MX::zeros(this->dep(1).sparsity())->get_nzassign(seed, inner, outer);
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
    casadi_int nnz = this->dep(2).nnz();
    casadi_int max_ind = this->dep(0).nnz();
    if (idata0 != odata) {
      copy(idata0, idata0+this->dep(0).nnz(), odata);
    }
    for (casadi_int k=0; k<nnz; ++k) {
      // Get index
      casadi_int index = static_cast<casadi_int>(*nz++);
      if (Add) {
        if (index>=0 && index<max_ind) odata[index] += *idata;
      } else {
        if (index>=0 && index<max_ind) odata[index] = *idata;
      }
      idata++;
    }
    return 0;
  }

  template<bool Add>
  int SetNonzerosParamSlice<Add>::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    const double* idata0 = arg[0];
    const double* idata = arg[1];
    const double* nz = arg[2];
    double* odata = res[0];
    // Dimensions
    casadi_int nnz = this->dep(2).nnz();
    casadi_int max_ind = this->dep(0).nnz();
    if (idata0 != odata) {
      copy(idata0, idata0+this->dep(0).nnz(), odata);
    }

    casadi_int* inner = iw; iw += nnz;
    for (casadi_int i=0; i<nnz; ++i) {
      // Get index
      inner[i] = static_cast<casadi_int>(*nz++);
    }
    for (casadi_int i=outer_.start;i<outer_.stop;i+= outer_.step) {
      // Get index
      for (casadi_int* inner_it=inner; inner_it!=inner+nnz; ++inner_it) {
        casadi_int index = i+*inner_it;
        if (Add) {
          if (index>=0 && index<max_ind) odata[index] += *idata;
        } else {
          if (index>=0 && index<max_ind) odata[index] = *idata;
        }
        idata++;
      }
    }
    return 0;
  }

  template<bool Add>
  int SetNonzerosSliceParam<Add>::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    const double* idata0 = arg[0];
    const double* idata = arg[1];
    const double* nz = arg[2];
    double* odata = res[0];
    // Dimensions
    casadi_int nnz = this->dep(2).nnz();
    casadi_int max_ind = this->dep(0).nnz();
    if (idata0 != odata) {
      copy(idata0, idata0+this->dep(0).nnz(), odata);
    }
    for (casadi_int k=0; k<nnz; ++k) {
      // Get index
      casadi_int ind = static_cast<casadi_int>(*nz++);
      for (casadi_int j=0;j<inner_.stop;j+= inner_.step) {
        casadi_int index = ind+j;
        if (Add) {
          if (index>=0 && index<max_ind) odata[index] += *idata;
        } else {
          if (index>=0 && index<max_ind) odata[index] = *idata;
        }
        idata++;
      }
    }
    return 0;
  }

  template<bool Add>
  int SetNonzerosParamParam<Add>::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
  const double* idata0 = arg[0];
    const double* idata = arg[1];
    const double* nz = arg[2];
    const double* nz2 = arg[3];
    double* odata = res[0];
    // Dimensions
    casadi_int nnz = this->dep(2).nnz();
    casadi_int nnz2 = this->dep(3).nnz();
    casadi_int max_ind = this->dep(0).nnz();
    if (idata0 != odata) {
      copy(idata0, idata0+this->dep(0).nnz(), odata);
    }

    casadi_int* inner = iw; iw += nnz;
    for (casadi_int i=0; i<nnz; ++i) {
      // Get index
      inner[i] = static_cast<casadi_int>(*nz++);
    }

    for (casadi_int k=0; k<nnz2; ++k) {
      // Get index
      casadi_int ind = static_cast<casadi_int>(*nz2++);
      for (casadi_int* inner_it=inner; inner_it!=inner+nnz; ++inner_it) {
        casadi_int index = ind+*inner_it;
        if (Add) {
          if (index>=0 && index<max_ind) odata[index] += *idata;
        } else {
          if (index>=0 && index<max_ind) odata[index] = *idata;
        }
        idata++;
      }
    }
    return 0;
  }

  template<bool Add>
  size_t SetNonzerosParamSlice<Add>::
  sz_iw() const {
    return this->dep(2).nnz();
  }

  template<bool Add>
  size_t SetNonzerosParamParam<Add>::
  sz_iw() const {
    return this->dep(2).nnz();
  }


  template<bool Add>
  int SetNonzerosParam<Add>::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // Parametric index -> any input disturbance propagates to any output
    bvec_t arg0 = bvec_or(arg[0], this->dep(0).nnz());
    bvec_t arg1 = bvec_or(arg[1], this->dep(1).nnz());

    bvec_t *r = res[0];
    std::fill(r, r+this->nnz(), arg0 | arg1);
    return 0;
  }

  template<bool Add>
  int SetNonzerosParam<Add>::
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
  std::string SetNonzerosParamVector<Add>::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << "(" << arg.at(0) << "[" << arg.at(2) << "]";
    ss << (Add ? " += " : " = ") << arg.at(1) << ")";
    return ss.str();
  }

  template<bool Add>
  std::string SetNonzerosParamSlice<Add>::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << "(" << arg.at(0) << "[(" << arg.at(2) << ";" << outer_ << ")]";
    ss << (Add ? " += " : " = ") << arg.at(1) << ")";
    return ss.str();
  }

  template<bool Add>
  std::string SetNonzerosSliceParam<Add>::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << "(" << arg.at(0) << "[(" << inner_ << ";" << arg.at(2) << ")]";
    ss << (Add ? " += " : " = ") << arg.at(1) << ")";
    return ss.str();
  }

  template<bool Add>
  std::string SetNonzerosParamParam<Add>::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << "(" << arg.at(0) << "[(" << arg.at(2) << ";" << arg.at(3) << ")]";
    ss << (Add ? " += " : " = ") << arg.at(1) << ")";
    return ss.str();
  }

  template<bool Add>
  void SetNonzerosParam<Add>::
  generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {
    // Copy first argument if not inplace
    if (arg[0]!=res[0]) {
      g << g.copy(g.work(arg[0], this->dep(0).nnz()), this->nnz(),
                          g.work(res[0], this->nnz())) << '\n';
    }
  }

  template<bool Add>
  void SetNonzerosParamVector<Add>::
  generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {
    SetNonzerosParam<Add>::generate(g, arg, res);

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
    SetNonzerosParam<Add>::generate(g, arg, res);

    casadi_int n = this->dep(1).nnz();
    casadi_int n_inner = this->dep(2).nnz();

    g.local("cii", "const casadi_int", "*");
    g.local("i", "casadi_int");
    g << "for (i=0;i<" << n_inner << ";++i) iw[i] = (int) "
      << g.work(arg[2], n_inner) << "[i];\n";

    g.local("cs", "const casadi_real", "*");
    g.local("k", "casadi_int");
    g << "for (cs=" << g.work(arg[1], n)
      << ", k=" << outer_.start << ";k<" << outer_.stop << ";k+=" << outer_.step << ") ";
    g << "for (cii=iw; cii!=iw" << "+" << n_inner << "; ++cii) { i=k+*cii; "
      << "if (i>=0 && i<" << this->dep(0).nnz() << ") "
      << g.work(res[0], this->nnz()) << "[i] " << (Add?"+= ":"= ")
      << "*cs; cs++; }\n";
  }

  template<bool Add>
  void SetNonzerosSliceParam<Add>::
  generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {
    SetNonzerosParam<Add>::generate(g, arg, res);

    casadi_int n = this->dep(1).nnz();
    casadi_int n_outer = this->dep(2).nnz();

    g.local("i", "casadi_int");
    g.local("j", "casadi_int");
    g.local("k", "casadi_int");
    g.local("cr", "const casadi_real", "*");
    g.local("cs", "const casadi_real", "*");
    g << "for (cr=" << g.work(arg[2], n_outer)
      << ", cs=" << g.work(arg[1], n)
      << "; cr!=" << g.work(arg[2], n_outer) << "+" << n_outer
      << "; ++cr) ";
    g << "for (j=(int) *cr, "
      << "k=" << inner_.start << ";k<" << inner_.stop << ";k+=" << inner_.step << ") ";
    g << "{ i=k+j; "
      << "if (i>=0 && i<" << this->dep(0).nnz() << ") "
      << g.work(res[0], this->nnz()) << "[i] " << (Add?"+= ":"= ")
      << "*cs; cs++; }\n";
  }

  template<bool Add>
  void SetNonzerosParamParam<Add>::
  generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {
    SetNonzerosParam<Add>::generate(g, arg, res);
    casadi_int n = this->dep(1).nnz();
    casadi_int n_outer = this->dep(3).nnz();
    casadi_int n_inner = this->dep(2).nnz();

    g.local("cii", "const casadi_int", "*");
    g.local("i", "casadi_int");
    g << "for (i=0;i<" << n_inner << ";++i) iw[i] = (int) "
      << g.work(arg[2], n_inner) << "[i];\n";

    g.local("j", "casadi_int");
    g.local("cr", "const casadi_real", "*");
    g.local("cs", "const casadi_real", "*");
    g << "for (cr=" << g.work(arg[3], n_outer)
      << ", cs=" << g.work(arg[1], n)
      << "; cr!=" << g.work(arg[3], n_outer) << "+" << n_outer
      << "; ++cr) ";
    g << "for (j=(int) *cr, cii=iw; cii!=iw" << "+" << n_inner << "; ++cii) { i=j+*cii; "
      << "if (i>=0 && i<" << this->dep(0).nnz() << ") "
      << g.work(res[0], this->nnz()) << "[i] " << (Add?"+= ":"= ")
      << "*cs; cs++; }\n";
  }

  template<bool Add>
  void SetNonzerosParamVector<Add>::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
  }

  template<bool Add>
  SetNonzerosParamVector<Add>::SetNonzerosParamVector(DeserializingStream& s) :
      SetNonzerosParam<Add>(s) {
  }

  template<bool Add>
  void SetNonzerosParamVector<Add>::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("SetNonzerosParam::type", 'a');
  }

  template<bool Add>
  void SetNonzerosParamSlice<Add>::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("SetNonzerosParamSlice::outer", outer_);
  }

  template<bool Add>
  SetNonzerosParamSlice<Add>::SetNonzerosParamSlice(DeserializingStream& s) :
      SetNonzerosParam<Add>(s) {
    s.unpack("SetNonzerosParamSlice::outer", outer_);
  }

  template<bool Add>
  void SetNonzerosParamSlice<Add>::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("SetNonzerosParam::type", 'b');
  }

  template<bool Add>
  void SetNonzerosSliceParam<Add>::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("SetNonzerosSliceParam::inner", inner_);
  }

  template<bool Add>
  SetNonzerosSliceParam<Add>::SetNonzerosSliceParam(DeserializingStream& s) :
      SetNonzerosParam<Add>(s) {
    s.unpack("SetNonzerosSliceParam::inner", inner_);
  }

  template<bool Add>
  void SetNonzerosSliceParam<Add>::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("SetNonzerosParam::type", 'c');
  }


  template<bool Add>
  void SetNonzerosParamParam<Add>::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("SetNonzerosParam::type", 'd');
  }

  template<bool Add>
  SetNonzerosParamParam<Add>::SetNonzerosParamParam(DeserializingStream& s) :
      SetNonzerosParam<Add>(s) {
  }

  template<bool Add>
  MXNode* SetNonzerosParam<Add>::deserialize(DeserializingStream& s) {
    char t;
    s.unpack("SetNonzerosParam::type", t);
    switch (t) {
      case 'a': return new SetNonzerosParamVector<Add>(s);
      case 'b': return new SetNonzerosParamSlice<Add>(s);
      case 'c': return new SetNonzerosSliceParam<Add>(s);
      case 'd': return new SetNonzerosParamParam<Add>(s);
      default: casadi_assert_dev(false);
    }
  }

} // namespace casadi

/// \endcond

#endif // CASADI_SETNONZEROS_PARAM_IMPL_HPP
