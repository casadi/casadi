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


#include "mx_layout.hpp"
#include "runtime/casadi_runtime.hpp"

using namespace std;
namespace casadi {

// Throw informative error message
#define THROW_ERROR(FNAME, WHAT) \
throw CasadiException("Error in PermuteLayout::" FNAME " "\
  "[" + this->class_name() + "] at " + CASADI_WHERE + ":\n"\
  + string(WHAT));


  std::string PermuteLayout::disp(const std::vector<std::string>& arg) const {
    return "permute(" + arg.at(0) + "," + relay_.label_ + ":" +  str(relay_.source()) + "--" + str(relay_.perms()) + "-->" + str(relay_.target()) + ")";
  }

  MX PermuteLayout::create(const MX& x, const Relayout& relay) {
    if (relay.is_trivial()) return x;
    //if (relay.is_trivial()) return x;
    return MX::create(new PermuteLayout(x, relay));
    //casadi_assert(x.is_dense(), "must be dense");
    if (x.is_dense()) {
      return MX::create(new PermuteLayout(x, relay));
    } else {
      return MX::create(new PermuteLayout(densify(x), relay));
    }
  }

  PermuteLayout::PermuteLayout(const MX& x, const Relayout& relay) : relay_(relay) {
    //casadi_assert(x.layout().nnz()==relay.source().nnz(), "impossible relayout");

    set_dep(x);
    set_sparsity(x.sparsity());
    sparsity().spy();
    uout() << "PermuteLayout" << std::endl;
    uout() << relay.perms() << std::endl;
    uout() << relay.source() << std::endl;
    uout() << relay.target() << std::endl;
    casadi_assert(x->nnz()==relay.source().nnz(), "Invalid permutation");
    if (x->sz_self()<relay.source().size()) {

      relay_ = Relayout(Layout(relay_.source().dims()),relay_.perms(), relay_.target(), relay.label());
    }
    casadi_assert(x->sz_self()==relay_.source().size(), "Invalid permutation sz_self" + str(x->sz_self()) + str(relay.source())  + str(relay_.source()) );   
    /*const Layout& source = dep(0).layout();
    if (source.has_padding() || target.has_padding()) {
      set_sparsity(Sparsity::dense(target.is_default() ? source.size() : target.size(), 1));
    } else {
      set_sparsity(x.sparsity());
    }*/

    //uout() << "my sparsity" << sparsity_ << relay.source() << relay.target() << std::endl;
    // Relax to repmat?
    //if (!x.is_dense()) {
    //  sparsity_.spy();
    //  uout() << x << std::endl;
    //}
    //casadi_assert(x.is_dense(), "Sparsity not supported");
    set_layout(relay_.target());
    relay_.source().assert_valid_permutation(relay_.target());
    //uout() << "okay" << std::endl;
  }

  void PermuteLayout::ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const {
    try {
      /*uout() << "forward" << std::endl;
      for (casadi_int d=0; d<fsens.size(); ++d) {
        fsens[d][0] = permute_layout(project(fseed[d][0], sparsity()), relay_);
      }
      uout() << "forward done" << std::endl;
      return;*/
      casadi_int n = fseed.size();
      //uout() << "ad_forward prepre" << std::endl;
      std::vector<MX> seeds;
      for (casadi_int d=0; d<n; ++d) {
        seeds.push_back(fseed[d][0]);//, Relayout(relay_.source(), relay_.target(),"fwdseed"));
        casadi_assert_dev(fseed[d][0]->sz_self()==fseed[d][0]->nnz());
        //uout() << "seed" << std::endl;
        //seeds[d].sparsity().spy();
      }
      //uout() << "ad_forward pre" << std::endl;
      //uout() << "ad_forward pre" << relay_.source().push_right(n) << std::endl;
      //uout() << "ad_forward pre" << relay_.target().push_right(n) << std::endl;
      MX all_seeds = horzcat(seeds);
      Sparsity target_sp = repmat(sparsity(), 1, n);
      casadi_assert((target_sp+all_seeds.sparsity()).nnz()==target_sp.nnz(), "seed sparsity must be subset");
      if (all_seeds.sparsity()!=target_sp) all_seeds = project(all_seeds, target_sp);
      uout() << "here" << std::endl;
      Relayout temp = relay_.push_right(n);
      uout() << temp << std::endl;
      MX sens = permute_layout(all_seeds, relay_.push_right(n));//Relayout(relay_.source().push_right(n), relay_.target().push_right(n), "fwdseed"));
      uout() << "here reached" << std::endl;
      //uout() << "ad_forward" << sens.size() << std::endl;
      std::vector<MX> senss = horzsplit(sens, sens.size2()/n);
      //uout() << "ad_forward post" << senss << std::endl;
      for (casadi_int d=0; d<n; ++d) {
        fsens[d][0] = senss[d];
      }
    } catch (exception& e) {
      THROW_ERROR("ad_forward", e.what());
    }
  }

  void PermuteLayout::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                        std::vector<std::vector<MX> >& asens) const {
    try {
      /*for (casadi_int d=0; d<asens.size(); ++d) {
        asens[d][0] += permute_layout(project(aseed[d][0], sparsity()), relay_.invert());
      }
      return;*/
      casadi_int n = aseed.size();      
      std::vector<MX> seeds;
      for (casadi_int d=0; d<n; ++d) {
        seeds.push_back(aseed[d][0]);
      }
      MX all_seeds = horzcat(seeds);
      Sparsity target_sp = repmat(sparsity(), 1, n);
      casadi_assert((target_sp+all_seeds.sparsity()).nnz()==target_sp.nnz(), "seed sparsity must be subset");
      if (all_seeds.sparsity()!=target_sp) all_seeds = project(all_seeds, target_sp);
      std::string name= "";
      if (relay_.push_right(n).invert().source().size()!=all_seeds->sz_self()) {
        name = "oops";
      }
      Relayout r = relay_.push_right(n).invert();
      MX sens = permute_layout(all_seeds, Relayout(r.source(), r.perms(), r.target(), name));//Relayout(relay_.target().push_right(n), relay_.source().push_right(n), "adjseed"));
      std::vector<MX> senss = horzsplit(sens, sens.size2()/n);
      for (casadi_int d=0; d<aseed.size(); ++d) {
        asens[d][0] += senss[d];
      }
    } catch (exception& e) {
      THROW_ERROR("ad_reverse", e.what());
    }
  }

  MX PermuteLayout::get_nzref(const Sparsity& sp, const vector<casadi_int>& nz) const {
    return as_nzref()->get_nzref(sp, nz);
  }

  MX PermuteLayout::as_nzref() const {
    const casadi_int* source = relay_.source();
    const casadi_int* target = relay_.target();
    const casadi_int* perms = get_ptr(relay_.perms());
    std::vector<casadi_int> arg = range(nnz());
    std::vector<casadi_int> res(nnz());
    std::vector<casadi_int> iw(sz_iw());
    casadi_relayout(get_ptr(arg), get_ptr(res), source, perms, target, get_ptr(iw));
    return dep(0)->get_nzref(sparsity(), res);
  }

  void PermuteLayout::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    //uout() << "eval_mx" << std::endl;
    //s//parsity().spy();
    //uout() << "eval_mx" << arg[0]->sz_self() << std::endl;
    //uout() << "eval_mx" << project(arg[0],sparsity())->sz_self() << std::endl;
    //uout() << "eval_mx" << relay_.source() << std::endl;
    //uout() << "eval_mx" << arg[0] << std::endl;
    //casadi_assert_dev(arg[0]->sz_self()==project(arg[0],sparsity())->sz_self());    
    //arg[0].sparsity().spy();
    res[0] = project(arg[0],sparsity())->get_permute_layout(Relayout(relay_.source(), relay_.perms(), relay_.target(), "eval_mx_"+relay_.label_));
  }

  size_t PermuteLayout::sz_iw() const {
    return relay_.sz_iw();
  }

  int PermuteLayout::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    const casadi_int* source = relay_.source();
    const casadi_int* target = relay_.target();
    const casadi_int* perms = get_ptr(relay_.perms());
    if (arg[0] && res[0]) casadi_relayout(arg[0], res[0], source, perms, target, iw);
    return 0;
  }

  int PermuteLayout::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    const casadi_int* source = relay_.source();
    const casadi_int* target = relay_.target();
    const casadi_int* perms = get_ptr(relay_.perms());
    if (arg[0] && res[0]) casadi_relayout(arg[0], res[0], source, perms, target, iw);
    return 0;
  }

  int PermuteLayout::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    const casadi_int* source = relay_.source();
    const casadi_int* target = relay_.target();
    const casadi_int* perms = get_ptr(relay_.perms());
    if (arg[0] && res[0]) casadi_relayout(arg[0], res[0], source, perms, target, iw);
    return 0;
  }

  template<typename T1>
  void casadi_relayout_rev(T1* arg, T1* res, const casadi_int* source, const casadi_int* perm, const casadi_int* target, casadi_int* iw) {
    casadi_int nnz, i, j, k, kk;
    casadi_int *counter, *target_strides_perm;
    const casadi_int *dims, *target_dims, *strides, *target_strides;
    casadi_int n_dims = source[1];
    nnz = source[2];
    dims = source+3;
    target_dims = target+3;
    strides = dims+2*n_dims;
    target_strides = target_dims+2*n_dims;
    
    counter = iw;
    target_strides_perm = iw+n_dims;
    for (j=0;j<n_dims;++j) {
      target_strides_perm[perm[j]] = target_strides[j];
      counter[j] = 0;
    }
    for (i=0;i<nnz;++i) {
      k = 0;
      for (j=0;j<n_dims;++j) {
        k += strides[j]*counter[j];
      }
      kk = 0;
      for (j=0;j<n_dims;++j) {
        kk += target_strides_perm[j]*counter[j];
      }
      arg[k] |= res[kk];
      res[kk] = 0;
      for (j=0;j<n_dims;++j) {
        counter[j]++;
        if (counter[j]<dims[j]) break;
        counter[j]=0;
      }
    }
  }

  int PermuteLayout::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    const casadi_int* source = relay_.source();
    const casadi_int* target = relay_.target();
    const casadi_int* perms = get_ptr(relay_.perms());
    if (arg[0] && res[0]) casadi_relayout_rev(arg[0], res[0], source, perms, target, iw);
    return 0;
  }

  

  MX PermuteLayout::get_permute_layout(const Relayout& relay) const {
    if (relay_.cancels(relay)) { // cancellation
      return dep(0);
    } else if (relay_.absorbs(relay)) {
      return permute_layout(dep(0), relay_.absorbed(relay));
    } else {
      return MXNode::get_permute_layout(relay);
    }
  }

  void PermuteLayout::generate(CodeGenerator& g,
                        const std::vector<casadi_int>& arg,
                        const std::vector<casadi_int>& res) const {
    relay_.generate(g, g.work(arg[0], nnz()), g.work(res[0], nnz()));

    //g << g.relayout(g.work(arg[0], nnz()), g.work(res[0], nnz()), relay_, "iw") << "\n";
  }

  void PermuteLayout::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("PermuteLayout::relay_source", relay_.source());
    s.pack("PermuteLayout::relay_perms", relay_.perms());
    s.pack("PermuteLayout::relay_target", relay_.target());
  }

  PermuteLayout::PermuteLayout(DeserializingStream& s) : MXNode(s) {
    Layout source, target;
    std::vector<casadi_int> perms;
    s.unpack("PermuteLayout::relay_source", source);
    s.unpack("PermuteLayout::relay_perms", perms);
    s.unpack("PermuteLayout::relay_target", target);
    relay_ = Relayout(source, perms, target);
  }


} // namespace casadi
