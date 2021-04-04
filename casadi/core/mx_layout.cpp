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

  std::string ReinterpretLayout::disp(const std::vector<std::string>& arg) const {
    return "reinterpret(" + arg.at(0) + "," + str(layout_) + ")";
  }

  MX PermuteLayout::create(const MX& x, const Relayout& relay) {
    if (relay.is_trivial()) return x;
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
    set_layout(relay.target());
    relay_.source().assert_valid_permutation(relay_.target());
    //uout() << "okay" << std::endl;
  }

  ReinterpretLayout::ReinterpretLayout(const MX& x, const Layout& target) {
    set_dep(x);
    set_sparsity(x.sparsity());
    // Relax to repmat?
    casadi_assert(x.is_dense(), "Sparsity not supported");
    set_layout(target);
    //uout() << target << target.size() << ":" << target.nnz() << std::endl;
    //casadi_assert(target.size()==x.nnz(), "reinterpret_layout cannot change total amount of data"); issue with padding and keeping sparsity
  }

  void PermuteLayout::ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const {
    try {
      //for (casadi_int d=0; d<fsens.size(); ++d) {
      //  fsens[d][0] = permute_layout(fseed[d][0], Relayout(relay_.source(), relay_.target(),"fwdseed"));
      //}
      casadi_int n = fseed.size();
      //uout() << "ad_forward prepre" << std::endl;
      std::vector<MX> seeds;
      for (casadi_int d=0; d<n; ++d) {
        seeds.push_back(fseed[d][0]);//, Relayout(relay_.source(), relay_.target(),"fwdseed"));
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
      MX sens = permute_layout(all_seeds, relay_.push_right(n));//Relayout(relay_.source().push_right(n), relay_.target().push_right(n), "fwdseed"));
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
      casadi_int n = aseed.size();      
      std::vector<MX> seeds;
      for (casadi_int d=0; d<n; ++d) {
        seeds.push_back(aseed[d][0]);
      }
      MX all_seeds = horzcat(seeds);
      Sparsity target_sp = repmat(sparsity(), 1, n);
      casadi_assert((target_sp+all_seeds.sparsity()).nnz()==target_sp.nnz(), "seed sparsity must be subset");
      if (all_seeds.sparsity()!=target_sp) all_seeds = project(all_seeds, target_sp);
      MX sens = permute_layout(all_seeds, relay_.push_right(n).invert());//Relayout(relay_.target().push_right(n), relay_.source().push_right(n), "adjseed"));
      std::vector<MX> senss = horzsplit(sens, sens.size2()/n);
      for (casadi_int d=0; d<aseed.size(); ++d) {
        asens[d][0] += senss[d];
      }
    } catch (exception& e) {
      THROW_ERROR("ad_reverse", e.what());
    }
  }

  void ReinterpretLayout::ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = reinterpret_layout(fseed[d][0], layout_);
    }
  }

  template<typename T>
  int ReinterpretLayout::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    // Perform operation
    if (arg[0]!=res[0]) {
      copy(arg[0], arg[0]+nnz(), res[0]);
    }
    return 0;
  }

  MX PermuteLayout::get_nzref(const Sparsity& sp, const vector<casadi_int>& nz) const {
    // TODO correct 
    //return dep(0)->get_nzref(sp, nz);
    return MXNode::get_nzref(sp, nz);
  }

  int ReinterpretLayout::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen(arg, res, iw, w);
  }

  int ReinterpretLayout::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen(arg, res, iw, w);
  }

  void PermuteLayout::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    uout() << "eval_mx" << std::endl;
    sparsity().spy();
    uout() << "foo" << std::endl;
    arg[0].sparsity().spy();
    res[0] = project(arg[0],sparsity())->get_permute_layout(Relayout(relay_.source(), relay_.perms(), relay_.target(), "eval_mx_"+relay_.label_));
  }

  void ReinterpretLayout::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0]->get_reinterpret_layout(layout_);
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

  MX PermuteLayout::get_permute_layout(const Relayout& relay) const {
    return MXNode::get_permute_layout(relay);
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

  void ReinterpretLayout::generate(CodeGenerator& g,
                        const std::vector<casadi_int>& arg,
                        const std::vector<casadi_int>& res) const {
    // Copy if not inplace
    if (arg[0]!=res[0]) {
      if (nnz()==1) {
        g << g.workel(res[0]) << " = " << g.workel(arg[0]) << ";\n";
      } else {
        g << g.copy(g.work(arg[0], nnz()), nnz(), g.work(res[0], nnz())) << "\n";
      }
    }
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
