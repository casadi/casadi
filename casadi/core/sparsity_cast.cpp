/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#include "sparsity_cast.hpp"
#include "casadi_misc.hpp"
#include "getnonzeros.hpp"

namespace casadi {

  SparsityCast::SparsityCast(const MX& x, Sparsity sp) {
    casadi_assert_dev(x.nnz()==sp.nnz());
    set_dep(x);
    set_sparsity(sp);
  }

  int SparsityCast::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int SparsityCast::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  int SparsityCast::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    if (arg[0]!=res[0]) std::copy(arg[0], arg[0]+nnz(), res[0]);
    return 0;
  }

  int SparsityCast::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    copy_fwd(arg[0], res[0], nnz());
    return 0;
  }

  int SparsityCast::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    copy_rev(arg[0], res[0], nnz());
    return 0;
  }

  std::string SparsityCast::disp(const std::vector<std::string>& arg) const {
    if (sparsity().is_dense() && sparsity().is_column()) {
      return "nonzeros(" + arg.at(0) + ")";
    } else {
      return "sparsity_cast(" + arg.at(0) + ")";
    }
  }

  void SparsityCast::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    casadi_assert(arg[0].sparsity().is_subset(dep().sparsity()), "SparsityCast: Invalid sparsity");
    Sparsity sp = arg[0].sparsity().sparsity_cast_mod(dep().sparsity(), sparsity());
    res[0] = sparsity_cast(arg[0], sp);
  }

  void SparsityCast::ad_forward(const std::vector<std::vector<MX> >& fseed,
                        std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d = 0; d<fsens.size(); ++d) {
      casadi_assert_dev(fseed[d][0].sparsity().is_subset(dep().sparsity()));
      Sparsity sp = fseed[d][0].sparsity().sparsity_cast_mod(dep().sparsity(), sparsity());
      fsens[d][0] = sparsity_cast(fseed[d][0], sp);
    }
  }

  void SparsityCast::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                        std::vector<std::vector<MX> >& asens) const {
    for (casadi_int d=0; d<aseed.size(); ++d) {
      MX seed = aseed[d][0];
      if (!seed.sparsity().is_subset(sparsity())) seed = seed(sparsity());
      Sparsity sp = seed.sparsity().sparsity_cast_mod(sparsity(), dep().sparsity());
      asens[d][0] += sparsity_cast(seed, sp);
    }
  }

  void SparsityCast::generate(CodeGenerator& g,
                         const std::vector<casadi_int>& arg,
                         const std::vector<casadi_int>& res,
                         const std::vector<bool>& arg_is_ref,
                         std::vector<bool>& res_is_ref) const {
    generate_copy(g, arg, res, arg_is_ref, res_is_ref, 0);
  }

  MX SparsityCast::get_reshape(const Sparsity& sp) const {
    if (sp.is_reshape(dep(0).sparsity())) {
      return reshape(dep(0), sp);
    } else {
      return MXNode::get_reshape(sp);
    }
  }

  MX SparsityCast::get_sparsity_cast(const Sparsity& sp) const {
    return sparsity_cast(dep(0), sp);
  }

  MX SparsityCast::get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const {
    return GetNonzeros::create(sp, dep(), nz);
  }

  MX SparsityCast::get_transpose() const {
    // For vectors, reshape is also a transpose
    if (sparsity().is_vector()) {
      return dep()->get_sparsity_cast(sparsity().T());
    } else {
      return MXNode::get_transpose();
    }
  }

  bool SparsityCast::is_valid_input() const {
    return dep()->is_valid_input();
  }

  casadi_int SparsityCast::n_primitives() const {
    return dep()->n_primitives();
  }

  void SparsityCast::primitives(std::vector<MX>::iterator& it) const {
    dep()->primitives(it);
  }

  template<typename T>
  void SparsityCast::split_primitives_gen(const T& x,
      typename std::vector<T>::iterator& it) const {
    dep()->split_primitives(sparsity_cast(project(x, sparsity()), dep().sparsity()), it);
  }

  void SparsityCast::split_primitives(const MX& x, std::vector<MX>::iterator& it) const {
    split_primitives_gen<MX>(x, it);
  }

  void SparsityCast::split_primitives(const SX& x, std::vector<SX>::iterator& it) const {
    split_primitives_gen<SX>(x, it);
  }

  void SparsityCast::split_primitives(const DM& x, std::vector<DM>::iterator& it) const {
    split_primitives_gen<DM>(x, it);
  }

  template<typename T>
  T SparsityCast::join_primitives_gen(typename std::vector<T>::const_iterator& it) const {
    return sparsity_cast(project(dep()->join_primitives(it), dep().sparsity()), sparsity());
  }

  MX SparsityCast::join_primitives(std::vector<MX>::const_iterator& it) const {
    return join_primitives_gen<MX>(it);
  }

  SX SparsityCast::join_primitives(std::vector<SX>::const_iterator& it) const {
    return join_primitives_gen<SX>(it);
  }

  DM SparsityCast::join_primitives(std::vector<DM>::const_iterator& it) const {
    return join_primitives_gen<DM>(it);
  }

  bool SparsityCast::has_duplicates() const {
    return dep()->has_duplicates();
  }

  void SparsityCast::reset_input() const {
    dep()->reset_input();
  }

} // namespace casadi
