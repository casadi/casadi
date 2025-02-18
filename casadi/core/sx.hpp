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

#ifndef CASADI_SX_HPP
#define CASADI_SX_HPP

#include "sx_fwd.hpp"
#include "matrix_decl.hpp"

namespace casadi {

template<> inline std::string matrixName<SXElem>() { return "SX"; }



  template<>
  bool SX::__nonzero__() const;
  template<>
  void SX::set_max_depth(casadi_int eq_depth);
  template<>
  casadi_int SX::get_max_depth();
  template<>
  SX SX::_sym(const std::string& name, const Sparsity& sp);

  template<>
  bool SX::is_regular() const;

  template<>
  bool SX::is_smooth() const;

  template<>
  casadi_int SX::element_hash() const;

  template<>
  bool SX::is_leaf() const;

  template<>
  bool SX::is_commutative() const;

  template<>
  bool SX::is_valid_input() const;

  template<>
  bool SX::is_symbolic() const;

  template<>
  bool SX::is_call() const;

  template<>
  std::vector<SXElem> SX::call(const Function& f, const std::vector<SXElem>& dep);

  template<>
  bool SX::is_output() const;

  template<>
  bool SX::has_output() const;

  template<>
  SX SX::get_output(casadi_int oind) const;

  template<>
  Function SX::which_function() const;

  template<>
  casadi_int SX::which_output() const;

  template<>
  casadi_int SX::op() const;

  template<>
  bool SX::is_op(casadi_int op) const;

  template<> bool SX::has_duplicates() const;

  template<> void SX::reset_input() const;

  template<>
  std::string SX::name() const;
  template<>
  SX SX::dep(casadi_int ch) const;

  template<>
  casadi_int SX::n_dep() const;

  template<>
  void SX::expand(const SX& ex2, SX& ww, SX& tt);

  template<>
  SX SX::pw_const(const SX& t, const SX& tval, const SX& val);

  template<>
  SX SX::pw_lin(const SX& t, const SX& tval, const SX& val);

  template<>
  SX SX::gauss_quadrature(const SX& f, const SX& x, const SX& a, const SX& b, casadi_int order,
                          const SX& w);

  template<>
  SX SX::simplify(const SX& x);

  template<>
  std::vector<SX>
  SX::substitute(const std::vector<SX>& ex, const std::vector<SX>& v, const std::vector<SX>& vdef);

  template<>
  SX SX::substitute(const SX& ex, const SX& v, const SX& vdef);

  template<>
  void SX::substitute_inplace(const std::vector<SX>& v, std::vector<SX>& vdef,
                             std::vector<SX>& ex, bool reverse);

  template<>
  void SX::extract_parametric(const SX &expr, const SX& par,
      SX& expr_ret, std::vector<SX>& symbols, std::vector<SX>& parametric,
      const Dict& opts);

  template<>
  void SX::separate_linear(const SX &expr,
    const SX &sym_lin, const SX &sym_const,
    SX& expr_const, SX& expr_lin, SX& expr_nonlin);

  template<>
  std::vector<SX> SX::cse(const std::vector<SX>& e);

  template<>
  bool SX::depends_on(const SX &x, const SX &arg);

  template<>
  bool SX::contains_all(const std::vector<SX>& v, const std::vector<SX> &n);

  template<>
  bool SX::contains_any(const std::vector<SX>& v, const std::vector<SX> &n);

  template<>
  SX SX::jacobian(const SX &f, const SX &x, const Dict& opts);
  template<>
  SX SX::hessian(const SX &ex, const SX &arg, SX &g, const Dict& opts);
  template<>
  SX SX::hessian(const SX &ex, const SX &arg, const Dict& opts);

  template<>
  std::vector<std::vector<SX> >
  SX::forward(const std::vector<SX> &ex, const std::vector<SX> &arg,
          const std::vector<std::vector<SX> > &v, const Dict& opts);

  template<>
  std::vector<std::vector<SX> >
  SX::reverse(const std::vector<SX> &ex, const std::vector<SX> &arg,
          const std::vector<std::vector<SX> > &v, const Dict& opts);

  template<>
  std::vector<bool> SX::which_depends(const SX &expr, const SX &var, casadi_int order, bool tr);

  template<>
  Sparsity SX::jacobian_sparsity(const SX &f, const SX &x);

  template<>
  SX SX::taylor(const SX& f, const SX& x,
                const SX& a, casadi_int order);

  template<>
  SX SX::mtaylor(const SX& f, const SX& x, const SX& a, casadi_int order,
                 const std::vector<casadi_int>& order_contributions);

  template<>
  SX SX::mtaylor(const SX& f, const SX& x, const SX& a, casadi_int order);

  template<>
  casadi_int SX::n_nodes(const SX& x);

  template<>
  std::string
  SX::print_operator(const SX& X, const std::vector<std::string>& args);

  template<>
  std::vector<SX> SX::symvar(const SX& x);

  template<>
  void SX::extract(std::vector<SX>& ex, std::vector<SX>& v_sx,
      std::vector<SX>& vdef_sx, const Dict& opts);

  template<>
  void SX::shared(std::vector<SX>& ex,
                         std::vector<SX>& v_sx,
                         std::vector<SX> & vdef_sx,
                         const std::string& v_prefix,
                         const std::string& v_suffix);

  template<>
  SX SX::poly_coeff(const SX& ex, const SX& x);

  template<>
  SX SX::poly_roots(const SX& p);

  template<>
  SX SX::eig_symbolic(const SX& m);

  template<>
  void SX::print_split(casadi_int nnz, const SXElem* nonzeros, std::vector<std::string>& nz,
                      std::vector<std::string>& inter);

  template<> std::vector<SX> SX::get_input(const Function& f);
  template<> std::vector<SX> SX::get_free(const Function& f);
  template<>
  Dict CASADI_EXPORT SX::info() const;

  template<>
  void CASADI_EXPORT SX::to_file(const std::string& filename, const Sparsity& sp,
    const SXElem* nonzeros, const std::string& format_hint);

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
  template<>
  std::mutex& SX::get_mutex_temp();
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

#ifndef CASADI_SX_INSTANTIATOR_CPP
  extern template class Matrix<SXElem>;
#endif // CASADI_SX_INSTANTIATOR_CPP

} // namespace casadi


#endif // CASADI_SX_HPP
