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


#ifndef CASADI_NLP_TOOLS_HPP
#define CASADI_NLP_TOOLS_HPP

#include "casadi/core/function.hpp"

namespace casadi {

  //@{
  /** \brief Detect simple bounds from general constraints
   *
   * Given parametric constraints:
   * \verbatim
   *   subject to lbg(p) <= g(x,p) <= ubg(p)
   * \endverbatim
   *
   * Returns an equivalent set
   * \verbatim
   *   subject to  lbg(p)(gi) <= g(x,p)(gi) <= ubg(p)(gi)
   *               lbx(p) <= x                 <= ubx(p)
   * \endverbatim
   *
   * \param[out] lam_forward (lam_g,p)->(lam_sg,lam_x)
   * \param[out] lam_backward (lam_sg,lam_x,p)->(lam_g)
   *
      \identifier{1sw} */
  CASADI_EXPORT void detect_simple_bounds(const SX& xX, const SX& p,
      const SX& g, const SX& lbg, const SX& ubg,
      std::vector<casadi_int>& SWIG_OUTPUT(gi),
      SX& SWIG_OUTPUT(lbx), SX& SWIG_OUTPUT(ubx),
      Function& SWIG_OUTPUT(lam_forward),
      Function& SWIG_OUTPUT(lam_backward));
  CASADI_EXPORT void detect_simple_bounds(const MX& xX, const MX& p,
      const MX& g, const MX& lbg, const MX& ubg,
      std::vector<casadi_int>& SWIG_OUTPUT(gi),
      MX& SWIG_OUTPUT(lbx), MX& SWIG_OUTPUT(ubx),
      Function& SWIG_OUTPUT(lam_forward),
      Function& SWIG_OUTPUT(lam_backward));
  //@}

/*
  CASADI_EXPORT void detect_simple_bounds(const SX& xX,
      const SX& g, const SX& lbg, const SX& ubg,
      std::vector<casadi_int>& SWIG_OUTPUT(gi),
      DM& SWIG_OUTPUT(lbx), DM& SWIG_OUTPUT(ubx),
      Function& SWIG_OUTPUT(lam_forward),
      Function& SWIG_OUTPUT(lam_backward));
*/


  /** \brief Check sos structure and generate defaults

      \identifier{1sx} */
  template <class T>
  void check_sos(casadi_int nx, const std::vector< std::vector<T> >& groups,
                  std::vector< std::vector<double> >& weights,
                  std::vector< casadi_int >& types) {
    // Checks on presence
    if (groups.empty()) {
      casadi_assert(weights.empty(), "Missing sos_groups.");
      casadi_assert(types.empty(),   "Missing sos_groups.");
    }

    // Checks on dimensions
    casadi_int sos_num = groups.size();

    casadi_assert(weights.empty() || weights.size()==sos_num,
      "sos_weights has incorrect size");

    // Set default types
    if (!groups.empty() && types.empty())
      types.resize(sos_num, 1);

    // Set default weights
    if (weights.empty()) {
      for (const auto& e : groups) {
        std::vector<double> w(e.size());
        for (casadi_int i=0;i<w.size();++i) w[i] = i;
        weights.push_back(w);
      }
    }

    casadi_assert(types.size()==sos_num,
      "sos_types has incorrect size");

    // Group-wise dimension check
    for (casadi_int i=0;i<weights.size();++i) {
      casadi_assert(weights[i].size()==groups[i].size(),
        "Dimension mismatch in weights for group " + str(i) + ": "
        "Expected " + str(groups[i].size()) + ", got " + str(weights[i].size()));
    }

    // Checks on contents
    for (casadi_int t : types) casadi_assert(t==1 || t==2, "SOS type must be either 1 or 2.");
    for (const auto& v : groups) casadi_assert(in_range(v, 0, nx), "Index out of bound");
  }

} // namespace casadi

#endif // CASADI_NLP_TOOLS_HPP
