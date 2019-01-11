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
   * */
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

} // namespace casadi

#endif // CASADI_NLP_TOOLS_HPP
