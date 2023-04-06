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


#ifndef CASADI_DETERMINANT_HPP
#define CASADI_DETERMINANT_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {
  /** \brief Matrix determinant

      \author Joel Andersson
      \date 2013

      \identifier{xf} */
  class CASADI_EXPORT Determinant : public MXNode {
  public:

    /// Constructor
    Determinant(const MX& x);

    /// Destructor
    ~Determinant() override {}

    /** \brief  Evaluate symbolically (MX)

        \identifier{xg} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{xh} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{xi} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Print expression

        \identifier{xj} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{xk} */
    casadi_int op() const override { return OP_DETERMINANT;}

    /** \brief Deserialize without type information

        \identifier{xl} */
    static MXNode* deserialize(DeserializingStream& s) { return new Determinant(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{xm} */
    explicit Determinant(DeserializingStream& s) : MXNode(s) {}

  };


} // namespace casadi
/// \endcond

#endif // CASADI_DETERMINANT_HPP
