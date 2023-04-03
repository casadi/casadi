/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
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


#ifndef CASADI_INVERSE_HPP
#define CASADI_INVERSE_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {
  /** \brief Matrix inverse

      \author Joel Andersson
      \date 2013

      \identifier{1q3} */
  class CASADI_EXPORT Inverse : public MXNode {
  public:

    /// Constructor
    Inverse(const MX& x);

    /// Destructor
    ~Inverse() override {}

    /** \brief  Evaluate symbolically (MX)

        \identifier{1q4} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{1q5} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{1q6} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Print expression

        \identifier{1q7} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{1q8} */
    casadi_int op() const override { return OP_INVERSE;}

    /** \brief Deserialize without type information

        \identifier{1q9} */
    static MXNode* deserialize(DeserializingStream& s) { return new Inverse(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{1qa} */
    explicit Inverse(DeserializingStream& s) : MXNode(s) {}
  };


} // namespace casadi
/// \endcond

#endif // CASADI_INVERSE_HPP
