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
#include "linsol_internal.hpp"
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
    Determinant(const MX& x, const Linsol& linsol);

    /// Destructor
    ~Determinant() override {}

    /** \brief Evaluate the function numerically */
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /** \brief Evaluate the function symbolically (SX) */
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief Length of w the generated code needs (eval uses Linsol memory) */
    size_t codegen_sz_w() const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{xg} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res,
        const std::vector<bool>& unique={}) const override;

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

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Deserialize without type information

        \identifier{xl} */
    static MXNode* deserialize(DeserializingStream& s) { return new Determinant(s); }

    // Linear solver factorizing the matrix (constructed by MX::det, plugin
    // selects the algorithm); the node owns it and serializes it directly
    Linsol linsol_;

  protected:
    /** \brief Deserializing constructor

        \identifier{xm} */
    explicit Determinant(DeserializingStream& s);

  };


} // namespace casadi
/// \endcond

#endif // CASADI_DETERMINANT_HPP
