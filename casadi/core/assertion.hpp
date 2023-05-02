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


#ifndef CASADI_ASSERTION_HPP
#define CASADI_ASSERTION_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL
namespace casadi {
  /** \brief Assertion

      \author Joris Gillis
      \date 2013

      \identifier{139} */
  class CASADI_EXPORT Assertion : public MXNode {
  public:

    /// Constructor
    Assertion(const MX& x, const MX& y, const std::string & fail_message);

    /// Destructor
    ~Assertion() override {}

    /** \brief  Evaluate symbolically (MX)

        \identifier{13a} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{13b} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{13c} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Propagate sparsity forward

        \identifier{13d} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{13e} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Generate code for the operation

        \identifier{13f} */
    void generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res) const override;

    /** \brief  Print expression

        \identifier{13g} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{13h} */
    casadi_int op() const override { return OP_ASSERTION;}

    /// Can the operation be performed inplace (i.e. overwrite the result)
    casadi_int n_inplace() const override { return 1;}

    /** \brief Serialize an object without type information

        \identifier{13i} */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Deserialize without type information

        \identifier{13j} */
    static MXNode* deserialize(DeserializingStream& s) { return new Assertion(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{13k} */
    explicit Assertion(DeserializingStream& s);

  private:
    std::string fail_message_;
  };


} // namespace casadi

/// \endcond

#endif // CASADI_ASSERTION_HPP
