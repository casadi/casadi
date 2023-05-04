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


#ifndef CASADI_REPMAT_HPP
#define CASADI_REPMAT_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {

  /** \brief Horizontal repmat

      \author Joris Gillis
      \date 2015

      \identifier{20k} */
  class CASADI_EXPORT HorzRepmat : public MXNode {
  public:

    /// Constructor
    HorzRepmat(const MX& x, casadi_int n);

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Destructor
    ~HorzRepmat() override {}

    /** \brief  Print expression

        \identifier{20l} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{20m} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief  Propagate sparsity forward

        \identifier{20n} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{20o} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{20p} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{20q} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Generate code for the operation

        \identifier{20r} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Get the operation

        \identifier{20s} */
    casadi_int op() const override { return OP_HORZREPMAT;}

    casadi_int n_;

    /** \brief Serialize an object without type information

        \identifier{20t} */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Deserialize without type information

        \identifier{20u} */
    static MXNode* deserialize(DeserializingStream& s) { return new HorzRepmat(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{20v} */
    explicit HorzRepmat(DeserializingStream& s);
  };

  /** \brief Horizontal repsum

      \author Joris Gillis
      \date 2015

      \identifier{20w} */
  class CASADI_EXPORT HorzRepsum : public MXNode {
  public:

    /// Constructor
    HorzRepsum(const MX& x, casadi_int n);

    /// Evaluate the function (template)
    template<typename T, typename R>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w, R reduction) const;

    /// Destructor
    ~HorzRepsum() override {}

    /** \brief  Print expression

        \identifier{20x} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{20y} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief  Propagate sparsity forward

        \identifier{20z} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{210} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{211} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{212} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Generate code for the operation

        \identifier{213} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Get the operation

        \identifier{214} */
    casadi_int op() const override { return OP_HORZREPSUM;}

    casadi_int n_;

    /** \brief Serialize an object without type information

        \identifier{215} */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Deserialize without type information

        \identifier{216} */
    static MXNode* deserialize(DeserializingStream& s) { return new HorzRepsum(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{217} */
    explicit HorzRepsum(DeserializingStream& s);
  };

} // namespace casadi
/// \endcond

#endif // CASADI_REPMAT_HPP
