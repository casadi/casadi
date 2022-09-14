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


#ifndef CASADI_LOW_HPP
#define CASADI_LOW_HPP

#include "mx_node.hpp"
/// \cond INTERNAL

namespace casadi {
  /** \brief Lows the first nonzero element in a vector
      \author Joel Andersson
      \date 2019
      \identifier{15e} */
  class CASADI_EXPORT Low : public MXNode {
  public:
    /** \brief  Constructor
        \identifier{15f} */
    Low(const MX& v, const MX& p, const Dict& opts);

    /** \brief  Destructor
        \identifier{15g} */
    ~Low() override {}

    /** \brief  Print expression
        \identifier{15h} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /** \brief  Evaluate symbolically (MX)
        \identifier{15i} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives
        \identifier{15j} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives
        \identifier{15k} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Propagate sparsity forward
        \identifier{15l} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards
        \identifier{15m} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Get the operation
        \identifier{15n} */
    casadi_int op() const override { return OP_LOW;}

    /** \brief Generate code for the operation
        \identifier{15o} */
    void generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res) const override;

    /** \brief Serialize an object without type information
        \identifier{15p} */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Deserialize without type information
        \identifier{15q} */
    static MXNode* deserialize(DeserializingStream& s) { return new Low(s); }

    static casadi_int interpret_lookup_mode(const std::string& lookup_mode, casadi_int n);
    static std::string lookup_mode_from_enum(casadi_int lookup_mode);

  protected:
    /** \brief Deserializing constructor
        \identifier{15r} */
    explicit Low(DeserializingStream& s);

  private:
    casadi_int lookup_mode_;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_LOW_HPP
