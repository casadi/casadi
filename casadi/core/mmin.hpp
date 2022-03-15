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


#ifndef CASADI_MMIN_HPP
#define CASADI_MMIN_HPP

#include "mx_node.hpp"

/// \cond INTERNAL
namespace casadi {

  /** \brief Matrix minimum

      \author Joris Gillis
      \date 2017
  */
  class CASADI_EXPORT MMin : public MXNode {
  public:

    /** \brief  Constructor */
    explicit MMin(const MX& x);

    /** \brief  Destructor */
    ~MMin() override {}

    /** \brief  Print expression */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation */
    Operation op() const override { return Operation::OP_MMIN;}

    /** \brief  Evaluate numerically */
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /** \brief  Evaluate symbolically (SX) */
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g,
                      const std::vector<casadi_int>& arg,
                      const std::vector<casadi_int>& res) const override;
    /** \brief Calculate forward mode directional derivatives */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Deserialize without type information */
    static MXNode* deserialize(DeserializingStream& s) { return new MMin(s); }
  protected:
    /** \brief Deserializing constructor */
    explicit MMin(DeserializingStream& s) : MXNode(s) {}
  };

  /** \brief Matrix maximum

      \author Joris Gillis
      \date 2017
  */
  class CASADI_EXPORT MMax : public MXNode {
  public:

    /** \brief  Constructor */
    explicit MMax(const MX& x);

    /** \brief  Destructor */
    ~MMax() override {}

    /** \brief  Print expression */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation */
    Operation op() const override { return Operation::OP_MMAX;}

    /** \brief  Evaluate numerically */
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /** \brief  Evaluate symbolically (SX) */
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g,
                      const std::vector<casadi_int>& arg,
                      const std::vector<casadi_int>& res) const override;
    /** \brief Calculate forward mode directional derivatives */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Deserialize without type information */
    static MXNode* deserialize(DeserializingStream& s) { return new MMax(s); }
  protected:
    /** \brief Deserializing constructor */
    explicit MMax(DeserializingStream& s) : MXNode(s) {}
  };

} // namespace casadi

/// \endcond

#endif // CASADI_MMIN_HPP
