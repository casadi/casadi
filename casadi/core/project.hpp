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


#ifndef CASADI_PROJECT_HPP
#define CASADI_PROJECT_HPP

#include "mx_node.hpp"
/// \cond INTERNAL

namespace casadi {
  /** \brief Change the sparsity of an expression

      \author Joel Andersson
      \date 2011-2013

      \identifier{1ih} */
  class CASADI_EXPORT Project : public MXNode {
  public:

    /** \brief  Constructor

        \identifier{1ii} */
    Project(const MX& x, const Sparsity& sp);

    /** \brief  Destructor

        \identifier{1ij} */
    ~Project() override {}

    /** \brief  Print expression

        \identifier{1ik} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{1il} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{1im} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{1in} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Generate code for the operation

        \identifier{1io} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    /** \brief  Propagate sparsity forward

        \identifier{1ip} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{1iq} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Get the operation

        \identifier{1ir} */
    casadi_int op() const override { return OP_PROJECT;}

    /** \brief Get required length of w field

        \identifier{1is} */
    size_t sz_w() const override { return size1();}

    /** \brief Serialize type information

        \identifier{1it} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserialize without type information

        \identifier{1iu} */
    static MXNode* deserialize(DeserializingStream& s);

  protected:
    /** \brief Deserializing constructor

        \identifier{1iv} */
    explicit Project(DeserializingStream& s) : MXNode(s) {}

  };


  /** \brief Densify

      \author Joris Gillis
      \date 2019

      \identifier{1iw} */
  class CASADI_EXPORT Densify : public Project {
  public:

    /// Constructor
    Densify(const MX& x, const Sparsity& sp) : Project(x, sp) {}

    /// Destructor
    ~Densify() override {}

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief Generate code for the operation

        \identifier{1ix} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    /** \brief Get required length of iw field

        \identifier{1iy} */
    size_t sz_w() const override { return 0;}

    /** \brief Serialize type information

        \identifier{1iz} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{1j0} */
    explicit Densify(DeserializingStream& s) : Project(s) {}
  };

  /** \brief Densify

      \author Joris Gillis
      \date 2019

      \identifier{1j1} */
  class CASADI_EXPORT Sparsify : public Project {
  public:

    /// Constructor
    Sparsify(const MX& x, const Sparsity& sp) : Project(x, sp) {}

    /// Destructor
    ~Sparsify() override {}

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief Generate code for the operation

        \identifier{1j2} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    /** \brief Get required length of iw field

        \identifier{1j3} */
    size_t sz_w() const override { return 0;}

    /** \brief Serialize type information

        \identifier{1j4} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{1j5} */
    explicit Sparsify(DeserializingStream& s) : Project(s) {}
  };

} // namespace casadi

/// \endcond

#endif // CASADI_PROJECT_HPP
