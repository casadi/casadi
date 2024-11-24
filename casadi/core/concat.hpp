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


#ifndef CASADI_CONCAT_HPP
#define CASADI_CONCAT_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {
  /** \brief Concatenation: Join multiple expressions stacking the nonzeros

      \author Joel Andersson
      \date 2014

      \identifier{146} */
  class CASADI_EXPORT Concat : public MXNode {
  public:

    /// Constructor
    Concat(const std::vector<MX>& x);

    /// Destructor
    ~Concat() override = 0;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Propagate sparsity forward

        \identifier{147} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{148} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Generate code for the operation

        \identifier{149} */
    void generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res,
                          const std::vector<bool>& arg_is_ref,
                          std::vector<bool>& res_is_ref) const override;

    /// Get the nonzeros of matrix
    MX get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const override;

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{14a} */
    bool is_equal(const MXNode* node, casadi_int depth) const override {
      return sameOpAndDeps(node, depth);
    }

    /** \brief  Check if valid function input

        \identifier{14b} */
    bool is_valid_input() const override;

    /** \brief Get the number of symbolic primitives

        \identifier{14c} */
    casadi_int n_primitives() const override;

    /** \brief Get symbolic primitives

        \identifier{14d} */
    void primitives(std::vector<MX>::iterator& it) const override;

    /** \brief Detect duplicate symbolic expressions

        \identifier{14e} */
    bool has_duplicates() const override;

    /** \brief Reset the marker for an input expression

        \identifier{14f} */
    void reset_input() const override;

  protected:
    /** \brief Deserializing constructor

        \identifier{14g} */
    explicit Concat(DeserializingStream& s) : MXNode(s) {}
  };


  /** \brief Horizontal concatenation

      \author Joel Andersson
      \date 2013

      \identifier{14h} */
  class CASADI_EXPORT Horzcat : public Concat {
  public:

    /// Constructor
    Horzcat(const std::vector<MX>& x);

    /// Destructor
    ~Horzcat() override {}

    /** \brief  Print expression

        \identifier{14i} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{14j} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{14k} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{14l} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Get the operation

        \identifier{14m} */
    casadi_int op() const override { return OP_HORZCAT;}

    /// Split up an expression along primitives (template)
    template<typename T>
    void split_primitives_gen(const T& x, typename std::vector<T>::iterator& it) const;

    /// @{
    /** \brief Split up an expression along symbolic primitives

        \identifier{14n} */
    void split_primitives(const MX& x, std::vector<MX>::iterator& it) const override;
    void split_primitives(const SX& x, std::vector<SX>::iterator& it) const override;
    void split_primitives(const DM& x, std::vector<DM>::iterator& it) const override;
    /// @}

    /// Join an expression along symbolic primitives (template)
    template<typename T>
    T join_primitives_gen(typename std::vector<T>::const_iterator& it) const;

    /// @{
    /** \brief Join an expression along symbolic primitives

        \identifier{14o} */
    MX join_primitives(std::vector<MX>::const_iterator& it) const override;
    SX join_primitives(std::vector<SX>::const_iterator& it) const override;
    DM join_primitives(std::vector<DM>::const_iterator& it) const override;
    /// @}

    /** \brief Get offsets for split

        \identifier{14p} */
    std::vector<casadi_int> off() const;

    /** \brief Deserialize without type information

        \identifier{14q} */
    static MXNode* deserialize(DeserializingStream& s) { return new Horzcat(s); }
  protected:
    /** \brief Deserializing constructor

        \identifier{14r} */
    explicit Horzcat(DeserializingStream& s) : Concat(s) {}
  };

  /** \brief Vertical concatenation of vectors

      \author Joel Andersson
      \date 2014

      \identifier{14s} */
  class CASADI_EXPORT Vertcat : public Concat {
  public:

    /// Constructor
    Vertcat(const std::vector<MX>& x);

    /// Destructor
    ~Vertcat() override {}

    /** \brief  Print expression

        \identifier{14t} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{14u} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{14v} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{14w} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Get the operation

        \identifier{14x} */
    casadi_int op() const override { return OP_VERTCAT;}

    /// Split up an expression along primitives (template)
    template<typename T>
    void split_primitives_gen(const T& x, typename std::vector<T>::iterator& it) const;

    /// @{
    /** \brief Split up an expression along symbolic primitives

        \identifier{14y} */
    void split_primitives(const MX& x, std::vector<MX>::iterator& it) const override;
    void split_primitives(const SX& x, std::vector<SX>::iterator& it) const override;
    void split_primitives(const DM& x, std::vector<DM>::iterator& it) const override;
    /// @}

    /// Join an expression along symbolic primitives (template)
    template<typename T>
    T join_primitives_gen(typename std::vector<T>::const_iterator& it) const;

    /// @{
    /** \brief Join an expression along symbolic primitives

        \identifier{14z} */
    MX join_primitives(std::vector<MX>::const_iterator& it) const override;
    SX join_primitives(std::vector<SX>::const_iterator& it) const override;
    DM join_primitives(std::vector<DM>::const_iterator& it) const override;
    /// @}

    /** \brief Get offsets for split

        \identifier{150} */
    std::vector<casadi_int> off() const;

    /** \brief Deserialize without type information

        \identifier{151} */
    static MXNode* deserialize(DeserializingStream& s) { return new Vertcat(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{152} */
    explicit Vertcat(DeserializingStream& s) : Concat(s) {}
  };

  /** \brief Diagonal concatenation of matrices

      \author Joris Gillis
      \date 2014

      \identifier{153} */
  class CASADI_EXPORT Diagcat : public Concat {
  public:

    /// Constructor
    Diagcat(const std::vector<MX>& x);

    /// Destructor
    ~Diagcat() override {}

    /** \brief  Print expression

        \identifier{154} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{155} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{156} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{157} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Get the operation

        \identifier{158} */
    casadi_int op() const override { return OP_DIAGCAT;}

    /// Split up an expression along primitives (template)
    template<typename T>
    void split_primitives_gen(const T& x, typename std::vector<T>::iterator& it) const;

    /// @{
    /** \brief Split up an expression along symbolic primitives

        \identifier{159} */
    void split_primitives(const MX& x, std::vector<MX>::iterator& it) const override;
    void split_primitives(const SX& x, std::vector<SX>::iterator& it) const override;
    void split_primitives(const DM& x, std::vector<DM>::iterator& it) const override;
    /// @}

    /// Join an expression along symbolic primitives (template)
    template<typename T>
    T join_primitives_gen(typename std::vector<T>::const_iterator& it) const;

    /// @{
    /** \brief Join an expression along symbolic primitives

        \identifier{15a} */
    MX join_primitives(std::vector<MX>::const_iterator& it) const override;
    SX join_primitives(std::vector<SX>::const_iterator& it) const override;
    DM join_primitives(std::vector<DM>::const_iterator& it) const override;
    /// @}

    /** \brief Get offsets for split

        \identifier{15b} */
    std::pair<std::vector<casadi_int>, std::vector<casadi_int> > off() const;

    /** \brief Deserialize without type information

        \identifier{15c} */
    static MXNode* deserialize(DeserializingStream& s) { return new Diagcat(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{15d} */
    explicit Diagcat(DeserializingStream& s) : Concat(s) {}
  };

} // namespace casadi
/// \endcond

#endif // CASADI_CONCAT_HPP
