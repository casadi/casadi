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
  */
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

    /** \brief  Propagate sparsity forward */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res) const override;

    /// Get the nonzeros of matrix
    MX get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const override;

    /** \brief Check if two nodes are equivalent up to a given depth */
    bool is_equal(const MXNode* node, casadi_int depth) const override {
      return sameOpAndDeps(node, depth);
    }

    /** \brief  Check if valid function input */
    bool is_valid_input() const override;

    /** \brief Get the number of symbolic primitives */
    casadi_int n_primitives() const override;

    /** \brief Get symbolic primitives */
    void primitives(std::vector<MX>::iterator& it) const override;

    /** \brief Detect duplicate symbolic expressions */
    bool has_duplicates() const override;

    /** \brief Reset the marker for an input expression */
    void reset_input() const override;

  protected:
    /** \brief Deserializing constructor */
    explicit Concat(DeserializingStream& s) : MXNode(s) {}
  };


  /** \brief Horizontal concatenation
      \author Joel Andersson
      \date 2013
  */
  class CASADI_EXPORT Horzcat : public Concat {
  public:

    /// Constructor
    Horzcat(const std::vector<MX>& x);

    /// Destructor
    ~Horzcat() override {}

    /** \brief  Print expression */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Get the operation */
    Operation op() const override { return Operation::OP_HORZCAT;}

    /** \brief Split up an expression along symbolic primitives */
    void split_primitives(const MX& x, std::vector<MX>::iterator& it) const override;

    /** \brief Join an expression along symbolic primitives */
    MX join_primitives(std::vector<MX>::const_iterator& it) const override;

    /** \brief Get offsets for split */
    std::vector<casadi_int> off() const;

    /** \brief Deserialize without type information */
    static MXNode* deserialize(DeserializingStream& s) { return new Horzcat(s); }
  protected:
    /** \brief Deserializing constructor */
    explicit Horzcat(DeserializingStream& s) : Concat(s) {}
  };

  /** \brief Vertical concatenation of vectors
      \author Joel Andersson
      \date 2014
  */
  class CASADI_EXPORT Vertcat : public Concat {
  public:

    /// Constructor
    Vertcat(const std::vector<MX>& x);

    /// Destructor
    ~Vertcat() override {}

    /** \brief  Print expression */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Get the operation */
    Operation op() const override { return Operation::OP_VERTCAT;}

    /** \brief Split up an expression along symbolic primitives */
    void split_primitives(const MX& x, std::vector<MX>::iterator& it) const override;

    /** \brief Join an expression along symbolic primitives */
    MX join_primitives(std::vector<MX>::const_iterator& it) const override;

    /** \brief Get offsets for split */
    std::vector<casadi_int> off() const;

    /** \brief Deserialize without type information */
    static MXNode* deserialize(DeserializingStream& s) { return new Vertcat(s); }

  protected:
    /** \brief Deserializing constructor */
    explicit Vertcat(DeserializingStream& s) : Concat(s) {}
  };

  /** \brief Diagonal concatenation of matrices
      \author Joris Gillis
      \date 2014
  */
  class CASADI_EXPORT Diagcat : public Concat {
  public:

    /// Constructor
    Diagcat(const std::vector<MX>& x);

    /// Destructor
    ~Diagcat() override {}

    /** \brief  Print expression */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Get the operation */
    Operation op() const override { return Operation::OP_DIAGCAT;}

    /** \brief Split up an expression along symbolic primitives */
    void split_primitives(const MX& x, std::vector<MX>::iterator& it) const override;

    /** \brief Join an expression along symbolic primitives */
    MX join_primitives(std::vector<MX>::const_iterator& it) const override;

    /** \brief Get offsets for split */
    std::pair<std::vector<casadi_int>, std::vector<casadi_int> > off() const;

    /** \brief Deserialize without type information */
    static MXNode* deserialize(DeserializingStream& s) { return new Diagcat(s); }

  protected:
    /** \brief Deserializing constructor */
    explicit Diagcat(DeserializingStream& s) : Concat(s) {}
  };

} // namespace casadi
/// \endcond

#endif // CASADI_CONCAT_HPP
