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


#ifndef CASADI_GETNONZEROS_HPP
#define CASADI_GETNONZEROS_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {
  /** \brief Get nonzeros of a matrix

      \author Joel Andersson
      \date 2013

      \identifier{i3} */
  class CASADI_EXPORT GetNonzeros : public MXNode {
  public:

    ///@{
    /// Create functions
    static MX create(const Sparsity& sp, const MX& x, const std::vector<casadi_int>& nz);
    static MX create(const Sparsity& sp, const MX& x, const Slice& s);
    static MX create(const Sparsity& sp, const MX& x, const Slice& inner, const Slice& outer);
    ///@}

    /// Constructor
    GetNonzeros(const Sparsity& sp, const MX& y);

    /// Destructor
    ~GetNonzeros() override {}

    /** \brief  Evaluate symbolically (MX)

        \identifier{i4} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{i5} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{i6} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /// Get an IM representation of a GetNonzeros or SetNonzeros node
    Matrix<casadi_int> mapping() const override;

    /// Get all the nonzeros
    virtual std::vector<casadi_int> all() const = 0;

    /** \brief Get the operation

        \identifier{i7} */
    casadi_int op() const override { return OP_GETNONZEROS;}

    /// Get the nonzeros of matrix
    MX get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const override;

    /** \brief Deserialize without type information

        \identifier{i8} */
    static MXNode* deserialize(DeserializingStream& s);

  protected:
    /** \brief Deserializing constructor

        \identifier{i9} */
    explicit GetNonzeros(DeserializingStream& s) : MXNode(s) {}
  };

  class CASADI_EXPORT GetNonzerosVector : public GetNonzeros {
  public:
    /// Constructor
    GetNonzerosVector(const Sparsity& sp, const MX& x,
                      const std::vector<casadi_int>& nz) : GetNonzeros(sp, x), nz_(nz) {}

    /// Destructor
    ~GetNonzerosVector() override {}

    /// Get all the nonzeros
    std::vector<casadi_int> all() const override { return nz_;}

    /** \brief  Propagate sparsity forward

        \identifier{ia} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{ib} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{ic} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Print expression

        \identifier{id} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{ie} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{if} */
    bool is_equal(const MXNode* node, casadi_int depth) const override;

    /** Obtain information about node */
    Dict info() const override { return {{"nz", nz_}}; }

    /// Operation sequence
    std::vector<casadi_int> nz_;

    /** \brief Serialize an object without type information

        \identifier{ig} */
    void serialize_body(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{ih} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{ii} */
    explicit GetNonzerosVector(DeserializingStream& s);
  };

  // Specialization of the above when nz_ is a Slice
  class CASADI_EXPORT GetNonzerosSlice : public GetNonzeros {
  public:

    /// Constructor
    GetNonzerosSlice(const Sparsity& sp, const MX& x, const Slice& s) : GetNonzeros(sp, x), s_(s) {}

    /// Destructor
    ~GetNonzerosSlice() override {}

    /// Get all the nonzeros
    std::vector<casadi_int> all() const override { return s_.all(s_.stop);}

    /** \brief  Propagate sparsity forward

        \identifier{ij} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{ik} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Print expression

        \identifier{il} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{im} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{in} */
    bool is_equal(const MXNode* node, casadi_int depth) const override;

    /** Obtain information about node */
    Dict info() const override { return {{"slice", s_.info()}}; }

    // Data member
    Slice s_;

    /** \brief Serialize an object without type information

        \identifier{io} */
    void serialize_body(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{ip} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{iq} */
    explicit GetNonzerosSlice(DeserializingStream& s);
  };

  // Specialization of the above when nz_ is a nested Slice
  class CASADI_EXPORT GetNonzerosSlice2 : public GetNonzeros {
  public:

    /// Constructor
    GetNonzerosSlice2(const Sparsity& sp, const MX& x, const Slice& inner,
                      const Slice& outer) : GetNonzeros(sp, x), inner_(inner), outer_(outer) {}

    /// Destructor
    ~GetNonzerosSlice2() override {}

    /// Get all the nonzeros
    std::vector<casadi_int> all() const override { return inner_.all(outer_, outer_.stop);}

    /** \brief  Propagate sparsity forward

        \identifier{ir} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{is} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Print expression

        \identifier{it} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{iu} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{iv} */
    bool is_equal(const MXNode* node, casadi_int depth) const override;

    /** Obtain information about node */
    Dict info() const override { return {{"inner", inner_.info()}, {"outer", outer_.info()}}; }

    // Data members
    Slice inner_, outer_;

    /** \brief Serialize an object without type information

        \identifier{iw} */
    void serialize_body(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{ix} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{iy} */
    explicit GetNonzerosSlice2(DeserializingStream& s);
  };


} // namespace casadi
/// \endcond

#endif // CASADI_GETNONZEROS_HPP
