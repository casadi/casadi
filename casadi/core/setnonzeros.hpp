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


#ifndef CASADI_SETNONZEROS_HPP
#define CASADI_SETNONZEROS_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {

  /** \brief Assign or add entries to a matrix

      \author Joel Andersson
      \date 2013

      \identifier{vx} */
  template<bool Add>
  class CASADI_EXPORT SetNonzeros : public MXNode {
  public:
    ///@{
    /** \brief Create functions

    *   returns y with y[nz]+=x

        \identifier{vy} */
    static MX create(const MX& y, const MX& x, const std::vector<casadi_int>& nz);
    static MX create(const MX& y, const MX& x, const Slice& s);
    static MX create(const MX& y, const MX& x, const Slice& inner, const Slice& outer);
    ///@}

    /// Constructor
    SetNonzeros(const MX& y, const MX& x);

    /// Destructor
    ~SetNonzeros() override = 0;

    /// Get all the nonzeros
    virtual std::vector<casadi_int> all() const = 0;

    /** \brief  Evaluate symbolically (MX)

        \identifier{vz} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{w0} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{w1} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Get the operation

        \identifier{w2} */
    casadi_int op() const override { return Add ? OP_ADDNONZEROS : OP_SETNONZEROS;}

    /// Get an IM representation of a GetNonzeros or SetNonzeros node
    Matrix<casadi_int> mapping() const override;

    /// Can the operation be performed inplace (i.e. overwrite the result)
    casadi_int n_inplace() const override { return 1;}

    /** \brief Deserialize with type disambiguation

        \identifier{w3} */
    static MXNode* deserialize(DeserializingStream& s);

  protected:
    /** \brief Deserializing constructor

        \identifier{w4} */
    explicit SetNonzeros(DeserializingStream& s) : MXNode(s) {}
  };


    /** \brief Add the nonzeros of a matrix to another matrix

      \author Joel Andersson
      \date 2013

        \identifier{w5} */
  template<bool Add>
  class CASADI_EXPORT SetNonzerosVector : public SetNonzeros<Add>{
  public:

    /// Constructor
    SetNonzerosVector(const MX& y, const MX& x, const std::vector<casadi_int>& nz);

    /// Destructor
    ~SetNonzerosVector() override {}

    /// Get all the nonzeros
    std::vector<casadi_int> all() const override { return nz_;}

    /** \brief  Evaluate symbolically (MX)

        \identifier{w6} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Propagate sparsity forward

        \identifier{w7} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{w8} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Print expression

        \identifier{w9} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{wa} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{wb} */
    bool is_equal(const MXNode* node, casadi_int depth) const override;

    /** Obtain information about node */
    Dict info() const override { return {{"nz", nz_}, {"add", Add}}; }

    /// Operation sequence
    std::vector<casadi_int> nz_;

    /** \brief Serialize an object without type information

        \identifier{wc} */
    void serialize_body(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{wd} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{we} */
    explicit SetNonzerosVector(DeserializingStream& s);
  };

  // Specialization of the above when nz_ is a Slice
  template<bool Add>
  class CASADI_EXPORT SetNonzerosSlice : public SetNonzeros<Add>{
  public:

    /// Constructor
    SetNonzerosSlice(const MX& y, const MX& x, const Slice& s) : SetNonzeros<Add>(y, x), s_(s) {}

    /// Destructor
    ~SetNonzerosSlice() override {}

    /// Get all the nonzeros
    std::vector<casadi_int> all() const override { return s_.all(s_.stop);}

    /** \brief  Propagate sparsity forward

        \identifier{wf} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{wg} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{wh} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Print expression

        \identifier{wi} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{wj} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{wk} */
    bool is_equal(const MXNode* node, casadi_int depth) const override;

    /** Obtain information about node */
    Dict info() const override { return {{"slice", s_.info()}, {"add", Add}}; }

    // Data member
    Slice s_;

    /** \brief Serialize an object without type information

        \identifier{wl} */
    void serialize_body(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{wm} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{wn} */
    explicit SetNonzerosSlice(DeserializingStream& s);
  };

  // Specialization of the above when nz_ is a nested Slice
  template<bool Add>
  class CASADI_EXPORT SetNonzerosSlice2 : public SetNonzeros<Add>{
  public:

    /// Constructor
    SetNonzerosSlice2(const MX& y, const MX& x, const Slice& inner, const Slice& outer) :
        SetNonzeros<Add>(y, x), inner_(inner), outer_(outer) {}

    /// Destructor
    ~SetNonzerosSlice2() override {}

    /// Get all the nonzeros
    std::vector<casadi_int> all() const override { return inner_.all(outer_, outer_.stop);}

    /** \brief  Propagate sparsity forward

        \identifier{wo} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{wp} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{wq} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Print expression

        \identifier{wr} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{ws} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{wt} */
    bool is_equal(const MXNode* node, casadi_int depth) const override;

    /** Obtain information about node */
    Dict info() const override { return {{"inner", inner_.info()}, {"outer", outer_.info()},
                                        {"add", Add}}; }

    // Data members
    Slice inner_, outer_;

    /** \brief Serialize an object without type information

        \identifier{wu} */
    void serialize_body(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{wv} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{ww} */
    explicit SetNonzerosSlice2(DeserializingStream& s);
  };

} // namespace casadi
/// \endcond

#endif // CASADI_SETNONZEROS_HPP
