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


#ifndef CASADI_TRANSPOSE_HPP
#define CASADI_TRANSPOSE_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {
  /** \brief Matrix transpose
      \author Joel Andersson
      \date 2013
      \identifier{13l} */
  class CASADI_EXPORT Transpose : public MXNode {
  public:

    /// Constructor
    Transpose(const MX& x);

    /// Destructor
    ~Transpose() override {}

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Evaluate symbolically (MX)
        \identifier{13m} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives
        \identifier{13n} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives
        \identifier{13o} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Propagate sparsity forward
        \identifier{13p} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards
        \identifier{13q} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Print expression
        \identifier{13r} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation
        \identifier{13s} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Get the operation
        \identifier{13t} */
    casadi_int op() const override { return OP_TRANSPOSE;}

    /** \brief Get required length of iw field
        \identifier{13u} */
    size_t sz_iw() const override { return size2()+1;}

    /// Transpose
    MX get_transpose() const override { return dep();}

    /// Solve for square linear system
    //virtual MX get_solve(const MX& r, bool tr, const Linsol& linear_solver) const {
    // return dep()->get_solve(r, !tr, linear_solver);} // FIXME #1001

    /// Solve a system of linear equations, upper triangular A
    MX get_solve_triu(const MX& r, bool tr) const override {
      return dep()->get_solve_tril(r, !tr);
    }

    /// Solve a system of linear equations, lower triangular A
    MX get_solve_tril(const MX& r, bool tr) const override {
      return dep()->get_solve_triu(r, !tr);
    }

    /// Solve a system of linear equations, upper triangular A, unity diagonal
    MX get_solve_triu_unity(const MX& r, bool tr) const override {
      return dep()->get_solve_tril_unity(r, !tr);
    }

    /// Solve a system of linear equations, lower triangular A, unity diagonal
    MX get_solve_tril_unity(const MX& r, bool tr) const override {
      return dep()->get_solve_triu_unity(r, !tr);
    }

    /** \brief Check if two nodes are equivalent up to a given depth
        \identifier{13v} */
    bool is_equal(const MXNode* node, casadi_int depth) const override {
      return sameOpAndDeps(node, depth);
    }

    /** \brief Serialize type information
        \identifier{13w} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserialize with type disambiguation
        \identifier{13x} */
    static MXNode* deserialize(DeserializingStream& s);

  protected:
    /** \brief Deserializing constructor
        \identifier{13y} */
    explicit Transpose(DeserializingStream& s) : MXNode(s) {}
  };

  /** \brief Matrix transpose (dense)
      \author Joel Andersson
      \date 2013
      \identifier{13z} */
  class CASADI_EXPORT DenseTranspose : public Transpose {
  public:

    /// Constructor
    DenseTranspose(const MX& x) : Transpose(x) {}

    /// Destructor
    ~DenseTranspose() override {}

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Propagate sparsity forward
        \identifier{140} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards
        \identifier{141} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Generate code for the operation
        \identifier{142} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Get required length of iw field
        \identifier{143} */
    size_t sz_iw() const override { return 0;}

    /** \brief Serialize type information
        \identifier{144} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor
        \identifier{145} */
    explicit DenseTranspose(DeserializingStream& s) : Transpose(s) {}
  };



} // namespace casadi

/// \endcond

#endif // CASADI_TRANSPOSE_HPP
