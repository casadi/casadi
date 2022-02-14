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


#ifndef CASADI_SOLVE_HPP
#define CASADI_SOLVE_HPP

#include "mx_node.hpp"
#include "casadi_call.hpp"

namespace casadi {
  /** \brief An MX atomic for linear solver solution: x = r * A^-1 or x = r * A^-T

      Forward derivatives:
      x_dot = (r_dot - x * A_dot) * A^-1

      Adjoint derivatives:
      r_bar = x_bar * A^-T
      A_bar = -x^T * r_bar

      \author Joel Andersson
      \date 2013
  */
  template<bool Tr>
  class CASADI_EXPORT Solve : public MXNode {
  public:
    /** \brief  Constructor */
    Solve(const MX& r, const MX& A);

    /** \brief  Destructor */
    ~Solve() override {}

    /** \brief  Print expression */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief  Modifier for linear system, before argument */
    virtual std::string mod_prefix() const {return "";}

    /** \brief  Modifier for linear system, after argument */
    virtual std::string mod_suffix() const {return "";}

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /// Can the operation be performed inplace (i.e. overwrite the result)
    casadi_int n_inplace() const override { return 1;}

    /** \brief  Propagate sparsity forward */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Get the operation */
    Operation op() const override { return Operation::OP_SOLVE;}

    /** Obtain information about function */
    Dict info() const override {
      return {{"tr", Tr}};
    }

    /// Solve another system with the same factorization
    virtual MX solve(const MX& A, const MX& B, bool tr) const = 0;

    /// Sparsity pattern for the linear system
    virtual const Sparsity& A_sp() const { return dep(1).sparsity();}

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Serialize type information */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserialize with type disambiguation */
    static MXNode* deserialize(DeserializingStream& s);

    /** \brief Deserializing constructor */
    explicit Solve(DeserializingStream& s);
  };

  /** \brief Linear solve operation with a linear solver instance

      \author Joel Andersson
      \date 2013
  */
  template<bool Tr>
  class CASADI_EXPORT LinsolCall : public Solve<Tr> {
  public:

    /** \brief  Constructor */
    LinsolCall(const MX& r, const MX& A, const Linsol& linear_solver);

    /** \brief  Destructor */
    ~LinsolCall() override {}

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief Get required length of w field */
    size_t sz_w() const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /// Linear solver (may be shared between multiple nodes)
    Linsol linsol_;

    /// Solve another system with the same factorization
    MX solve(const MX& A, const MX& B, bool tr) const override {
      return linsol_.solve(A, B, tr);
    }

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Serialize type information */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserialize with type disambiguation */
    static MXNode* deserialize(DeserializingStream& s);

    /** \brief Deserializing constructor */
    explicit LinsolCall(DeserializingStream& s);
  };

  /** \brief Linear solve with an upper triangular matrix

      \author Joel Andersson
      \date 2020
  */
  template<bool Tr>
  class CASADI_EXPORT TriuSolve : public Solve<Tr> {
  public:

    /** \brief  Constructor */
    TriuSolve(const MX& r, const MX& A);

    /** \brief  Destructor */
    ~TriuSolve() override {}

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /// Solve another system with the same factorization
    MX solve(const MX& A, const MX& B, bool tr) const override {
      return A->get_solve_triu(B, tr);
    }

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g, const std::vector<casadi_int>& arg,
      const std::vector<casadi_int>& res) const override;
  };

  /** \brief Linear solve with an upper triangular matrix

      \author Joel Andersson
      \date 2020
  */
  template<bool Tr>
  class CASADI_EXPORT TrilSolve : public Solve<Tr> {
  public:

    /** \brief  Constructor */
    TrilSolve(const MX& r, const MX& A);

    /** \brief  Destructor */
    ~TrilSolve() override {}

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /// Solve another system with the same factorization
    MX solve(const MX& A, const MX& B, bool tr) const override {
      return A->get_solve_tril(B, tr);
    }

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g, const std::vector<casadi_int>& arg,
      const std::vector<casadi_int>& res) const override;
  };

  /** \brief Linear solve with unity diagonal added

      \author Joel Andersson
      \date 2020
  */
  template<bool Tr>
  class CASADI_EXPORT SolveUnity : public Solve<Tr> {
  public:

    /** \brief  Constructor */
    SolveUnity(const MX& r, const MX& A);

    /** \brief  Destructor */
    ~SolveUnity() override {}

    /** \brief  Modifier for linear system, before argument */
    std::string mod_prefix() const override {return "(I - ";}

    /** \brief  Modifier for linear system, after argument */
    std::string mod_suffix() const override {return ")";}

    /// Sparsity pattern for the linear system
    const Sparsity& A_sp() const override;

    // Sparsity pattern of linear system, cached
    mutable Sparsity A_sp_;
  };

  /** \brief Linear solve with an upper triangular matrix, unity diagonal

      \author Joel Andersson
      \date 2020
  */
  template<bool Tr>
  class CASADI_EXPORT TriuSolveUnity : public SolveUnity<Tr> {
  public:

    /** \brief  Constructor */
    TriuSolveUnity(const MX& r, const MX& A);

    /** \brief  Destructor */
    ~TriuSolveUnity() override {}

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /// Solve another system with the same factorization
    MX solve(const MX& A, const MX& B, bool tr) const override {
      return A->get_solve_triu_unity(B, tr);
    }

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g, const std::vector<casadi_int>& arg,
      const std::vector<casadi_int>& res) const override;
  };

  /** \brief Linear solve with an upper triangular matrix

      \author Joel Andersson
      \date 2020
  */
  template<bool Tr>
  class CASADI_EXPORT TrilSolveUnity : public SolveUnity<Tr> {
  public:

    /** \brief  Constructor */
    TrilSolveUnity(const MX& r, const MX& A);

    /** \brief  Destructor */
    ~TrilSolveUnity() override {}

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /// Solve another system with the same factorization
    MX solve(const MX& A, const MX& B, bool tr) const override {
      return A->get_solve_tril_unity(B, tr);
    }

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g, const std::vector<casadi_int>& arg,
      const std::vector<casadi_int>& res) const override;
  };

} // namespace casadi

#endif // CASADI_SOLVE_HPP
