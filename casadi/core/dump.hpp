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


#ifndef CASADI_DUMP_HPP
#define CASADI_DUMP_HPP

#include "mx_node.hpp"

/// \cond INTERNAL
namespace casadi {
  /** \brief Dump

      \author Joris Gillis
      \date 2026

      \identifier{2fa} */
  class CASADI_EXPORT Dump : public MXNode {
  public:

    /// Constructor
    Dump(const MX& x, const std::string& base_filename,
         const std::string& dir, const std::string& format, bool verbose);

    /// Destructor
    ~Dump() override {}

    /// Validate options (called from constructors)
    void finalize();

    /// Reset the dump counter
    void reset_dump_count();

    /** \brief  Evaluate symbolically (MX)

        \identifier{2fb} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res,
        const std::vector<bool>& unique={}) const override;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief Evaluate the MX node on a const/linear/nonlinear partition

        \identifier{2fc} */
    void eval_linear(const std::vector<std::array<MX, 3> >& arg,
            std::vector<std::array<MX, 3> >& res) const override {
        eval_linear_rearrange(arg, res);
    }

    /** \brief  Propagate sparsity forward

        \identifier{2fd} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{2fe} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Add dependencies for code generation

        \identifier{2ff} */
    void add_dependency(CodeGenerator& g) const override;

    /** \brief Generate code for the operation

        \identifier{2fg} */
    void generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res,
                          const std::vector<bool>& arg_is_ref,
                          std::vector<bool>& res_is_ref,
                          bool prefer_inline=false) const override;

    /** \brief  Print expression

        \identifier{2fh} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{2fi} */
    casadi_int op() const override { return OP_DUMP;}

    /// Can the operation be performed inplace (i.e. overwrite the result)
    casadi_int n_inplace() const override { return 1;}

    /** \brief Serialize an object without type information

        \identifier{2fj} */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Deserialize without type information

        \identifier{2fk} */
    static MXNode* deserialize(DeserializingStream& s) { return new Dump(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{2fl} */
    explicit Dump(DeserializingStream& s);

  private:
    std::string base_filename_;
    std::string dir_;
    std::string format_;
    bool verbose_;

    // Counter for unique dump filenames (per-node)
#ifdef CASADI_WITH_THREAD
    mutable std::atomic<casadi_int> dump_count_{0};
#else
    mutable casadi_int dump_count_{0};
#endif // CASADI_WITH_THREAD
  };


} // namespace casadi

/// \endcond

#endif // CASADI_DUMP_HPP
