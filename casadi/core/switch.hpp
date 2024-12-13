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


#ifndef CASADI_SWITCH_HPP
#define CASADI_SWITCH_HPP

#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** Switch statement
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT Switch : public FunctionInternal {
  public:

    /** \brief Constructor (generic switch)

        \identifier{1ko} */
    Switch(const std::string& name,
                   const std::vector<Function>& f, const Function& f_def);

    /** \brief  Destructor

        \identifier{1kp} */
    ~Switch() override;

    /** \brief Get type name

        \identifier{1kq} */
    std::string class_name() const override {return "Switch";}

    ///@{
    /** \brief Number of function inputs and outputs

        \identifier{1kr} */
    size_t get_n_in() override;
    size_t get_n_out() override;
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs

        \identifier{1ks} */
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    /// @}

    /** \brief  Initialize

        \identifier{1kt} */
    void init(const Dict& opts) override;

    /** \brief  Evaluate numerically, work vectors given

        \identifier{1ku} */
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    /** \brief  evaluate symbolically while also propagating directional derivatives

        \identifier{1kv} */
    int eval_sx(const SXElem** arg, SXElem** res,
                casadi_int* iw, SXElem* w, void* mem,
                bool always_inline, bool never_inline) const override;

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives

        \identifier{1kw} */
    bool has_forward(casadi_int nfwd) const override { return true;}
    Function get_forward(casadi_int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives

        \identifier{1kx} */
    bool has_reverse(casadi_int nadj) const override { return true;}
    Function get_reverse(casadi_int nadj, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    /** \brief  Print description

        \identifier{1ky} */
    void disp_more(std::ostream& stream) const override;

    /** \brief Generate code for the declarations of the C function

        \identifier{1kz} */
    void codegen_declarations(CodeGenerator& g) const override;

    /** \brief Is codegen supported?

        \identifier{1l0} */
    bool has_codegen() const override { return true;}

    /** \brief Generate code for the body of the C function

        \identifier{1l1} */
    void codegen_body(CodeGenerator& g) const override;

    // Function to be evaluated for each case
    std::vector<Function> f_;

    // Default case;
    Function f_def_;

    // Sparsity projection needed?
    bool project_in_, project_out_;

    /** \brief Serialize an object without type information

        \identifier{1l2} */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize without type information

        \identifier{1l3} */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new Switch(s); }

    /** Obtain information about node */
    Dict info() const override;

    // Get all embedded functions, recursively
    void find(std::map<FunctionInternal*, Function>& all_fun, casadi_int max_depth) const override;

  protected:
    /** \brief Deserializing constructor

        \identifier{1l4} */
    explicit Switch(DeserializingStream& s);
  };

} // namespace casadi
/// \endcond

#endif // CASADI_SWITCH_HPP
