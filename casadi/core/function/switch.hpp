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

    /** \brief Constructor (generic switch) */
    Switch(const std::string& name,
                   const std::vector<Function>& f, const Function& f_def);

    /** \brief  Destructor */
    virtual ~Switch();

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in();
    virtual size_t get_n_out();
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int i);
    virtual Sparsity get_sparsity_out(int i);
    /// @}

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /** \brief  Evaluate numerically, work vectors given */
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives */
    virtual Function get_forward(const std::string& name, int nfwd,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names,
                                 const Dict& opts);
    virtual int get_n_forward() const { return 64;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives */
    virtual Function get_reverse(const std::string& name, int nadj,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names,
                                 const Dict& opts);
    virtual int get_n_reverse() const { return 64;}
    ///@}

    /** \brief  Print description */
    virtual void print(std::ostream &stream) const;

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(CodeGenerator& g) const;

    /** \brief Is codegen supported? */
    virtual bool has_codegen() const { return true;}

    /** \brief Generate code for the body of the C function */
    virtual void generateBody(CodeGenerator& g) const;

    // Function to be evaluated for each case
    std::vector<Function> f_;

    // Default case;
    Function f_def_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_SWITCH_HPP
