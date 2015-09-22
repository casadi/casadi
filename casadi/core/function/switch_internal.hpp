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


#ifndef CASADI_SWITCH_INTERNAL_HPP
#define CASADI_SWITCH_INTERNAL_HPP

#include "switch.hpp"
#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** Switch statement
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT SwitchInternal : public FunctionInternal {
    friend class Switch;
  public:

    /** \brief Constructor (generic switch) */
    SwitchInternal(const std::vector<Function>& f, const Function& f_def);

    /** \brief  clone function */
    virtual SwitchInternal* clone() const { return new SwitchInternal(*this);}

    /** \brief  Destructor */
    virtual ~SwitchInternal();

    /** \brief  Initialize */
    virtual void init();

    /** \brief  Evaluate numerically, work vectors given */
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives */
    virtual Function getDerForward(const std::string& name, int nfwd, Dict& opts);
    virtual int numDerForward() const { return 64;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives */
    virtual Function getDerReverse(const std::string& name, int nadj, Dict& opts);
    virtual int numDerReverse() const { return 64;}
    ///@}

    /** \brief  Print description */
    virtual void print(std::ostream &stream) const;

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(CodeGenerator& g) const;

    /** \brief Generate code for the body of the C function */
    virtual void generateBody(CodeGenerator& g) const;

    // Function to be evaluated for each case
    std::vector<Function> f_;

    // Default case;
    Function f_def_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_SWITCH_INTERNAL_HPP
