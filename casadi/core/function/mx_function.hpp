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


#ifndef CASADI_MX_FUNCTION_HPP
#define CASADI_MX_FUNCTION_HPP

#include <set>
#include <iostream>
#include "../mx/mx.hpp"
#include "sx_function.hpp"

namespace casadi {

  /** \brief  Forward declaration of internal class */
  class MXFunctionInternal;

  /** \brief  General function mapping from/to MX
      \author Joel Andersson
      \date 2010-2015
  */
  class CASADI_EXPORT MXFunction : public Function {
  public:

    /** \brief  Default constructor */
    MXFunction();

    /** \brief  Attempt to form an MXFunction out of an Function */
    explicit MXFunction(const Function& function);

    /** \brief Construct from vectors (new syntax, includes initialization) */
    MXFunction(const std::string& name, const std::vector<MX>& arg,
               const std::vector<MX>& res, const Dict& opts=Dict());

    /** \brief Construct from vectors (new syntax, includes initialization) */
    MXFunction(const std::string& name, const std::pair< MXDict, std::vector<std::string> >& arg,
               const std::vector<MX>& res, const Dict& opts=Dict());

    /** \brief Construct from vectors (new syntax, includes initialization) */
    MXFunction(const std::string& name, const std::vector<MX>& arg,
               const std::pair< MXDict, std::vector<std::string> >& res, const Dict& opts=Dict());

    /** \brief Construct from vectors (new syntax, includes initialization) */
    MXFunction(const std::string& name, const std::pair< MXDict, std::vector<std::string> >& arg,
               const std::pair< MXDict, std::vector<std::string> >& res, const Dict& opts=Dict());

#ifndef SWIG
#ifdef USE_CXX11
    /** \brief Construct from initializer lists (new syntax, includes initialization) */
    MXFunction(const std::string& name,
               std::initializer_list<MX> arg,
               std::initializer_list<MX> res,
               const Dict& opts=Dict());

    /** \brief Construct from vector & initializer list (new syntax, includes initialization) */
    MXFunction(const std::string& name,
               std::vector<MX> arg,
               std::initializer_list<MX> res,
               const Dict& opts=Dict());

    /** \brief Construct from initializer list & vector (new syntax, includes initialization) */
    MXFunction(const std::string& name,
               std::initializer_list<MX> arg,
               std::vector<MX> res,
               const Dict& opts=Dict());
#endif // USE_CXX11
#endif // SWIG

    /// \cond INTERNAL
    /** \brief  Access functions of the node */
    MXFunctionInternal* operator->();

    /** \brief  Const access functions of the node */
    const MXFunctionInternal* operator->() const;

#ifndef SWIG
    /** \brief Called from constructor */
    void construct(const std::string& name, const std::vector<MX>& arg,
                   const std::vector<MX>& res, const Dict& opts,
                   const std::vector<std::string>& ischeme=std::vector<std::string>(),
                   const std::vector<std::string>& oscheme=std::vector<std::string>());
#endif // SWIG
    /// \endcond

    /** \brief Number of nodes in the algorithm */
    int countNodes() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

#ifdef WITH_DEPRECATED_FEATURES
    ///@{
    /** \brief Jacobian expression */
    MX jac(int iind=0, int oind=0, bool compact=false, bool symmetric=false);
    MX jac(const std::string & iname, int oind=0, bool compact=false, bool symmetric=false)
    { return jac(index_in(iname), oind, compact, symmetric); }
    MX jac(int iind, const std::string & oname, bool compact=false, bool symmetric=false)
    { return jac(iind, index_out(oname), compact, symmetric); }
    MX jac(const std::string & iname, const std::string & oname,
           bool compact=false, bool symmetric=false)
    { return jac(index_in(iname), index_out(oname), compact, symmetric); }
    ///@}

    ///@{
    /** \brief Gradient expression */
    MX grad(int iind=0, int oind=0);
    MX grad(const std::string & iname, int oind=0) { return grad(index_in(iname), oind); }
    MX grad(int iind, const std::string & oname) { return grad(iind, index_out(oname)); }
    MX grad(const std::string & iname, const std::string & oname)
    { return grad(index_in(iname), index_out(oname)); }
    ///@}

    ///@{
    /** \brief Tangent expression */
    MX tang(int iind=0, int oind=0);
    MX tang(const std::string & iname, int oind=0) { return tang(index_in(iname), oind); }
    MX tang(int iind, const std::string & oname) { return tang(iind, index_out(oname)); }
    MX tang(const std::string & iname, const std::string & oname)
    { return tang(index_in(iname), index_out(oname)); }
    ///@}

    /** \brief Expand the matrix valued graph into a scalar valued graph */
    SXFunction expand(const std::vector<SX>& inputv = std::vector<SX>());
#endif // WITH_DEPRECATED_FEATURES

    /** \brief Get all the free variables of the function */
    std::vector<MX> getFree() const;

    /// \cond INTERNAL
    /** \brief Extract the functions needed for the Lifted Newton method */
    void generateLiftingFunctions(MXFunction& SWIG_OUTPUT(vdef_fcn),
                                  MXFunction& SWIG_OUTPUT(vinit_fcn));
    /// \endcond

    /** \brief Get the corresponding matrix type */
    typedef MX MatType;
  };

} // namespace casadi


#endif // CASADI_MX_FUNCTION_HPP

