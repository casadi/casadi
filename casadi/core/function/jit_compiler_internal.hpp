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


#ifndef CASADI_JIT_COMPILER_INTERNAL_HPP
#define CASADI_JIT_COMPILER_INTERNAL_HPP

#include "jit_compiler.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"


/// \cond INTERNAL
namespace casadi {

/** \brief NLP solver storage class

  @copydoc JitCompiler_doc
  \author Joel Andersson
  \date 2010-2013
*/
  class CASADI_EXPORT
  JitCompilerInternal : public FunctionInternal,
                      public PluginInterface<JitCompilerInternal> {

  public:
    /// Constructor
    JitCompilerInternal(const Function& f);

    /// Destructor
    virtual ~JitCompilerInternal() = 0;

    /// Initialize
    virtual void init();

    /// The Function to be JIT'ed
    Function f_;

    // Creator function for internal class
    typedef JitCompilerInternal* (*Creator)(const Function& f);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "jitcompiler";}

    /// Queery plugin name
    virtual const char* plugin_name() const = 0;

    /// Get a function pointer for numerical evaluation
    virtual void* getFunction(const std::string& symname) = 0;

  };
} // namespace casadi
/// \endcond
#endif // CASADI_JIT_COMPILER_INTERNAL_HPP
