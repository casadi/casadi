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


#ifndef CASADI_SHELL_INTERFACE_HPP
#define CASADI_SHELL_INTERFACE_HPP

#include "casadi/core/function/compiler_internal.hpp"
#include <casadi/solvers/casadi_compiler_shell_export.h>

/** \defgroup plugin_Compiler_shell
      Interface to the JIT compiler SHELL
*/

/** \pluginsection{Compiler,shell} */

/// \cond INTERNAL
namespace casadi {
  /** \brief \pluginbrief{Compiler,shell}


   \author Joris Gillis
   \date 2015
   *
   @copydoc Compiler_doc
   @copydoc plugin_Compiler_shell
   * */
  class CASADI_COMPILER_SHELL_EXPORT ShellCompiler : public CompilerInternal {
  public:

    /** \brief Constructor */
    explicit ShellCompiler(const std::string& name);

    /** \brief  Create a new JIT function */
    static CompilerInternal* creator(const std::string& name) {
      return new ShellCompiler(name);
    }

    /** \brief Destructor */
    virtual ~ShellCompiler();

    /** \brief Initialize */
    virtual void init();

    /// A documentation string
    static const std::string meta_doc;

    /// Get name of plugin
    virtual const char* plugin_name() const { return "shell";}

    /// Get a function pointer for numerical evaluation
    virtual void* getFunction(const std::string& symname);
  protected:
    /// Temporary file
    std::string bin_name_;

    // Shared library handle
    typedef void* handle_t;
    handle_t handle_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_SHELL_INTERFACE_HPP
