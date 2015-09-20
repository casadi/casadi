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


#ifndef CASADI_COMMANDLINE_INTERFACE_HPP
#define CASADI_COMMANDLINE_INTERFACE_HPP

#include "casadi/core/function/jit_compiler_internal.hpp"
#include <casadi/solvers/casadi_jitcompiler_commandline_export.h>

/** \defgroup plugin_JitCompiler_commandline
      Interface to the JIT compiler COMMANDLINE
*/

/** \pluginsection{JitCompiler,commandline} */

/// \cond INTERNAL
namespace casadi {
  /** \brief \pluginbrief{JitCompiler,commandline}


   \author Joris Gillis
   \date 2015
   *
   @copydoc JitCompiler_doc
   @copydoc plugin_JitCompiler_commandline
   * */
  class CASADI_JITCOMPILER_COMMANDLINE_EXPORT CommandlineJitCompiler : public JitCompilerInternal {
  public:

    /** \brief Constructor */
    explicit CommandlineJitCompiler(const std::string& name);

    /** \brief Clone */
    virtual CommandlineJitCompiler* clone() const;

    /** \brief  Create a new JIT function */
    static JitCompilerInternal* creator(const std::string& name) {
      return new CommandlineJitCompiler(name);
    }

    /** \brief Destructor */
    virtual ~CommandlineJitCompiler();

    /** \brief Initialize */
    virtual void init();

    /// A documentation string
    static const std::string meta_doc;

    /// Get name of plugin
    virtual const char* plugin_name() const { return "commandline";}

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

#endif // CASADI_COMMANDLINE_INTERFACE_HPP
