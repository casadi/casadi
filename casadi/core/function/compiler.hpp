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


#ifndef CASADI_COMPILER_HPP
#define CASADI_COMPILER_HPP

#include "function.hpp"


namespace casadi {

  class CompilerInternal;

  /** \brief Compiler

      Just-in-time compilation of code

      \generalsection{Compiler}
      \pluginssection{Compiler}

      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT Compiler : public OptionsFunctionality {
  public:

    /// Default constructor
    Compiler();

    /// Compiler factory (new syntax, includes initialization)
    explicit Compiler(const std::string& name,
                         const std::string& compiler,
                         const Dict& opts=Dict());

    /// Access functions of the node
    CompilerInternal* operator->();
    const CompilerInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

    /// Check if a plugin is available
    static bool hasPlugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /// Query plugin name
    std::string plugin_name() const;

#ifndef SWIG
    /// Get a function pointer for numerical evaluation
    void* getFunction(const std::string& symname);
#endif // SWIG
  };

} // namespace casadi

#endif // CASADI_COMPILER_HPP

