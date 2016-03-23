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
#include "../casadi_file.hpp"

namespace casadi {

  // Forward declaration of internal class
  class CompilerInternal;

  /** \brief Compiler

      Just-in-time compilation of code

      \generalsection{Compiler}
      \pluginssection{Compiler}

      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT Compiler : public SharedObject {
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
    static bool test_cast(const SharedObjectNode* ptr);

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
    void* get_function(const std::string& symname);

    /// Get meta information
    const ParsedFile& meta() const;
#endif // SWIG
  };

  // Forward declaration of internal class
  class LibraryInternal;

  /** \brief Library, either just-in-time compiled or dynamically loaded
  */
  class CASADI_EXPORT Library : public SharedObject {
  public:
    /// Default constructor
    Library();

    // Constructor, DLL
    explicit Library(const std::string& bin_name);

    // Constructor, JIT
    explicit Library(const Compiler& compiler);

    /// Access functions of the node
    LibraryInternal* operator->();
    const LibraryInternal* operator->() const;

    // Check if symbol exists
    bool has(const std::string& sym) const;

#ifndef SWIG
    // Dummy type
    signal_t get(const std::string& sym);

    // Get meta
    const ParsedFile& meta() const;
#endif // SWIG
  };

} // namespace casadi

#endif // CASADI_COMPILER_HPP

