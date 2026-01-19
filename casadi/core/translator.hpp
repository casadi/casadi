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


#ifndef CASADI_TRANSLATOR_HPP
#define CASADI_TRANSLATOR_HPP

#include "function.hpp"
#include "printable.hpp"

namespace casadi {

  /** \defgroup main_translator Title
      \par

      Create a graph translator for importing/exporting computational graphs

      \generalsection{Translator}
      \pluginssection{Translator}

      \author Joris Gillis
      \date 2025

      \identifier{translator} */

  /** \defgroup translator Title
  * @copydoc main_translator
  *  @{
  */

  // Forward declaration of internal class
  class TranslatorInternal;

  /** \brief Translator

      Graph import/export for various formats (ONNX, etc.)

      \generalsection{Translator}
      \pluginssection{Translator}

      \author Joris Gillis
      \date 2025

      \identifier{translator} */
  class CASADI_EXPORT Translator
    : public SharedObject,
      public SWIG_IF_ELSE(PrintableCommon, Printable<Translator>) {
  public:
    /** \brief Get type name

        \identifier{translator_type} */
    static std::string type_name() {return "Translator";}

    /// Default constructor
    Translator();

    /// Translator factory
    explicit Translator(const std::string& name,
                       const Dict& opts=Dict());

    /// Access functions of the node
    TranslatorInternal* operator->();
    const TranslatorInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);

    /// Check if a plugin is available
    static bool has_plugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void load_plugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /// Query plugin name
    std::string plugin_name() const;

    /** \brief Load a graph from file

        \identifier{translator_load_file} */
    void load(const std::string& filename);

    /** \brief Load a CasADi Function

        \identifier{translator_load_function} */
    void load(const Function& f);

    /** \brief Set dimension for a symbolic variable

        Some formats (like ONNX) allow symbolic dimensions.
        This method allows fixing those dimensions.

        \identifier{translator_set_dimension} */
    void set_dimension(const std::string& name, casadi_int dim);

    /** \brief Create a CasADi Function from the loaded graph

        \identifier{translator_create} */
    Function create(const std::string& name);

    /** \brief Save the loaded graph/function to file

        \identifier{translator_save} */
    void save(const std::string& filename);

    /// \cond INTERNAL
#ifndef SWIG
    /** \brief  Create from node

        \identifier{translator_create_node} */
    static Translator create(TranslatorInternal* node);

    /** \brief  Create from node and initialize

        \identifier{translator_create_init} */
    static Translator create(TranslatorInternal* node, const Dict& opts);
#endif // SWIG
    /// \endcond

  };

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_translator(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_translator(const std::string& name);

  /// Get solver specific documentation
  CASADI_EXPORT std::string doc_translator(const std::string& name);

} // namespace casadi

#endif // CASADI_TRANSLATOR_HPP
