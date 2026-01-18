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


#ifndef CASADI_TRANSLATOR_INTERNAL_HPP
#define CASADI_TRANSLATOR_INTERNAL_HPP

#include "translator.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"


/// \cond INTERNAL
namespace casadi {

/** \brief Translator internal class

  @copydoc Translator_doc
  \author Joris Gillis
  \date 2025

    \identifier{translator_internal} */
  class CASADI_EXPORT
  TranslatorInternal : public SharedObjectInternal,
                       public PluginInterface<TranslatorInternal> {

  public:
    /// Constructor
    explicit TranslatorInternal(const std::string& name);

    /// Destructor
    ~TranslatorInternal() override;

    /** \brief Get type name

        \identifier{translator_internal_type} */
    std::string class_name() const override { return "TranslatorInternal";}

    /** \brief Print

        \identifier{translator_internal_print} */
    void disp(std::ostream& stream, bool more) const override;

    // Creator function for internal class
    typedef TranslatorInternal* (*Creator)();

    /** \brief Construct

        Prepares the translator for use

        \identifier{translator_internal_construct} */
    void construct(const Dict& opts);

    ///@{
    /** \brief Options

        \identifier{translator_internal_options} */
    static const Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /** \brief Initialize

        \identifier{translator_internal_init} */
    virtual void init(const Dict& opts);

    virtual void finalize() {}

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "translator";}

    /// Query plugin name
    const char* plugin_name() const override = 0;

    /** \brief Load a graph from file

        \identifier{translator_internal_load_file} */
    virtual void load(const std::string& filename) = 0;

    /** \brief Load a CasADi Function

        \identifier{translator_internal_load_function} */
    virtual void load(const Function& f) = 0;

    /** \brief Set dimension for a symbolic variable

        \identifier{translator_internal_set_dimension} */
    virtual void set_dimension(const std::string& name, casadi_int dim) = 0;

    /** \brief Create a CasADi Function from the loaded graph

        \identifier{translator_internal_create} */
    virtual Function create(const std::string& name) = 0;

    /** \brief Save the loaded graph/function to file

        \identifier{translator_internal_save} */
    virtual void save(const std::string& filename) = 0;

  protected:
    /// Translator name
    std::string name_;

    /** \brief  Verbose -- for debugging purposes

        \identifier{translator_internal_verbose} */
    bool verbose_;
  };

} // namespace casadi
/// \endcond
#endif // CASADI_TRANSLATOR_INTERNAL_HPP
