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


#ifndef CASADI_LACE_MODELICA_HPP
#define CASADI_LACE_MODELICA_HPP

/** \defgroup plugin_ModelicaParser_lacemodelica Title
    \par

 * ModelicaParser using LaceModelica

    \identifier{lace_modelica_plugin} */

/** \pluginsection{ModelicaParser,lacemodelica} */

/// \cond INTERNAL
#include "casadi/core/modelica_parser_internal.hpp"
#include <casadi/interfaces/lacemodelica/casadi_modelicaparser_lacemodelica_export.h>

namespace casadi {

  /** \brief \pluginbrief{ModelicaParser,lacemodelica}
   * @copydoc ModelicaParser_doc
   * @copydoc plugin_ModelicaParser_lacemodelica
   */
  class CASADI_MODELICAPARSER_LACEMODELICA_EXPORT LaceModelica : public ModelicaParserInternal {
  public:

    // Create a Modelica parser
    LaceModelica();

    /** \brief  Create a new ModelicaParser */
    static ModelicaParserInternal* creator()
    { return new LaceModelica();}

    // Get name of the plugin
    const char* plugin_name() const override { return "lacemodelica";}

    // Get name of the class
    std::string class_name() const override { return "LaceModelica";}

    // Parse a Modelica file and generate output in the given directory
    void parse(const std::string& filename, const std::string& output_dir) override;

    // Destructor
    ~LaceModelica() override;

    /// A documentation string
    static const std::string meta_doc;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_LACE_MODELICA_HPP
