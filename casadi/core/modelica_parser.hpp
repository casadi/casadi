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


#ifndef CASADI_MODELICA_PARSER_HPP
#define CASADI_MODELICA_PARSER_HPP

#include "modelica_parser.hpp"
#include "shared_object.hpp"
#include "printable.hpp"

namespace casadi {

  /** Forward declaration of internal class */
  class ModelicaParserInternal;

  /** \brief Modelica parser

      Can be used for parsing Modelica files into CasADi data structures.

      \author Joris Gillis
      \date 2026

      \identifier{modelica_parser} */
  class CASADI_EXPORT ModelicaParser
    : public SharedObject,
      public SWIG_IF_ELSE(PrintableCommon, Printable<ModelicaParser>) {
  public:
    /** \brief Get type name

        \identifier{modelica_parser_type} */
    static std::string type_name() {return "ModelicaParser";}

    // Default constructor
    ModelicaParser();

    // Constructor
    ModelicaParser(const std::string& name);

    // Destructor
    ~ModelicaParser();

    /// Load a plugin dynamically
    static void load_plugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

#ifndef SWIG
    /** \brief  Access functions of the node

        \identifier{modelica_parser_access} */
    ModelicaParserInternal* operator->();

    /** \brief  Const access functions of the node

        \identifier{modelica_parser_const_access} */
    const ModelicaParserInternal* operator->() const;

    // Parse a Modelica file
    void parse(const std::string& filename);

#endif // SWIG
  };

} // namespace casadi

#endif // CASADI_MODELICA_PARSER_HPP
