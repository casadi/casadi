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


#ifndef CASADI_XML_FILE_INTERNAL_HPP
#define CASADI_XML_FILE_INTERNAL_HPP

#include "xml_file.hpp"
#include "../function/plugin_interface.hpp"

namespace casadi {

  class CASADI_EXPORT
  XmlFileInternal : public OptionsFunctionalityNode,
                    public PluginInterface<XmlFileInternal> {
  public:
    // Constructor
    XmlFileInternal();

    // Destructor
    virtual ~XmlFileInternal();

    /** \brief  Print */
    virtual void print(std::ostream &stream) const;

    // Parse an XML file
    virtual XmlNode parse(const std::string& filename);

    // Creator function for internal class
    typedef XmlFileInternal* (*Creator)();

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Short name
    static std::string shortname() { return "xmlfile";}

    /// Infix

  };

} // namespace casadi

#endif // CASADI_XML_FILE_INTERNAL_HPP
