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


#ifndef CASADI_FMU_FUNCTION_HPP
#define CASADI_FMU_FUNCTION_HPP

#include "function.hpp"

namespace casadi {

/** \brief  Load a function from an FMU DLL
  \param name    Name assigned to the resulting function object
  \param path    Path to the unpacked FMU, where the modelDescription.xml file is located
  \param id_in   Identifiers of all the inputs
  \param id_out  Identifiers of all the outputs
  \param name_in   Names of all the inputs
  \param name_out  Names of all the outputs
  \param guid    Global unique identifier, from the corresponding modelDescription.xml file
  \param opts    Optional settings
*/
CASADI_EXPORT Function fmu_function(const std::string& name,
    const std::string& path,
    const std::vector<std::vector<casadi_int>>& id_in,
    const std::vector<std::vector<casadi_int>>& id_out,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out,
    const std::string& guid,
    const Dict& opts=Dict());

} // namespace casadi

#endif // CASADI_FMU_FUNCTION_HPP
