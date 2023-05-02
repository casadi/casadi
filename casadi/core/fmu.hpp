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


#ifndef CASADI_FMU_HPP
#define CASADI_FMU_HPP

#ifdef WITH_FMU

#include "shared_object.hpp"
#include "printable.hpp"

/// \cond INTERNAL

namespace casadi {

// Forward declarations
class DaeBuilderInternal;
struct FmuMemory;
struct InputStruct;

// Forward declaration of internal class
class Fmu2;

/// Which C API
enum class FmuApi {FMI2, NUMEL};

/// Convert to string
CASADI_EXPORT std::string to_string(FmuApi v);

/** \brief Interface to binary FMU

    Can be shared between multiple CasADi functions.

    \author Joel Andersson
    \date 2023
*/
class CASADI_EXPORT Fmu
  : public SharedObject,
    public SWIG_IF_ELSE(PrintableCommon, Printable<Fmu>) {
 public:
  /** \brief Get type name */
  static std::string type_name() {return "Fmu";}

  /// Default constructor
  Fmu();

  /// Importer factory
  Fmu(const std::string& name, FmuApi api, const DaeBuilderInternal* dae,
    const std::vector<std::string>& scheme_in,
    const std::vector<std::string>& scheme_out,
    const std::map<std::string, std::vector<size_t>>& scheme,
    const std::vector<std::string>& aux);

  /** \brief Name of the instance */
  const std::string& name() const;

  /// Access functions of the node
  Fmu2* operator->();
  const Fmu2* operator->() const;

  // Index lookup for input
  size_t index_in(const std::string& n) const;

  // Index lookup for output
  size_t index_out(const std::string& n) const;

  /// Does the interface support analytic derivatives?
  bool has_ad() const;

  /** \brief Initalize memory block */
  int init_mem(FmuMemory* m) const;

  // Free FMU instance
  void free_instance(void* c) const;
};

} // namespace casadi

/// \endcond

#endif  // WITH_FMU

#endif // CASADI_FMU_HPP
