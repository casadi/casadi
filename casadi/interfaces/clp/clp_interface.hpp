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

#ifndef CASADI_CLP_INTERFACE_HPP
#define CASADI_CLP_INTERFACE_HPP

#include "casadi/core/function/conic_impl.hpp"
#include <casadi/interfaces/clp/casadi_conic_clp_export.h>

#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"
#include "ClpNetworkMatrix.hpp"

#include <string>

/** \defgroup plugin_Conic_clp

      Interface to Clp solver for sparse Quadratic Programs
*/

/** \pluginsection{Conic,clp} */

/// \cond INTERNAL

namespace casadi {

  struct CASADI_CONIC_CLP_EXPORT ClpMemory {
    /// Constructor
    ClpMemory();

    /// Destructor
    ~ClpMemory();
  };

  /** \brief \pluginbrief{Conic,clp}

      @copydoc Conic_doc
      @copydoc plugin_Conic_clp

      \author Attila Kozma, Joel Andersson
      \date 2012
  */
  class CASADI_CONIC_CLP_EXPORT ClpInterface : public Conic {
  public:
    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                          const std::map<std::string, Sparsity>& st) {
      return new ClpInterface(name, st);
    }

    /// Constructor using sparsity patterns
    explicit ClpInterface(const std::string& name,
                            const std::map<std::string, Sparsity>& st);

    /// Destructor
    virtual ~ClpInterface();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "clp";}

    // Initialize the solver
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual void* alloc_memory() const { return new ClpMemory();}

    /** \brief Free memory block */
    virtual void free_memory(void *mem) const { delete static_cast<ClpMemory*>(mem);}

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    // Solve the QP
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /// A documentation string
    static const std::string meta_doc;

  };
} // end namespace casadi
/// \endcond
#endif // CASADI_CLP_INTERFACE_HPP
