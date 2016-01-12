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


#ifndef CASADI_NEWTON_HPP
#define CASADI_NEWTON_HPP

#include "casadi/core/function/rootfinder.hpp"

#include <casadi/solvers/casadi_rootfinder_newton_export.h>

/** \defgroup plugin_Rootfinder_newton
     Implements simple newton iterations to solve an implicit function.
*/

/** \pluginsection{Rootfinder,newton} */

/// \cond INTERNAL
namespace casadi {

  // Memory
  struct CASADI_ROOTFINDER_NEWTON_EXPORT NewtonMemory : public Memory {
    /// Constructor
    NewtonMemory();

    /// Destructor
    virtual ~NewtonMemory() {}

    /// Stats
    const char* return_status;
    int iter;
  };

  /** \brief \pluginbrief{Rootfinder,newton}

      @copydoc Rootfinder_doc
      @copydoc plugin_Rootfinder_newton

      \author Joris Gillis
      \date 2012
  */
  class CASADI_ROOTFINDER_NEWTON_EXPORT Newton : public Rootfinder {
  public:
    /** \brief  Constructor */
    explicit Newton(const std::string& name, const Function& f);

    /** \brief  Destructor */
    virtual ~Newton();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "newton";}

    /** \brief  Create a new Rootfinder */
    static Rootfinder* creator(const std::string& name, const Function& f) {
      return new Newton(name, f);
    }

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /** \brief Allocate memory block */
    virtual Memory* memory() const { return new NewtonMemory();}

    /// Solve the system of equations and calculate derivatives
    virtual void eval(Memory& mem, const double** arg, double** res,
                      int* iw, double* w) const;

    /// A documentation string
    static const std::string meta_doc;

  protected:
    /// Maximum number of Newton iterations
    int max_iter_;

    /// Absolute tolerance that should be met on residual
    double abstol_;

    /// Absolute tolerance that should be met on step
    double abstolStep_;

    /// If true, each iteration will be printed
    bool print_iteration_;

    /// Print iteration header
    void printIteration(std::ostream &stream) const;

    /// Print iteration
    void printIteration(std::ostream &stream, int iter,
                        double abstol, double abstolStep) const;
  };

} // namespace casadi
/// \endcond
#endif // CASADI_NEWTON_HPP
