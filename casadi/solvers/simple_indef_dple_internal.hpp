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


#ifndef CASADI_SIMPLE_INDEF_DPLE_INTERNAL_HPP
#define CASADI_SIMPLE_INDEF_DPLE_INTERNAL_HPP

#include "../core/function/dple_internal.hpp"
#include <casadi/solvers/casadi_dplesolver_simple_export.h>

/** \defgroup plugin_DpleSolver_simple
 Solving the Discrete Periodic Lyapunov Equations with a regular LinearSolver

*/
/** \pluginsection{DpleSolver,simple} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{DpleSolver,simple}

   @copydoc DPLE_doc
   @copydoc plugin_DpleSolver_simple

       \author Joris Gillis
      \date 2014

  */
  class CASADI_DPLESOLVER_SIMPLE_EXPORT SimpleIndefDpleInternal : public DpleInternal,
    public Wrapper<SimpleIndefDpleInternal> {
  public:
    /** \brief  Constructor
     * \param st \structargument{Dple}
     */
    SimpleIndefDpleInternal(const DpleStructure& st);

    /** \brief  Destructor */
    virtual ~SimpleIndefDpleInternal();

    /** \brief  Clone */
    virtual SimpleIndefDpleInternal* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual SimpleIndefDpleInternal* create(const DpleStructure& st) const {
        return new SimpleIndefDpleInternal(st);}

    /** \brief  Create a new DPLE Solver */
    static DpleInternal* creator(const DpleStructure& st)
    { return new SimpleIndefDpleInternal(st);}

    /** \brief  Print solver statistics */
    virtual void printStats(std::ostream &stream) const {}

    /** \brief  evaluate */
    virtual void evaluate();

    /** \brief  Initialize */
    virtual void init();

    /** \brief Generate a function that calculates \a nfwd forward derivatives
     and \a nadj adjoint derivatives
    */
    virtual Function getDerivative(int nfwd, int nadj);

    /// A documentation string
    static const std::string meta_doc;

  private:

    /// State space dimension
    int n_;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_SIMPLE_INDEF_DPLE_INTERNAL_HPP
