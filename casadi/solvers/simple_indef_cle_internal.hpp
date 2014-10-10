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


#ifndef CASADI_SIMPLE_INDEF_CLE_INTERNAL_HPP
#define CASADI_SIMPLE_INDEF_CLE_INTERNAL_HPP

#include "../core/function/cle_internal.hpp"
#include <casadi/solvers/casadi_clesolver_simple_export.h>

/** \defgroup plugin_CleSolver_simple
 Solving the Discrete Lyapunov Equations with a regular LinearSolver

*/
/** \pluginsection{CleSolver,simple} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{CleSolver,simple}

   @copydoc CLE_doc
   @copydoc plugin_CleSolver_simple

       \author Joris Gillis
      \date 2014

  */
  class CASADI_CLESOLVER_SIMPLE_EXPORT SimpleIndefCleInternal : public CleInternal {
  public:
    /** \brief  Constructor
     * \param st \structargument{Cle}
     */
    SimpleIndefCleInternal(const CleStructure& st);

    /** \brief  Destructor */
    virtual ~SimpleIndefCleInternal();

    /** \brief  Clone */
    virtual SimpleIndefCleInternal* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual SimpleIndefCleInternal* create(const CleStructure& st) const {
        return new SimpleIndefCleInternal(st);}

    /** \brief  Create a new CLE Solver */
    static CleInternal* creator(const CleStructure& st)
    { return new SimpleIndefCleInternal(st);}

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
    /// Main implementation as MXFunction
    Function f_;

    /// State space dimension
    int n_;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_SIMPLE_INDEF_CLE_INTERNAL_HPP
