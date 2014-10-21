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


#ifndef CASADI_CLE_INTERNAL_HPP
#define CASADI_CLE_INTERNAL_HPP

#include "cle_solver.hpp"
#include "cle_internal.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL

namespace casadi {

  /** \brief Internal storage for CleSolver related data

      @copydoc CLE_doc
     \author Joris Gillis
      \date 2014
  */
  class CASADI_CORE_EXPORT
  CleInternal : public FunctionInternal,
                 public PluginInterface<CleInternal> {
  public:
    /** \brief  Constructor
     *  \param st \structargument{Cle}
     */
    CleInternal(const CleStructure& st,
                 int nrhs=1, bool transp=false);

    /** \brief  Destructor */
    virtual ~CleInternal()=0;

    /** \brief  Clone */
    virtual CleInternal* clone() const=0;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual CleInternal* create(const CleStructure& st) const = 0;

    /** \brief  Print solver statistics */
    virtual void printStats(std::ostream &stream) const {}

    /** \brief  evaluate */
    virtual void evaluate()=0;

    /** \brief  Initialize */
    virtual void init();

    /** \brief Generate a function that calculates \a nfwd forward derivatives
     and \a nadj adjoint derivatives
     */
    virtual Function getDerivative(int nfwd, int nadj)=0;

    /// Structure of Cle
    CleStructure st_;

    /// Sparsities of A
    Sparsity A_;

    /// Sparsities of V
    Sparsity V_;

    /// Assume positive definiteness of P_i
    bool pos_def_;

    /// Throw an error when system is unstable
    bool error_unstable_;

    /// Margin for instability detection
    double eps_unstable_;

    /// Number of right hand sides
    int nrhs_;

    /// Tranpose the system?
    bool transp_;

    // Creator function for internal class
    typedef CleInternal* (*Creator)(const CleStructure& st);

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "cle";}
  };

} // namespace casadi
/// \endcond
#endif // CASADI_CLE_INTERNAL_HPP
