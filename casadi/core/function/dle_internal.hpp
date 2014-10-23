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


#ifndef CASADI_DLE_INTERNAL_HPP
#define CASADI_DLE_INTERNAL_HPP

#include "dle_solver.hpp"
#include "dle_internal.hpp"
#include "lr_dle_internal.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL

namespace casadi {

  /** \brief Internal storage for DleSolver related data

      @copydoc DLE_doc
     \author Joris Gillis
      \date 2014
  */
  class CASADI_CORE_EXPORT
  DleInternal : public FunctionInternal,
                 public PluginInterface<DleInternal> {
  public:
    /** \brief  Constructor
     *  \param st \structargument{Dle}
     */
    DleInternal(const DleStructure& st,
                 int nrhs=1, bool transp=false);

    /** \brief  Destructor */
    virtual ~DleInternal()=0;

    /** \brief  Clone */
    virtual DleInternal* clone() const=0;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual DleInternal* create(const DleStructure& st) const = 0;

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

    /// Problem structure
    DleStructure st_;

    /// Sparsity of A
    Sparsity A_;

    /// List of sparsities of V_i
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
    typedef DleInternal* (*Creator)(const DleStructure& st);

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "dle";}

    /// Get the resulting sparsity
    static Sparsity getSparsity(const DleStructure& st);

  };

} // namespace casadi
/// \endcond
#endif // CASADI_DLE_INTERNAL_HPP
