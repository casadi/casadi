/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef DSDP_INTERNAL_HPP
#define DSDP_INTERNAL_HPP

#include "symbolic/function/sdp_solver_internal.hpp"

#include <dsdp5.h>

/// \cond INTERNAL
namespace CasADi{

  /** \brief Internal class for DSDPSolver
   * 
   @copydoc SDPSolver_doc
   * */
  class DSDPInternal : public SDPSolverInternal {
    friend class DSDPSolver;
  public:

    /** \brief Constructor */
    explicit DSDPInternal(const std::vector<Sparsity> &st);

    /** \brief Clone */
    virtual DSDPInternal* clone() const;
  
    /** \brief Destructor */
    virtual ~DSDPInternal();

    /** \brief Initialize */
    virtual void init();
  
    /** \brief Solve the SDP */
    virtual void evaluate();
  
    // Data members (all public)
    DSDP dsdp_;
    SDPCone sdpcone_;
    LPCone lpcone_;
    BCone bcone_;
    
    /** Get termination reason from flag */
    static const char* terminationReason(int flag);

    /** Get solution type from flag */
    static const char* solutionType(int flag);
    
    
    std::vector< std::vector< std::vector<int> > > pattern_;
    std::vector< std::vector< std::vector<double> > > values_;
    
    /// Temporary work vector of size n*(n+1)/2
    std::vector< std::vector<double> > store_X_;
    std::vector< std::vector<double> > store_P_;
    
    /// Mapping to get [A LBA]'
    Function mappingA_;
  };

} // namespace CasADi
/// \endcond

#endif //DSDP_INTERNAL_HPP
