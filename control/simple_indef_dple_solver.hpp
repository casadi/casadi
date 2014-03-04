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

#ifndef SIMPLE_INDEF_DPLE_SOLVER_HPP
#define SIMPLE_INDEF_DPLE_SOLVER_HPP

#include "dple_solver.hpp"

namespace CasADi{

  /// Forward declaration of internal class
  class SimpleIndefDpleInternal;

  /**  @copydoc DPLE_doc
  
       Uses Periodic Schur Decomposition (simple) and does not assume positive definiteness.
       Based on Periodic Lyapunov equations: some applications and new algorithms. Int. J. Control, vol. 67, pp. 69-87, 1997. 
  
       \author Joris gillis
      \date 2014
      
  */
  class SimpleIndefDpleSolver : public DpleSolver {
  public:
    /// Default constructor
    SimpleIndefDpleSolver();
    
    /** \brief  Constructor
     *  \param[in] A  List of sparsities of A_i 
     *  \param[in] V  List of sparsities of V_i 
     */
    explicit SimpleIndefDpleSolver(const std::vector< Sparsity > & A, const std::vector< Sparsity > &V);
    
    /// Access functions of the node
    SimpleIndefDpleInternal* operator->();

    /// Access functions of the node
    const SimpleIndefDpleInternal* operator->() const;
 
    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;
  
    /// Static creator function
    #ifdef SWIG
    %callback("%s_cb");
    #endif
    static DpleSolver creator(const std::vector< Sparsity > & A, const std::vector< Sparsity > &V){ return SimpleIndefDpleSolver(A,V);}
    #ifdef SWIG
    %nocallback;
    #endif
  
  };

} // namespace CasADi

#endif // DPLE_SOLVER_HPP
