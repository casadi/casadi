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

#ifndef LAPACK_QR_NULLSPACE_HPP
#define LAPACK_QR_NULLSPACE_HPP

#include "../../symbolic/fx/nullspace.hpp"

namespace CasADi{

  // Forward declaration of internal class
  class LapackQRNullspaceInternal;

  /** @copydoc Nullspace_doc
      \author Joris Gillis 
      \date 2014
  */

  class LapackQRNullspace : public Nullspace{
  public:

    /// Default constructor 
    LapackQRNullspace();

    /// Default constructor 
    LapackQRNullspace(const CRSSparsity& A_sp);
    
    /// Access functions of the node.
    LapackQRNullspaceInternal* operator->();

    /// Const access functions of the node.
    const LapackQRNullspaceInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;
         
  };
  
} // namespace CasADi

#endif //LAPACK_NULLSPACE_HPP
