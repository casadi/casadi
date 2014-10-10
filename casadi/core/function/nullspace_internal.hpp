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


#ifndef CASADI_NULLSPACE_INTERNAL_HPP
#define CASADI_NULLSPACE_INTERNAL_HPP

#include "nullspace.hpp"
#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** @copydoc Nullspace_doc
      \author Joris Gillis
      \date 2014
  */
  class CASADI_CORE_EXPORT NullspaceInternal : public FunctionInternal {
  public:

    /** \brief  Constructor */
    NullspaceInternal(const Sparsity & A_sp);

    /** \brief  Destructor */
    virtual ~NullspaceInternal();

    /** \brief  Clone */
    virtual NullspaceInternal* clone() const { return new NullspaceInternal(A_sp_);}

    /** \brief  initialize */
    virtual void init();

    /** \brief  Integrate */
    virtual void evaluate();

  protected:

    /// Should the subspace be dense?
    bool dense_;

    /// Storage for input structure
    Sparsity A_sp_;

    /// number of rows
    int m_;

    /// number of columns
    int n_;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_NULLSPACE_INTERNAL_HPP
