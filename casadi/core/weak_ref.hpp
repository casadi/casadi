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


#ifndef CASADI_WEAK_REF_HPP
#define CASADI_WEAK_REF_HPP

#include "shared_object.hpp"


/// \cond INTERNAL
namespace casadi {

  // Forward declaration
  class WeakRefInternal;


  /** \brief Weak reference type
      A weak reference to a SharedObject
      \author Joel Andersson
      \date 2013
  */
  class CASADI_EXPORT WeakRef : public SharedObject {
  public:
    friend class SharedObjectNode;

    /** \brief Default constructor */
    WeakRef(int dummy=0);

    /** \brief Construct from a shared object (also implicit type conversion) */
    WeakRef(SharedObject shared);

    /** \brief Get a shared (owning) reference */
    SharedObject shared();

    /** \brief Check if alive */
    bool alive() const;

    /** \brief  Access functions of the node */
    WeakRefInternal* operator->();

    /** \brief  Const access functions of the node */
    const WeakRefInternal* operator->() const;

#ifndef SWIG
  private:
    /** \brief Construct from a shared object (internal) */
    explicit WeakRef(SharedObjectNode* raw);

    /** \brief The shared object has been deleted */
    void kill();
#endif // SWIG
 };

#ifndef SWIG
  class CASADI_EXPORT WeakRefInternal : public SharedObjectNode {
  public:
    // Constructor
    WeakRefInternal(SharedObjectNode* raw);

    // Destructor
    ~WeakRefInternal();

    // Raw pointer to the cached object
    SharedObjectNode* raw_;
  };

#endif // SWIG

} // namespace casadi

/// \endcond


#endif // CASADI_WEAK_REF_HPP
