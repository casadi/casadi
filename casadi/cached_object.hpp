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

#ifndef CACHED_OBJECT_HPP
#define CACHED_OBJECT_HPP

#include "shared_object.hpp"
#include <stack>

namespace CasADi{
    
  // Forward declaration
  class CachedObjectNode;
  
/** \brief A cached reference counted object

  \author Joel Andersson 
  \date 2012
*/
class CachedObject : public SharedObject{
  public:
    /// Default constructor
    CachedObject();
    
    /// Destructor
    ~CachedObject();
    
    /// Access a member function or object
    CachedObjectNode* operator->();
    
    /// Const access a member function or object
    const CachedObjectNode* operator->() const;

    /// Assert that the node is pointing to the right type of object
    virtual bool checkNode() const;
};

#ifndef SWIG

/// Forward declaration of the WeakRef class
class WeakRef;

/** \brief Internal class
  \author Joel Andersson 
  \date 2012
*/
class CachedObjectNode : public SharedObjectNode{
  friend class CachedObject;
  friend class WeakRef;
  public:
  
    /// Constructor
    CachedObjectNode();

    /// Destructor (marks the element as deleted in any lingering references)
    virtual ~CachedObjectNode();
    
    /// Register a weak reference
    int regRef(WeakRef* ref);
    
    /// Unregister a weak reference
    void unregRef(int weak_loc);
    
    /// List of weak references
    std::vector<WeakRef*> weak_;
    
    /// List of unused locations in the list of weak references
    std::stack<int,std::vector<int> > unused_;
};

/** \brief A weak reference to a cached object
  \author Joel Andersson 
  \date 2012
*/
class WeakRef{
  friend class CachedObjectNode;
  public:
    /// Default constructor
    WeakRef();

    /// Destructor
    ~WeakRef();

    /// Copy constructor
    WeakRef(const WeakRef& obj);
    
    /// Construct weak reference from shared object
    WeakRef(CachedObject& obj);

    /// Assignment operator (needed since copy operator overloaded)
    WeakRef& operator=(const WeakRef& obj);
    
    /// Construct weak reference from shared object through assignment operation
    WeakRef& operator=(CachedObject& obj);
    
    /// Get a shared pointer from the weak reference through an (implicit) type conversion
    template<class B>
    operator B(){
      if(isNull()){
	return B();
      } else {
	return owner_->shared_from_this<B>();
      }
    }

    /// Check if the weak object is null
    bool isNull() const { return owner_==0;}

    /// Clear the weak pointer by setting the owner to null
    void clear();
    
  private:

    /// The owner of the reference
    CachedObjectNode* owner_;
    
    /// The place in the list of weak refences
    int weak_loc_;
};

#endif // SWIG


} // namespace CasADi


#endif // OPTIONS_FUNCTIONALITY
