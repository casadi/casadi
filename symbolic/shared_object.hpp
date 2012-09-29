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

#ifndef SHARED_OBJECT_HPP
#define SHARED_OBJECT_HPP

#include "printable_object.hpp"
#include "casadi_exception.hpp"
#include <map>
#include <vector>

namespace CasADi{

// Forward declaration of internal class
class SharedObjectNode;

/** \brief SharedObject implements a reference counting framework simular for effient and easily-maintained memory management.
  
  To use the class, both the SharedObject class (the public class), and the SharedObjectNode class (the internal class)
  must be inherited from. It can be done in two different files and together with memory management, this approach
  provides a clear destinction of which methods of the class are to be considered "public", i.e. methods for public use 
  that can be considered to remain over time with small changes, and the internal memory.
  
  When interfacing a software, which typically includes including some header file, this is best done only in the file 
  where the internal class is defined, to avoid polluting the global namespace and other side effects.

  The default constructor always means creating a null pointer to an internal class only. To allocate an internal class
  (this works only when the internal class isn't abstract), use the constructor with arguments.
  
  The copy constructor and the assignment operator perform shallow copies only, to make a deep copy you must use the
  clone method explictly. This will give a shared pointer instance.
  
  In an inheritance hierarchy, you can cast down automatically, e.g. (SXFunction is a child class of FX):
  SXFunction derived(...);
  FX base = derived;
  
  To cast up, use the shared_cast template function, which works analogously to dynamic_cast, static_cast, const_cast etc, e.g.:
  SXFunction derived(...);
  FX base = derived;
  SXFunction derived_from_base = shared_cast<SXFunction>(base);
  
  A failed shared_cast will result in a null pointer (cf. dynamic_cast)

  \author Joel Andersson 
  \date 2010	
*/
class SharedObject : public PrintableObject{
#ifndef SWIG
  template<class B> friend B shared_cast(SharedObject& A);
  template<class B> friend const B shared_cast(const SharedObject& A);
#endif // SWIG
  
  public:
#ifndef SWIG
    /// Default constructor
    SharedObject();
    
    /// Copy constructor (shallow copy)
    SharedObject(const SharedObject& ref);

    /// Deep copy
    SharedObject clone() const;
    
    /// Destructor
    ~SharedObject();
    
    /// Assignment operator
    SharedObject& operator=(const SharedObject& ref);
    
    /// Assign the node to a node class pointer (or null)
    void assignNode(SharedObjectNode* node);
    
    /// Assign the node to a node class pointer without reference counting: inproper use will cause memory leaks!
    void assignNodeNoCount(SharedObjectNode* node);
    
    /// Get a const pointer to the node
    const SharedObjectNode* get() const;

    /// Get a pointer to the node
    SharedObjectNode* get();

    /// Get the reference count
    int getCount() const;
    
    /// Swap content with another instance
    void swap(SharedObject& other);

    /// Access a member function or object
    SharedObjectNode* operator->();

    /// Const access a member function or object
    const SharedObjectNode* operator->() const;

    /// Print a representation of the object
    virtual void repr(std::ostream &stream) const;

    /// Print a destription of the object
    virtual void print(std::ostream &stream=std::cout) const;
    
    #endif // SWIG
    
    /// Initialize the object: more documentation in the node class (SharedObjectNode and derived classes)
    void init();

    /// Is initialized?
    bool isInit() const;
    
    /// Assert that it is initialized
    void assertInit() const;
    
    /// Is a null pointer?
    bool isNull() const;

    /// Assert that the node is pointing to the right type of object
    virtual bool checkNode() const;
    
    //@{
    /// If there are other references to the object, then make a deep copy of it and point to this new object
    void makeUnique(bool clone_members=true);
#ifndef SWIG
    void makeUnique(std::map<SharedObjectNode*,SharedObject>& already_copied, bool clone_members=true);
    //@}
    
  private:
    SharedObjectNode *node;
    void count_up(); // increase counter of the node
    void count_down(); // decrease counter of the node
#endif // SWIG
};

#ifndef SWIG

/// Internal class for the reference counting framework, see comments on the public class.
class SharedObjectNode{
  friend class SharedObject;
  public:
  
  /// Default constructor
  SharedObjectNode();

  /// Copy constructor
  SharedObjectNode(const SharedObjectNode& node);
  
  /// Assignment operator
  SharedObjectNode& operator=(const SharedObjectNode& node);
  
  /// Destructor
  virtual ~SharedObjectNode() = 0;

  /// Make a deep copy of the instance  
  virtual SharedObjectNode* clone() const=0;

  /// Deep copy data members
  virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

  /// Get the reference count
  int getCount() const;

  /// Initialize the object
  virtual void init();
  
  /// Check if the object has been initialized
  bool isInit() const;
  
  /// Assert that the object has been initialized
  void assertInit() const;

  /// Print a representation of the object
  virtual void repr(std::ostream &stream) const;

  /// Print a destription of the object
  virtual void print(std::ostream &stream) const;

  protected:
    /// Get a shared object from the current internal object
    template<class B>
    B shared_from_this();

    /// Get a shared object from the current internal object
    template<class B>
    const B shared_from_this() const;

    /// Has the function been initialized?
    bool is_init_;

  private:
    /// Number of references pointing to the object
    unsigned int count;
};

/// Typecast a shared object to a base class to a shared object to a derived class, cf. dynamic_cast
template<class B>
B shared_cast(SharedObject& A){
  
  /// Get a pointer to the node
  SharedObjectNode* ptr = A.get();
  
  /// Create a return object
  B ret;
  
  /// Assign node of B and return
  ret.assignNode(ptr);
      
  /// Null pointer if not pointing towards the right type of object
  if(ptr && !ret.checkNode()) 
    ret.assignNode(0);

  return ret;
}

/// Typecast a shared object to a base class to a shared object to a derived class, cf. dynamic_cast (const)
template<class B>
const B shared_cast(const SharedObject& A){
  SharedObject A_copy = A;
  return shared_cast<B>(A_copy);
}

/// Check if a shared object is of a certain type
template<class B>
bool is_a(const SharedObject& A){
  casadi_assert(!A.isNull());
  return !shared_cast<B>(A).isNull();
}

//@{
/// Make a deep copy of an object (Note: default is a shallow copy!)
template<class A>
A deepcopy(const A& a){
  A ret = a;
  ret.makeUnique();
  return ret;
}

template<class A>
A deepcopy(const A& a, std::map<SharedObjectNode*,SharedObject>& already_copied){
  A ret = a;
  ret.makeUnique(already_copied);
  return ret;
}

template<class A>
std::vector<A> deepcopy(const std::vector<A>& a, std::map<SharedObjectNode*,SharedObject>& already_copied){
  std::vector<A> ret = a;
  for(typename std::vector<A>::iterator it=ret.begin(); it!=ret.end(); ++it){
    it->makeUnique(already_copied);
  }
  return ret;
}
//@}

// Template function implementations
template<class B>
B SharedObjectNode::shared_from_this(){
  B ret;
  ret.assignNode(this);
  
  // Assert that the object is valid
  casadi_assert(ret.checkNode());
  
  return ret;
}

template<class B>
const B SharedObjectNode::shared_from_this() const{
  B ret;
  ret.assignNode(const_cast<SharedObjectNode*>(this));
  
  // Assert that the object is valid
  casadi_assert(ret.checkNode());
  
  return ret;
}


#endif // SWIG


} // namespace CasADi


#endif // SHARED_OBJECT_HPP

