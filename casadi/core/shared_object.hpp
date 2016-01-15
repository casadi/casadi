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


#ifndef CASADI_SHARED_OBJECT_HPP
#define CASADI_SHARED_OBJECT_HPP

#include "printable_object.hpp"
#include "exception.hpp"
#include <map>
#include <vector>

namespace casadi {

  // Forward declaration of weak reference class
  class WeakRef;

  /// \cond INTERNAL
  // Forward declaration of internal class
  class SharedObjectNode;
  /// \endcond

  /** \brief SharedObject implements a reference counting framework similar for efficient and
      easily-maintained memory management.

      To use the class, both the SharedObject class (the public class), and the SharedObjectNode
      class (the internal class) must be inherited from. It can be done in two different files
      and together with memory management, this approach provides a clear distinction of which
      methods of the class are to be considered "public", i.e. methods for public use that can
      be considered to remain over time with small changes, and the internal memory.

      When interfacing a software, which typically includes including some header file,
      this is best done only in the file where the internal class is defined, to avoid polluting
      the global namespace and other side effects.

      The default constructor always means creating a null pointer to an internal class only.
      To allocate an internal class (this works only when the internal class isn't abstract),
      use the constructor with arguments.

      The copy constructor and the assignment operator perform shallow copies only,
      to make a deep copy you must use the clone method explicitly.
      This will give a shared pointer instance.

      In an inheritance hierarchy, you can cast down automatically,
      e.g. (SXFunction is a child class of Function):
      SXFunction derived(...);
      Function base = derived;

      To cast up, use the shared_cast template function, which works analogously to
      dynamic_cast, static_cast, const_cast etc, e.g.:
      SXFunction derived(...);
      Function base = derived;
      SXFunction derived_from_base = shared_cast<SXFunction>(base);

      A failed shared_cast will result in a null pointer (cf. dynamic_cast)

      \author Joel Andersson
      \date 2010
  */
  class CASADI_EXPORT SharedObject : public PrintableObject<SharedObject> {
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

    /// Destructor
    ~SharedObject();

    /// Assignment operator
    SharedObject& operator=(const SharedObject& ref);

    /// \cond INTERNAL
    /// Assign the node to a node class pointer (or null)
    void assignNode(SharedObjectNode* node);

    /** \brief Assign the node to a node class pointer without reference counting
     *
     * improper use will cause memory leaks!
     */
    void assignNodeNoCount(SharedObjectNode* node);

    /// Get a const pointer to the node
    SharedObjectNode* get() const;

    /// Get the reference count
    int getCount() const;

    /// Swap content with another instance
    void swap(SharedObject& other);

    /// Access a member function or object
    SharedObjectNode* operator->() const;
    /// \endcond
#endif // SWIG

    /// Print a representation of the object
    void repr(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// Print a description of the object
    void print(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// \cond INTERNAL
    /// Print the pointer to the internal class
    void printPtr(std::ostream &stream=casadi::userOut()) const;
    /// \endcond

    /// Is a null pointer?
    bool is_null() const;

/// \cond INTERNAL
#ifndef SWIG
    /** \brief Get a weak reference to the object */
    WeakRef* weak();
#endif // SWIG
/// \endcond

    /** \brief Returns a number that is unique for a given Node.
     * If the Object does not point to any node, "0" is returned.
     */
    size_t __hash__() const;

  protected:
    void count_up(); // increase counter of the node
    void count_down(); // decrease counter of the node
  private:
    SharedObjectNode *node;

  };

#ifndef SWIG
  /// \cond INTERNAL
  /// Internal class for the reference counting framework, see comments on the public class.
  class CASADI_EXPORT SharedObjectNode {
    friend class SharedObject;
    friend class Memory;
  public:

    /// Default constructor
    SharedObjectNode();

    /// Copy constructor
    SharedObjectNode(const SharedObjectNode& node);

    /// Assignment operator
    SharedObjectNode& operator=(const SharedObjectNode& node);

    /// Destructor
    virtual ~SharedObjectNode() = 0;

    /// Get the reference count
    int getCount() const;

    /// Print a representation of the object
    virtual void repr(std::ostream &stream) const;

    /// Print a description of the object
    virtual void print(std::ostream &stream) const;

    /** \brief Get a weak reference to the object */
    WeakRef* weak();

  protected:
    /** Called in the constructor of singletons to avoid that the counter reaches zero */
    void initSingleton() {
      casadi_assert(count==0);
      count++;
    }

    /** Called in the destructor of singletons */
    void destroySingleton() {
      count--;
    }

    /// Get a shared object from the current internal object
    template<class B>
    B shared_from_this();

    /// Get a shared object from the current internal object
    template<class B>
    const B shared_from_this() const;

  private:
    /// Number of references pointing to the object
    unsigned int count;

    /// Weak pointer (non-owning) object for the object
    WeakRef* weak_ref_;
  };
  /// \endcond

  /// \cond INTERNAL
  /** \brief Typecast a shared object to a base class to a shared object to a derived class,
   * cf. dynamic_cast
   */
  template<class B>
  B shared_cast(SharedObject& A) {

    /// Get a pointer to the node
    SharedObjectNode* ptr = A.get();

    /// Create a return object
    B ret;

    /// Quick return if not allowed
    if (!B::test_cast(ptr)) return ret;

    /// Assign node of B and return
    ret.assignNode(ptr);
    return ret;
  }

  /** \brief Typecast a shared object to a base class to a shared object to a derived class,
   * cf. dynamic_cast (const)
   */
  template<class B>
  const B shared_cast(const SharedObject& A) {
    SharedObject A_copy = A;
    return shared_cast<B>(A_copy);
  }
  /// \endcond

  ///@{
  /// \cond INTERNAL
  template<class A>
  A getcopy(const A& a, std::map<SharedObjectNode*, SharedObject>& already_copied) {
    A ret;
    if (!a.is_null()) {
      std::map<SharedObjectNode*, SharedObject>::iterator it =
          already_copied.find(const_cast<SharedObjectNode*>(a.get()));
      if (it!=already_copied.end()) {
        ret.assignNode(it->second.get());
      }
    }
    return ret;
  }
  /// \endcond
  ///@}

  /// \cond INTERNAL
  /// Template function implementations
  template<class B>
  B SharedObjectNode::shared_from_this() {
    casadi_assert(B::test_cast(this));
    B ret;
    ret.assignNode(this);
    return ret;
  }

  template<class B>
  const B SharedObjectNode::shared_from_this() const {
    casadi_assert(B::test_cast(this));
    B ret;
    ret.assignNode(const_cast<SharedObjectNode*>(this));
    return ret;
  }
  /// \endcond

#endif // SWIG


} // namespace casadi


#endif // CASADI_SHARED_OBJECT_HPP

