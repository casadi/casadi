/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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

#include "casadi_common.hpp"
#include "exception.hpp"
#include <map>
#include <vector>
#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.mutex.h>
#else // CASADI_WITH_THREAD_MINGW
#include <mutex>
#endif // CASADI_WITH_THREAD_MINGW
#endif //CASADI_WITH_THREAD

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
#include <memory>
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

namespace casadi {

  // Forward declaration of weak reference class
  class WeakRef;

  /// \cond INTERNAL
  // Forward declaration of internal classes
  class SharedObjectInternal;
  class WeakRefInternal;
  /// \endcond

  /** \brief SharedObject implements a reference counting framework similar for efficient and

      easily-maintained memory management.

      To use the class, both the SharedObject class (the public class), and the SharedObjectInternal
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

      \identifier{as} */
  class CASADI_EXPORT SharedObject {
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
    void own(SharedObjectInternal* node);

    /** \brief Assign the node to a node class pointer without reference counting
     *
     * improper use will cause memory leaks!

        \identifier{at} */
    void assign(SharedObjectInternal* node);

    /// Get a const pointer to the node
    SharedObjectInternal* get() const;

    /// Get the reference count
    casadi_int getCount() const;

    /// Swap content with another instance
    void swap(SharedObject& other);

    /// Access a member function or object
    SharedObjectInternal* operator->() const;
    /// \endcond
#endif // SWIG

    /** \brief Get class name

        \identifier{au} */
    std::string class_name() const;

    /// Print a description of the object
    void disp(std::ostream& stream, bool more=false) const;

    /// Get string representation
    std::string get_str(bool more=false) const {
      std::stringstream ss;
      disp(ss, more);
      return ss.str();
    }

    /// \cond INTERNAL
    /// Print the pointer to the internal class
    void print_ptr(std::ostream &stream=casadi::uout()) const;
    /// \endcond

    /// Is a null pointer?
    bool is_null() const;

    /** \brief Returns a number that is unique for a given Node.

     * If the Object does not point to any node, "0" is returned.

        \identifier{av} */
    casadi_int __hash__() const;

/// \cond INTERNAL
#ifndef SWIG
    /** \brief Get a weak reference to the object

        \identifier{aw} */
    WeakRef* weak();
  protected:
    void count_up(); // increase counter of the node
    void count_down(); // decrease counter of the node
  private:
    SharedObjectInternal *node;
#endif // SWIG
/// \endcond
  };

/** \brief Weak reference type

    A weak reference to a SharedObject
    \author Joel Andersson
    \date 2013

    \identifier{ax} */
  class CASADI_EXPORT WeakRef : public SharedObject {
  public:
    friend class SharedObjectInternal;

    /** \brief Default constructor

        \identifier{ay} */
    WeakRef(int dummy=0);

    /** \brief Construct from a shared object (also implicit type conversion)

        \identifier{az} */
    WeakRef(SharedObject shared);

    /** \brief Get a shared (owning) reference

        \identifier{b0} */
    SharedObject shared() const;

    /** \brief Check if alive

        \identifier{b1} */
    bool alive() const;

    /** \brief Thread-safe alternative to alive()/shared() */
    bool shared_if_alive(SharedObject& shared) const;

    /** \brief  Access functions of the node

        \identifier{b2} */
    WeakRefInternal* operator->();

    /** \brief  Const access functions of the node

        \identifier{b3} */
    const WeakRefInternal* operator->() const;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    std::shared_ptr<std::mutex> get_mutex() const;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

#ifndef SWIG
  private:
    /** \brief Construct from a shared object (internal)

        \identifier{b4} */
    explicit WeakRef(SharedObjectInternal* raw);

    /** \brief The shared object has been deleted

        \identifier{b5} */
    void kill();
#endif // SWIG
  };

#ifndef SWIG

  /** \brief Typecast a shared object to a base class to a shared object to a derived class,

   * cf. dynamic_cast

      \identifier{b6} */
  template<class B>
  B shared_cast(SharedObject& A) {

    /// Get a pointer to the node
    SharedObjectInternal* ptr = A.get();

    /// Create a return object
    B ret;

    /// Quick return if not allowed
    if (!B::test_cast(ptr)) return ret;

    /// Assign node of B and return
    ret.own(ptr);
    return ret;
  }

  /** \brief Typecast a shared object to a base class to a shared object to a derived class,

   * cf. dynamic_cast (const)

      \identifier{b7} */
  template<class B>
  const B shared_cast(const SharedObject& A) {
    SharedObject A_copy = A;
    return shared_cast<B>(A_copy);
  }

#endif // SWIG


template<typename K, typename T>
class CASADI_EXPORT WeakCache {
  public:
    void tocache(const K& key, const T& f, bool needs_lock=true) {
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
      // Safe access to cache_
      casadi::conditional_lock_guard<std::mutex> lock(mtx_, needs_lock);
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
      // Add to cache
      cache_.insert(std::make_pair(key, f));
      // Remove a lost reference, if any, to prevent uncontrolled growth
      for (auto it = cache_.begin(); it!=cache_.end(); ++it) {
        if (!it->second.alive()) {
          cache_.erase(it);
          break; // just one dead reference is enough
        }
      }
    }
    /* \brief Thread-safe unique caching
    * While an incache/tocache pair in multi-threaded context is safe
    * it may lead to fresh cache entries being overwritten.
    *
    * A mutex lock_guard on the scope of an incache/tocache pair
    * may lead to deadlocks.
    *
    */
    void tocache_if_missing(const K& key, T& f) {
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
      // Safe access to cache_
      std::lock_guard<std::mutex> lock(mtx_);
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
      if (!incache(key, f, false)) {
        tocache(key, f, false);
      }
    }
    bool incache(const K& key, T& f, bool needs_lock=true) const {
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
      // Safe access to cache_
      casadi::conditional_lock_guard<std::mutex> lock(mtx_, needs_lock);
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
      auto it = cache_.find(key);
      SharedObject temp;
      if (it!=cache_.end() && it->second.shared_if_alive(temp)) {
        f = shared_cast<T>(temp);
        return true;
      } else {
        return false;
      }
    }
    void cache(std::vector<K>& keys, std::vector<T>& entries) const {
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
      // Safe access to cache_
      std::lock_guard<std::mutex> lock(mtx_);
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
      keys.clear();
      entries.clear();
      // Add all entries that haven't been deleted
      for (auto&& cf : cache_) {
        SharedObject temp;
        if (cf.second.shared_if_alive(temp)) {
          keys.push_back(cf.first);
          entries.push_back(shared_cast<T>(temp));
        }
      }
    }
  private:
    std::map<K, WeakRef> cache_;
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    mutable std::mutex mtx_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
};

} // namespace casadi


#endif // CASADI_SHARED_OBJECT_HPP
