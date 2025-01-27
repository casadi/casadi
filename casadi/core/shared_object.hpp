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

#include "generic_shared_internal.hpp"

namespace casadi {

  /// \cond INTERNAL
  // Forward declaration of internal classes
  class SharedObjectInternal;
  class WeakRefInternal;
  /// \endcond

  /** \brief GenericShared implements a reference counting framework similar for efficient and 

      easily-maintained memory management.

      To use the class, both the GenericShared class (the public class), and the GenericSharedInternal
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
  class CASADI_EXPORT SharedObject :
      public GenericShared<SharedObject, SharedObjectInternal> {
    public:

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

      using internal_base_type = SharedObjectInternal;
      using base_type = SharedObject;

  };

/** \brief Weak reference type

    A weak reference to a GenericShared
    \author Joel Andersson
    \date 2013

    \identifier{ax} */
  class CASADI_EXPORT WeakRef :
      public GenericWeakRef<SharedObject, SharedObjectInternal> {
  public:
    WeakRef(int dummy=0) : GenericWeakRef<SharedObject, SharedObjectInternal>(dummy) {
    }
    WeakRef(SharedObject shared) : GenericWeakRef<SharedObject, SharedObjectInternal>(shared) {
    }
  /*private:
    explicit WeakRef(SharedObjectInternal* raw) : GenericWeakRef<SharedObject, SharedObjectInternal>(raw) {
    };*/
  };

#ifndef SWIG
  class CASADI_EXPORT SharedObjectInternal :
    public GenericSharedInternal<SharedObject, SharedObjectInternal> {
    friend class GenericShared<SharedObject, SharedObjectInternal>;
    friend class SharedObject;
    friend class GenericWeakRef<SharedObject, SharedObjectInternal>;
    friend class GenericSharedInternal<SharedObject, SharedObjectInternal>;
    friend class Memory;
    friend class UniversalNodeOwner;
  public:
    /// NOTE: these two constructors added because defaults are ill-formed defaults
    /// This may hint at a bug
    /// Default constructor
    SharedObjectInternal() : GenericSharedInternal<SharedObject, SharedObjectInternal>() {
    }
    /// Copy constructor
    SharedObjectInternal(const SharedObjectInternal& node) :
      GenericSharedInternal<SharedObject, SharedObjectInternal>(node) {
    }

    /// Readable name of the internal class
    virtual std::string class_name() const = 0;

    /// Print a description of the object
    virtual void disp(std::ostream& stream, bool more) const = 0;

    using weak_ref_type = WeakRefInternal;

  private:
    /// Number of references pointing to the object
#ifdef CASADI_WITH_THREAD
    std::atomic<casadi_int> count;
#else // CASADI_WITH_THREAD
    casadi_int count;
#endif

  };

  class CASADI_EXPORT WeakRefInternal :
      public GenericWeakRefInternal<SharedObject, SharedObjectInternal> {
    public:
      /// Print a description of the object
    void disp(std::ostream& stream, bool more) const override;

    /// Readable name of the class
    std::string class_name() const override {return "WeakRefInternal";}

    using GenericWeakRefInternal<SharedObject, SharedObjectInternal>::GenericWeakRefInternal;

  };


#endif // SWIG

} // namespace casadi

#endif // CASADI_SHARED_OBJECT_HPP
