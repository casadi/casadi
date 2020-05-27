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


#ifndef CASADI_SERIALIZING_STREAM_HPP
#define CASADI_SERIALIZING_STREAM_HPP

#include <set>
#include <sstream>
#include <unordered_map>

namespace casadi {
  class Slice;
  class Linsol;
  class Sparsity;
  class Function;
  class MX;
  class SXElem;
  class GenericType;
  class Importer;
  class SharedObject;
  class SharedObjectInternal;
  class SXNode;
  class SerializingStream;
  class UniversalNodeOwner {
  public:
    UniversalNodeOwner() = delete;
    UniversalNodeOwner(const UniversalNodeOwner&) = delete;
    UniversalNodeOwner(UniversalNodeOwner&& rhs) noexcept;
    UniversalNodeOwner(SharedObjectInternal* obj);
    UniversalNodeOwner(SXNode* obj);
    UniversalNodeOwner& operator=(const UniversalNodeOwner& other) = delete;
    UniversalNodeOwner& operator=(UniversalNodeOwner&& other) noexcept;
    ~UniversalNodeOwner();
    void* get() { return node; }
  private:
    void* node;
    bool is_sx;
  };
  typedef std::map<std::string, GenericType> Dict;

  /** \brief Helper class for Serialization
      \author Joris Gillis
      \date 2018
  */
  class CASADI_EXPORT DeserializingStream {
    friend class SerializingStream;
  public:
    /// Constructor
    DeserializingStream(std::istream &in_s);
    DeserializingStream(const DeserializingStream&) = delete;

    //@{
    /** \brief Reconstruct an object from the input stream
    *
    * If the reference is not of the same type as the object encoded in the stream.
    * an error will be raised.
    */
    void unpack(Sparsity& e);
    void unpack(MX& e);
    void unpack(SXElem& e);
    void unpack(Linsol& e);
    template <class T>
    void unpack(Matrix<T>& e) {
      e = Matrix<T>::deserialize(*this);
    }
    void unpack(Function& e);
    void unpack(Importer& e);
    void unpack(GenericType& e);
    void unpack(std::ostream& s);
    void unpack(Slice& e);
    void unpack(int& e);
    void unpack(bool& e);
    void unpack(casadi_int& e);
    void unpack(size_t& e);
    void unpack(std::string& e);
    void unpack(double& e);
    void unpack(char& e);
    template <class T>
    void unpack(std::vector<T>& e) {
      assert_decoration('V');
      casadi_int s;
      unpack(s);
      e.resize(s);
      for (T& i : e) unpack(i);
    }

    template <class K, class V>
    void unpack(std::map<K, V>& e) {
      assert_decoration('D');
      casadi_int s;
      unpack(s);
      e.clear();
      for (casadi_int i=0;i<s;++i) {
        K k;
        V v;
        unpack(k);
        unpack(v);
        e[k] = v;
      }
    }

    template <class A, class B>
    void unpack(std::pair<A, B>& e) {
      assert_decoration('p');
      unpack(e.first);
      unpack(e.second);
    }

    template <class T>
    void unpack(const std::string& descr, T& e) {
      if (debug_) {
        std::string d;
        unpack(d);
        casadi_assert(d==descr, "Mismatch: '" + descr + "' expected, got '" + d + "'.");
      }
      unpack(e);
    }
    //@}

    void version(const std::string& name, int v);
    int version(const std::string& name);
    int version(const std::string& name, int min, int max);

    void connect(SerializingStream & s);
    void reset();

  private:

    /* \brief Unpacks a shared object
    * 
    * Also treats SXNode, which is not actually a SharedObjectInternal
    */
    template <class T, class M>
    void shared_unpack(T& e) {
      char i;
      unpack("Shared::flag", i);
      switch (i) {
        case 'd': // definition
          e = T::deserialize(*this);
          if (shared_map_) (*shared_map_)[e.get()] = nodes_.size();
          nodes_.emplace_back(e.get());
          break;
        case 'r': // reference
          {
            casadi_int k;
            unpack("Shared::reference", k);
            UniversalNodeOwner& t = nodes_.at(k);
            e = T::create(static_cast<M*>(t.get()));
          }
          break;
        default:
          casadi_assert_dev(false);
      }
    }

    /** \brief Primitive typecheck during deserialization
     *
     * No-op unless in debug mode
     */
    void assert_decoration(char e);

    /// Collection of all shared pointer deserialized so far
    std::vector<UniversalNodeOwner> nodes_;
    std::unordered_map<void*, casadi_int>* shared_map_ = nullptr;
    /// Input stream
    std::istream& in;
    /// Debug mode?
    bool debug_;
  };

  /** \brief Helper class for Serialization


      \author Joris Gillis
      \date 2018
  */
  class CASADI_EXPORT SerializingStream {
    friend class DeserializingStream;
  public:
    /// Constructor
    SerializingStream(std::ostream& out);
    SerializingStream(std::ostream& out, const Dict& opts);

    // @{
    /** \brief Serializes an object to the output stream  */
    void pack(const Sparsity& e);
    void pack(const MX& e);
    void pack(const SXElem& e);
    void pack(const Linsol& e);
    template <class T>
    void pack(const Matrix<T>& e) {
      e.serialize(*this);
    }
    void pack(const Function& e);
    void pack(const Importer& e);
    void pack(const Slice& e);
    void pack(const GenericType& e);
    void pack(std::istream& s);
    void pack(int e);
    void pack(bool e);
    void pack(casadi_int e);
    void pack(size_t e);
    void pack(double e);
    void pack(const std::string& e);
    void pack(char e);
    template <class T>
    void pack(const std::vector<T>& e) {
      decorate('V');
      pack(static_cast<casadi_int>(e.size()));
      for (const T & i : e) pack(i);
    }
    template <class K, class V>
    void pack(const std::map<K, V>& e) {
      decorate('D');
      pack(static_cast<casadi_int>(e.size()));
      for (const auto & i : e) {
        pack(i.first);
        pack(i.second);
      }
    }
    template <class A, class B>
    void pack(const std::pair<A, B>& e) {
      decorate('p');
      pack(e.first);
      pack(e.second);
    }
    template <class T>
    void pack(const std::string& descr, const T& e) {
      if (debug_) pack(descr);
      pack(e);
    }
    template <class T>
    void pack(const std::string& descr, T& e) {
      if (debug_) pack(descr);
      pack(e);
    }
    //@}

    void version(const std::string& name, int v);

    void connect(DeserializingStream & s);
    void reset();

  private:
    /** \brief Insert information for a primitive typecheck during deserialization
     *
     * No-op unless in debug mode
     */
    void decorate(char e);

    /* \brief Packs a shared object
    * 
    * Also treats SXNode, which is not actually a SharedObjectInternal
    */
    template <class T>
    void shared_pack(const T& e) {
      auto it = shared_map_.find(e.get());
      if (it==shared_map_.end()) {
        // Not found
        pack("Shared::flag", 'd'); // definition
        e.serialize(*this);
        casadi_int r = shared_map_.size();
        shared_map_[e.get()] = r;
        if (nodes_) nodes_->emplace_back(e.get());
      } else {
        pack("Shared::flag", 'r'); // reference
        pack("Shared::reference", it->second);
      }
    }

    /// Mapping from shared pointers to running counter
    std::unordered_map<void*, casadi_int> shared_map_;
    std::vector<UniversalNodeOwner>* nodes_ = nullptr;
    /// Output stream
    std::ostream& out;
    /// Debug mode?
    bool debug_;
  };

  template <>
  CASADI_EXPORT void DeserializingStream::unpack(std::vector<bool>& e);

} // namespace casadi

#endif // CASADI_SERIALIZING_STREAM_HPP
