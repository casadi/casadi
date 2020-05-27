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


#include "function.hpp"
#include "serializing_stream.hpp"
#include "slice.hpp"
#include "linsol.hpp"
#include "importer.hpp"
#include "generic_type.hpp"
#include "shared_object_internal.hpp"
#include "sx_node.hpp"
#include "sparsity_internal.hpp"
#include "mx_node.hpp"
#include "function_internal.hpp"
#include <iomanip>

using namespace std;
namespace casadi {

    static casadi_int serialization_protocol_version = 3;
    static casadi_int serialization_check = 123456789012345;

    DeserializingStream::DeserializingStream(std::istream& in_s) : in(in_s), debug_(false) {

      casadi_assert(in_s.good(), "Invalid input stream. If you specified an input file, "
        "make sure it exists relative to the current directory.");

      // Sanity check
      casadi_int check;
      unpack(check);
      casadi_assert(check==serialization_check,
        "DeserializingStream sanity check failed. "
        "Expected " + str(serialization_check) + ", but got " + str(check) + ".");

      // API version check
      casadi_int v;
      unpack(v);
      casadi_assert(v==serialization_protocol_version,
        "Serialization protocol is not compatible. "
        "Got version " + str(v) + ", while " +
        str(serialization_protocol_version) + " was expected.");

      bool debug;
      unpack(debug);
      debug_ = debug;

    }

    SerializingStream::SerializingStream(std::ostream& out_s) :
      SerializingStream(out_s, Dict()) {
    }

    SerializingStream::SerializingStream(std::ostream& out_s, const Dict& opts) :
        out(out_s), debug_(false) {
      // Sanity check
      pack(serialization_check);
      // API version check
      pack(casadi_int(serialization_protocol_version));

      bool debug = false;

      // Read options
      for (auto&& op : opts) {
        if (op.first=="debug") {
          debug = op.second;
        } else {
          casadi_error("Unknown option: '" + op.first + "'.");
        }
      }

      pack(debug);
      debug_ = debug;
    }

    void SerializingStream::decorate(char e) {
      if (debug_) pack(e);
    }

    void DeserializingStream::assert_decoration(char e) {
      if (debug_) {
        char t;
        unpack(t);
        casadi_assert(t==e, "DeserializingStream error '" + str(e) + "' vs '" + str(t) + "'.");
      }
    }

    void DeserializingStream::unpack(casadi_int& e) {
      assert_decoration('J');
      int64_t n;
      char* c = reinterpret_cast<char*>(&n);

      for (int j=0;j<8;++j) unpack(c[j]);
      e = n;
    }

    void SerializingStream::pack(casadi_int e) {
      decorate('J');
      int64_t n = e;
      const char* c = reinterpret_cast<const char*>(&n);
      for (int j=0;j<8;++j) pack(c[j]);
    }

    void SerializingStream::pack(size_t e) {
      decorate('K');
      uint64_t n = e;
      const char* c = reinterpret_cast<const char*>(&n);
      for (int j=0;j<8;++j) pack(c[j]);
    }

    void DeserializingStream::unpack(size_t& e) {
      assert_decoration('K');
      uint64_t n;
      char* c = reinterpret_cast<char*>(&n);

      for (int j=0;j<8;++j) unpack(c[j]);
      e = n;
    }

    void DeserializingStream::unpack(int& e) {
      assert_decoration('i');
      int32_t n;
      char* c = reinterpret_cast<char*>(&n);

      for (int j=0;j<4;++j) unpack(c[j]);
      e = n;
    }

    void SerializingStream::pack(int e) {
      decorate('i');
      int32_t n = e;
      const char* c = reinterpret_cast<const char*>(&n);
      for (int j=0;j<4;++j) pack(c[j]);
    }

    void DeserializingStream::unpack(bool& e) {
      assert_decoration('b');
      char n;
      unpack(n);
      e = n;
    }

    void SerializingStream::pack(bool e) {
      decorate('b');
      pack(static_cast<char>(e));
    }

    void DeserializingStream::unpack(char& e) {
      unsigned char ref = 'a';
      in.get(e);
      char t;
      in.get(t);
      e = (reinterpret_cast<unsigned char&>(e)-ref) +
          ((reinterpret_cast<unsigned char&>(t)-ref) << 4);
    }

    void SerializingStream::pack(char e) {
      unsigned char ref = 'a';
      // Note: outputstreams work neatly with std::hex,
      // but inputstreams don't
      out.put(ref + (reinterpret_cast<unsigned char&>(e) % 16));
      out.put(ref + (reinterpret_cast<unsigned char&>(e) >> 4));
    }

    void SerializingStream::pack(const std::string& e) {
      decorate('s');
      int s = e.size();
      pack(s);
      const char* c = e.c_str();
      for (int j=0;j<s;++j) pack(c[j]);
    }

    void DeserializingStream::unpack(std::string& e) {
      assert_decoration('s');
      int s;
      unpack(s);
      e.resize(s);
      for (int j=0;j<s;++j) unpack(e[j]);
    }

    void DeserializingStream::unpack(double& e) {
      assert_decoration('d');
      char* c = reinterpret_cast<char*>(&e);
      for (int j=0;j<8;++j) unpack(c[j]);
    }

    void SerializingStream::pack(double e) {
      decorate('d');
      const char* c = reinterpret_cast<const char*>(&e);
      for (int j=0;j<8;++j) pack(c[j]);
    }

    void SerializingStream::pack(const Sparsity& e) {
      decorate('S');
      shared_pack(e);
    }

    void DeserializingStream::unpack(Sparsity& e) {
      assert_decoration('S');
      shared_unpack<Sparsity, SparsityInternal>(e);
    }

    void SerializingStream::pack(const MX& e) {
      decorate('X');
      shared_pack(e);
    }

    void DeserializingStream::unpack(MX& e) {
      assert_decoration('X');
      shared_unpack<MX, MXNode>(e);
    }

    void SerializingStream::pack(const Function& e) {
      decorate('F');
      shared_pack(e);
    }

    void DeserializingStream::unpack(Function& e) {
      assert_decoration('F');
      shared_unpack<Function, FunctionInternal>(e);
    }

    void SerializingStream::pack(const Importer& e) {
      decorate('M');
      shared_pack(e);
    }

    void DeserializingStream::unpack(Importer& e) {
      assert_decoration('M');
      shared_unpack<Importer, ImporterInternal>(e);
    }

    void SerializingStream::pack(const Linsol& e) {
      decorate('L');
      shared_pack(e);
    }

    void DeserializingStream::unpack(Linsol& e) {
      assert_decoration('L');
      shared_unpack<Linsol, LinsolInternal>(e);
    }

    void SerializingStream::pack(const GenericType& e) {
      decorate('G');
      shared_pack(e);
    }

    void DeserializingStream::unpack(GenericType& e) {
      assert_decoration('G');
      shared_unpack<GenericType, SharedObjectInternal>(e);
    }

    void SerializingStream::pack(std::istream& s) {
      decorate('B');
      s.seekg(0, std::ios::end);
      size_t len = s.tellg();
      s.seekg(0, std::ios::beg);
      pack(len);
      char buffer[1024];
      for (size_t i=0;i<len;++i) {
        s.read(buffer, 1024);
        size_t c = s.gcount();
        for (size_t j=0;j<c;++j) {
          pack(buffer[j]);
        }
        if (s.rdstate() & std::ifstream::eofbit) break;
      }
    }

    void DeserializingStream::unpack(std::ostream& s) {
      assert_decoration('B');
      size_t len;
      unpack(len);
      for (size_t i=0;i<len;++i) {
        char c;
        unpack(c);
        s.put(c);
      }
    }

    void SerializingStream::pack(const Slice& e) {
      decorate('S');
      e.serialize(*this);
    }

    void DeserializingStream::unpack(Slice& e) {
      assert_decoration('S');
      e = Slice::deserialize(*this);
    }

    void SerializingStream::pack(const SXElem& e) {
      decorate('E');
      shared_pack(e);
    }

    void DeserializingStream::unpack(SXElem& e) {
      assert_decoration('E');
      shared_unpack<SXElem, SXNode>(e);
    }

  template<>
  void DeserializingStream::unpack(std::vector<bool>& e) {
    assert_decoration('V');
    casadi_int s;
    unpack(s);
    e.resize(s);
    for (casadi_int i=0;i<s;++i) {
      bool b;
      unpack(b);
      e[i] = b;
    }
  }

  int DeserializingStream::version(const std::string& name) {
    int load_version;
    unpack(name+"::serialization::version", load_version);
    return load_version;
  }

  int DeserializingStream::version(const std::string& name, int min, int max) {
    int load_version = version(name);
    casadi_assert(load_version>=min && load_version<=max,
      "DeSerialization of " + name + " failed. "
      "Object written in version " + str(load_version) +
      " but can only read version " + str(min) + "..." + str(max) + ".");
    return load_version;
  }

  void DeserializingStream::version(const std::string& name, int v) {
    int load_version = version(name);
    casadi_assert(load_version==v,
      "DeSerialization of " + name + " failed. "
      "Object written in version " + str(load_version) +
      " but can only read in version " + str(v) + ".");
  }

  void SerializingStream::version(const std::string& name, int v) {
    pack(name+"::serialization::version", v);
  }

  UniversalNodeOwner::UniversalNodeOwner(SharedObjectInternal* obj) :
      node(obj), is_sx(false) {
    if (node) obj->count++;
  }

  UniversalNodeOwner::UniversalNodeOwner(SXNode* obj) :
      node(obj), is_sx(true) {
    if (node) obj->count++;
  }

  UniversalNodeOwner::UniversalNodeOwner(UniversalNodeOwner&& rhs) noexcept :
    node(std::move(rhs.node)), is_sx(std::move(rhs.is_sx)) {
    rhs.node = nullptr;
  }

  UniversalNodeOwner& UniversalNodeOwner::operator=(UniversalNodeOwner&& other) noexcept {
    std::swap(node, other.node);
    std::swap(is_sx, other.is_sx);
    return *this;
  }

  UniversalNodeOwner::~UniversalNodeOwner() {
    if (!node) return;
    if (is_sx) {
      if (--static_cast<SXNode*>(node)->count == 0) {
        delete static_cast<SXNode*>(node);
      }
    } else {
      if (--static_cast<SharedObjectInternal*>(node)->count == 0) {
        delete static_cast<SharedObjectInternal*>(node);
      }
    }
  }

  void SerializingStream::connect(DeserializingStream & s) {
    nodes_ = &s.nodes_;
  }

  void DeserializingStream::connect(SerializingStream & s) {
    shared_map_ = &s.shared_map_;
  }

  void SerializingStream::reset() {
    shared_map_.clear();
  }

  void DeserializingStream::reset() {
    nodes_.clear();
  }

} // namespace casadi
