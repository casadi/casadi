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
#include "serializer.hpp"

#include "function_internal.hpp"
#include "slice.hpp"
#include "linsol.hpp"
#include <iomanip>

using namespace std;
namespace casadi {

    DeSerializer::DeSerializer(std::istream& in_s) : in(in_s) {

    }

    Serializer::Serializer(std::ostream& out_s, const Dict& /*opts*/) : out(out_s) {

    }

    casadi_int Serializer::add(const Function& f) {
      // Quick return if already added
      for (auto&& e : added_functions_) if (e==f) return 0;

      added_functions_.push_back(f);

      f.serialize(*this);

      return 0;
    }

    void Serializer::decorate(char e) {
      pack(e);
    }

    void DeSerializer::assert_decoration(char e) {
      char t;
      unpack(t);
      casadi_assert(t==e, "Serializer error '" + str(e) + "' vs '" + str(t) + "'.");
    }

    void DeSerializer::unpack(casadi_int& e) {
      assert_decoration('J');
      int64_t n;
      char* c = reinterpret_cast<char*>(&n);

      for (int j=0;j<8;++j) unpack(c[j]);
      e = n;
    }

    void Serializer::pack(casadi_int e) {
      decorate('J');
      int64_t n = e;
      const char* c = reinterpret_cast<const char*>(&n);
      for (int j=0;j<8;++j) pack(c[j]);
    }

    void DeSerializer::unpack(int& e) {
      assert_decoration('i');
      int32_t n;
      char* c = reinterpret_cast<char*>(&n);

      for (int j=0;j<4;++j) unpack(c[j]);
      e = n;
    }

    void Serializer::pack(int e) {
      decorate('i');
      int32_t n = e;
      const char* c = reinterpret_cast<const char*>(&n);
      for (int j=0;j<4;++j) pack(c[j]);
    }

    void DeSerializer::unpack(bool& e) {
      assert_decoration('b');
      char n;
      unpack(n);
      e = n;
    }

    void Serializer::pack(bool e) {
      decorate('b');
      pack(static_cast<char>(e));
    }

    void DeSerializer::unpack(char& e) {
      unsigned char ref = 'a';
      in.get(e);
      char t;
      in.get(t);
      e = (reinterpret_cast<unsigned char&>(e)-ref) +
          ((reinterpret_cast<unsigned char&>(t)-ref) << 4);
    }

    void Serializer::pack(char e) {
      unsigned char ref = 'a';
      // Note: outputstreams work neatly with std::hex,
      // but inputstreams don't
      out.put(ref + (reinterpret_cast<unsigned char&>(e) % 16));
      out.put(ref + (reinterpret_cast<unsigned char&>(e) >> 4));
    }

    void Serializer::pack(const std::string& e) {
      decorate('s');
      int s = e.size();
      pack(s);
      const char* c = e.c_str();
      for (int j=0;j<s;++j) pack(c[j]);
    }

    void DeSerializer::unpack(std::string& e) {
      assert_decoration('s');
      int s;
      unpack(s);
      e.resize(s);
      for (int j=0;j<s;++j) unpack(e[j]);
    }

    void DeSerializer::unpack(double& e) {
      assert_decoration('d');
      char* c = reinterpret_cast<char*>(&e);
      for (int j=0;j<8;++j) unpack(c[j]);
    }

    void Serializer::pack(double e) {
      decorate('d');
      const char* c = reinterpret_cast<const char*>(&e);
      for (int j=0;j<8;++j) pack(c[j]);
    }

    void Serializer::pack(const Sparsity& e) {
      decorate('S');
      shared_pack(e, sparsities_);
    }

    void DeSerializer::unpack(Sparsity& e) {
      assert_decoration('S');
      shared_unpack(e, sparsities);
    }

    void Serializer::pack(const MX& e) {
      decorate('X');
      shared_pack(e, MX_nodes_);
    }

    void DeSerializer::unpack(MX& e) {
      assert_decoration('X');
      shared_unpack(e, nodes);
    }

    void Serializer::pack(const Function& e) {
      decorate('X');
      shared_pack(e, functions_);
    }

    void DeSerializer::unpack(Function& e) {
      assert_decoration('X');
      shared_unpack(e, functions);
    }

    void Serializer::pack(const Linsol& e) {
      decorate('L');
      shared_pack(e, linsols_);
    }

    void DeSerializer::unpack(Linsol& e) {
      assert_decoration('L');
      shared_unpack(e, linsols);
    }

    void Serializer::pack(const Slice& e) {
      decorate('S');
      e.serialize(*this);
    }

    void DeSerializer::unpack(Slice& e) {
      assert_decoration('S');
      e = Slice::deserialize(*this);
    }

    void Serializer::pack(const SXElem& e) {
      decorate('E');
      shared_pack(e, SX_nodes_);
    }

    void DeSerializer::unpack(SXElem& e) {
      assert_decoration('E');
      shared_unpack(e, sx_nodes);
    }

} // namespace casadi
