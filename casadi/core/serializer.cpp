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


#include "function.hpp"
#include "serializer.hpp"
#include "serializing_stream.hpp"
#include "slice.hpp"
#include "linsol.hpp"
#include "importer.hpp"
#include "generic_type.hpp"
#include <iomanip>

namespace casadi {

    StringSerializer::StringSerializer(const Dict& opts) :
        SerializerBase(std::unique_ptr<std::ostream>(new std::stringstream()), opts) {
    }

    FileSerializer::FileSerializer(const std::string& fname, const Dict& opts) :
        SerializerBase(
          std::unique_ptr<std::ostream>(
            new std::ofstream(fname, std::ios_base::binary | std::ios::out)),
          opts) {
      if ((sstream_->rdstate() & std::ifstream::failbit) != 0) {
        casadi_error("Could not open file '" + fname + "' for writing.");
      }
    }

    SerializerBase::SerializerBase(std::unique_ptr<std::ostream> stream, const Dict& opts) :
        sstream_(std::move(stream)),
        serializer_(new SerializingStream(*sstream_, opts)) {
    }

    std::string SerializerBase::type_to_string(SerializationType type) {
      switch (type) {
        case SERIALIZED_SPARSITY: return "sparsity";
        case SERIALIZED_MX: return "mx";
        case SERIALIZED_MX_v1: return "mx_v1";
        case SERIALIZED_DM: return "dm";
        case SERIALIZED_SX: return "sx";
        case SERIALIZED_SX_v1: return "sx_v1";
        case SERIALIZED_LINSOL: return "linsol";
        case SERIALIZED_FUNCTION: return "function";
        case SERIALIZED_GENERICTYPE: return "generictype";
        case SERIALIZED_INT: return "int";
        case SERIALIZED_DOUBLE: return "double";
        case SERIALIZED_STRING: return "string";
        case SERIALIZED_SPARSITY_VECTOR: return "sparsity_vector";
        case SERIALIZED_MX_VECTOR: return "mx_vector";
        case SERIALIZED_MX_VECTOR_v1: return "mx_vector_v1";
        case SERIALIZED_DM_VECTOR: return "dm_vector";
        case SERIALIZED_SX_VECTOR: return "sx_vector";
        case SERIALIZED_SX_VECTOR_v1: return "sx_vector_v1";
        case SERIALIZED_LINSOL_VECTOR: return "linsol_vector";
        case SERIALIZED_FUNCTION_VECTOR: return "function_vector";
        case SERIALIZED_GENERICTYPE_VECTOR: return "generictype_vector";
        case SERIALIZED_INT_VECTOR: return "int_vector";
        case SERIALIZED_DOUBLE_VECTOR: return "double_vector";
        case SERIALIZED_STRING_VECTOR: return "string_vector";
        default: casadi_error("Unknown type" + str(type));
      }
    }

    FileSerializer::~FileSerializer() {
    }

    std::string StringSerializer::encode() {
      std::string ret = static_cast<std::stringstream*>(sstream_.get())->str();
      static_cast<std::stringstream*>(sstream_.get())->str("");
      sstream_->clear();
      return ret;
    }
    void StringDeserializer::decode(const std::string& string) {
      casadi_assert(dstream_->peek() == std::char_traits<char>::eof(),
        "StringDeserializer::decode does not apply: current string not fully consumed yet.");
      static_cast<std::stringstream*>(dstream_.get())->str(string);
      dstream_->clear(); // reset error flags
    }

    SerializerBase::~SerializerBase() { }
    StringSerializer::~StringSerializer() { }

    DeserializerBase::DeserializerBase(std::unique_ptr<std::istream> stream) :
      dstream_(std::move(stream)),
      deserializer_(new DeserializingStream(*dstream_)) {
    }

    FileDeserializer::FileDeserializer(const std::string& fname) :
        DeserializerBase(std::unique_ptr<std::istream>(
          new std::ifstream(fname, std::ios_base::binary | std::ios::in))) {
      if ((dstream_->rdstate() & std::ifstream::failbit) != 0) {
        casadi_error("Could not open file '" + fname + "' for reading.");
      }
    }

    StringDeserializer::StringDeserializer(const std::string& string) :
        DeserializerBase(std::unique_ptr<std::istream>(
          new std::stringstream(string))) {
    }

    DeserializerBase::~DeserializerBase() { }
    StringDeserializer::~StringDeserializer() { }
    FileDeserializer::~FileDeserializer() { }

    SerializingStream& SerializerBase::serializer() {
      return *serializer_;
    }

    DeserializingStream& DeserializerBase::deserializer() {
      casadi_assert(dstream_->peek() != std::char_traits<char>::eof(),
        "Deserializer reached end of stream. Nothing left to unpack.");
      return *deserializer_;
    }

    SerializerBase::SerializationType DeserializerBase::pop_type() {
      char type;
      deserializer().unpack(type);
      return static_cast<SerializerBase::SerializationType>(type);
    }

    void SerializerBase::pack(const MX& e) {
      serializer().pack(static_cast<char>(SERIALIZED_MX));
      serializer().pack(Function::order({e}));
      serializer().pack(e);
    }
    void SerializerBase::pack(const std::vector<MX>& e) {
      serializer().pack(static_cast<char>(SERIALIZED_MX_VECTOR));
      serializer().pack(Function::order(e));
      serializer().pack(e);
    }
    void SerializerBase::pack(const SX& e) {
      serializer().pack(static_cast<char>(SERIALIZED_SX));
      serializer().pack(Function::order({e}));
      serializer().pack(e);
    }
    void SerializerBase::pack(const std::vector<SX>& e) {
      serializer().pack(static_cast<char>(SERIALIZED_SX_VECTOR));
      serializer().pack(Function::order(e));
      serializer().pack(e);
    }
    MX DeserializerBase::blind_unpack_mx() {
      std::vector<MX> sorted;
      deserializer().unpack(sorted);
      MX ret;
      deserializer().unpack(ret);
      return ret;
    }
    SX DeserializerBase::blind_unpack_sx() {
      std::vector<SX> sorted;
      deserializer().unpack(sorted);
      SX ret;
      deserializer().unpack(ret);
      return ret;
    }
    std::vector<MX> DeserializerBase::blind_unpack_mx_vector() {
      std::vector<MX> sorted;
      deserializer().unpack(sorted);
      std::vector<MX> ret;
      deserializer().unpack(ret);
      return ret;
    }
    std::vector<SX> DeserializerBase::blind_unpack_sx_vector() {
      std::vector<SX> sorted;
      deserializer().unpack(sorted);
      std::vector<SX> ret;
      deserializer().unpack(ret);
      return ret;
    }
    MX DeserializerBase::blind_unpack_mx_v1() {
      Function f;
      deserializer().unpack(f);
      MX ret;
      deserializer().unpack(ret);
      return ret;
    }
    SX DeserializerBase::blind_unpack_sx_v1() {
      Function f;
      deserializer().unpack(f);
      SX ret;
      deserializer().unpack(ret);
      return ret;
    }
    std::vector<MX> DeserializerBase::blind_unpack_mx_vector_v1() {
      Function f;
      deserializer().unpack(f);
      std::vector<MX> ret;
      deserializer().unpack(ret);
      return ret;
    }
    std::vector<SX> DeserializerBase::blind_unpack_sx_vector_v1() {
      Function f;
      deserializer().unpack(f);
      std::vector<SX> ret;
      deserializer().unpack(ret);
      return ret;
    }
    MX DeserializerBase::unpack_mx() {
      SerializerBase::SerializationType t = pop_type();
      if (t==SerializerBase::SerializationType::SERIALIZED_MX_v1) {
        return blind_unpack_mx_v1();
      }
      casadi_assert(t==SerializerBase::SerializationType::SERIALIZED_MX,
        "Expected to find a '" + SerializerBase::type_to_string(
          SerializerBase::SerializationType::SERIALIZED_MX)+
        "', but encountered a '" + SerializerBase::type_to_string(t) + "' instead.");
      return blind_unpack_mx();
    }
    SX DeserializerBase::unpack_sx() {
      SerializerBase::SerializationType t = pop_type();
      if (t==SerializerBase::SerializationType::SERIALIZED_SX_v1) {
        return blind_unpack_sx_v1();
      }
      casadi_assert(t==SerializerBase::SerializationType::SERIALIZED_SX,
        "Expected to find a '" + SerializerBase::type_to_string(
          SerializerBase::SerializationType::SERIALIZED_SX)+
        "', but encountered a '" + SerializerBase::type_to_string(t) + "' instead.");
      return blind_unpack_sx();
    }
    std::vector<MX> DeserializerBase::unpack_mx_vector() {
      SerializerBase::SerializationType t = pop_type();
      if (t==SerializerBase::SerializationType::SERIALIZED_MX_VECTOR_v1) {
        return blind_unpack_mx_vector_v1();
      }
      casadi_assert(t==SerializerBase::SerializationType::SERIALIZED_MX_VECTOR, \
        "Expected to find a '" + SerializerBase::type_to_string(
          SerializerBase::SerializationType::SERIALIZED_MX_VECTOR)+
        "', but encountered a '" + SerializerBase::type_to_string(t) + "' instead.");
      return blind_unpack_mx_vector();
    }
    std::vector<SX> DeserializerBase::unpack_sx_vector() {
      SerializerBase::SerializationType t = pop_type();
      if (t==SerializerBase::SerializationType::SERIALIZED_SX_VECTOR_v1) {
        return blind_unpack_sx_vector_v1();
      }
      casadi_assert(t==SerializerBase::SerializationType::SERIALIZED_SX_VECTOR, \
        "Expected to find a '" + SerializerBase::type_to_string(\
          SerializerBase::SerializationType::SERIALIZED_SX_VECTOR)+
        "', but encountered a '" + SerializerBase::type_to_string(t) + "' instead.");
      return blind_unpack_sx_vector();
    }

#define SERIALIZE(TYPE, Type, type) \
    void SerializerBase::pack(const Type& e) { \
      serializer().pack(static_cast<char>(SERIALIZED_ ## TYPE));\
      serializer().pack(e); \
    } \
    \
    Type DeserializerBase::blind_unpack_ ## type() { \
      Type ret;\
      deserializer().unpack(ret);\
      return ret;\
    }\
    Type DeserializerBase::unpack_ ## type() { \
      SerializerBase::SerializationType t = pop_type();\
      casadi_assert(t==SerializerBase::SerializationType::SERIALIZED_ ## TYPE, \
        "Expected to find a '" + SerializerBase::type_to_string(\
          SerializerBase::SerializationType::SERIALIZED_ ## TYPE) +\
        "', but encountered a '" + SerializerBase::type_to_string(t) + "' instead.");\
      return blind_unpack_ ## type();\
    }

#define SERIALIZE_ALL(TYPE, Type, type)\
  SERIALIZE(TYPE, Type, type)\
  SERIALIZE(TYPE ## _VECTOR, std::vector< Type >, type ## _vector)

SERIALIZE_ALL(SPARSITY, Sparsity, sparsity)
SERIALIZE_ALL(DM, DM, dm)
SERIALIZE_ALL(LINSOL, Linsol, linsol)
SERIALIZE_ALL(FUNCTION, Function, function)
SERIALIZE_ALL(GENERICTYPE, GenericType, generictype)
SERIALIZE_ALL(INT, casadi_int, int)
SERIALIZE_ALL(DOUBLE, double, double)
SERIALIZE_ALL(STRING, std::string , string)

  void SerializerBase::connect(DeserializerBase & s) {
    serializer_->connect(*s.deserializer_);
  }
  void SerializerBase::reset() {
    serializer_->reset();
  }
  void DeserializerBase::connect(SerializerBase & s) {
    deserializer_->connect(*s.serializer_);
  }
  void DeserializerBase::reset() {
    deserializer_->reset();
  }

} // namespace casadi
