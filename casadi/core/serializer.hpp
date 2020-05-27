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


#ifndef CASADI_SERIALIZER_HPP
#define CASADI_SERIALIZER_HPP

#include <memory>
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
  class DeserializerBase;



  class CASADI_EXPORT SerializerBase {
    friend class DeserializerBase;
  public:
#ifndef SWIG
    SerializerBase(std::unique_ptr<std::ostream> stream, const Dict& opts = Dict());
#endif // SWIG
    ~SerializerBase();
    void pack(const Sparsity& e);
    void pack(const MX& e);
    void pack(const Matrix<double>& e);
    void pack(const Matrix<SXElem>& e);
    void pack(const Linsol& e);
    void pack(const Function& e);
    void pack(const GenericType& e);
    void pack(const casadi_int& e);
    void pack(const double& e);
    void pack(const std::string& e);
    void pack(const std::vector<Sparsity>& e);
    void pack(const std::vector<MX>& e);
    void pack(const std::vector< Matrix<double> >& e);
    void pack(const std::vector< Matrix<SXElem> >& e);
    void pack(const std::vector<Linsol>& e);
    void pack(const std::vector<Function>& e);
    void pack(const std::vector<GenericType>& e);
    void pack(const std::vector<casadi_int>& e);
    void pack(const std::vector<double>& e);
    void pack(const std::vector<std::string>& e);



    enum SerializationType {
      SERIALIZED_SPARSITY,
      SERIALIZED_MX,
      SERIALIZED_DM,
      SERIALIZED_SX,
      SERIALIZED_LINSOL,
      SERIALIZED_FUNCTION,
      SERIALIZED_GENERICTYPE,
      SERIALIZED_INT,
      SERIALIZED_DOUBLE,
      SERIALIZED_STRING,
      SERIALIZED_SPARSITY_VECTOR,
      SERIALIZED_MX_VECTOR,
      SERIALIZED_DM_VECTOR,
      SERIALIZED_SX_VECTOR,
      SERIALIZED_LINSOL_VECTOR,
      SERIALIZED_FUNCTION_VECTOR,
      SERIALIZED_GENERICTYPE_VECTOR,
      SERIALIZED_INT_VECTOR,
      SERIALIZED_DOUBLE_VECTOR,
      SERIALIZED_STRING_VECTOR,
    };

    static std::string type_to_string(SerializationType type);

    void connect(DeserializerBase & s);
    void reset();

  protected:
    SerializingStream& serializer();
    std::unique_ptr<std::ostream> sstream_;
    std::unique_ptr<SerializingStream> serializer_;
  };

  class CASADI_EXPORT DeserializerBase {
    friend class SerializerBase;
  public:
#ifndef SWIG
    DeserializerBase(std::unique_ptr<std::istream> stream);
#endif // SWIG
    ~DeserializerBase();

    SerializerBase::SerializationType pop_type();

    Sparsity blind_unpack_sparsity();
    MX blind_unpack_mx();
    Matrix<double> blind_unpack_dm();
    Matrix<SXElem> blind_unpack_sx();
    Linsol blind_unpack_linsol();
    Function blind_unpack_function();
    GenericType blind_unpack_generictype();
    casadi_int blind_unpack_int();
    double blind_unpack_double();
    std::string blind_unpack_string();
    std::vector<Sparsity> blind_unpack_sparsity_vector();
    std::vector<MX> blind_unpack_mx_vector();
    std::vector< Matrix<double> > blind_unpack_dm_vector();
    std::vector< Matrix<SXElem> > blind_unpack_sx_vector();
    std::vector<Linsol> blind_unpack_linsol_vector();
    std::vector<Function> blind_unpack_function_vector();
    std::vector<GenericType> blind_unpack_generictype_vector();
    std::vector<casadi_int> blind_unpack_int_vector();
    std::vector<double> blind_unpack_double_vector();
    std::vector<std::string> blind_unpack_string_vector();

    Sparsity unpack_sparsity();
    MX unpack_mx();
    Matrix<double> unpack_dm();
    Matrix<SXElem> unpack_sx();
    Linsol unpack_linsol();
    Function unpack_function();
    GenericType unpack_generictype();
    casadi_int unpack_int();
    double unpack_double();
    std::string unpack_string();
    std::vector<Sparsity> unpack_sparsity_vector();
    std::vector<MX> unpack_mx_vector();
    std::vector< Matrix<double> > unpack_dm_vector();
    std::vector< Matrix<SXElem> > unpack_sx_vector();
    std::vector<Linsol> unpack_linsol_vector();
    std::vector<Function> unpack_function_vector();
    std::vector<GenericType> unpack_generictype_vector();
    std::vector<casadi_int> unpack_int_vector();
    std::vector<double> unpack_double_vector();
    std::vector<std::string> unpack_string_vector();

    void connect(SerializerBase & s);
    void reset();

  protected:
    DeserializingStream& deserializer();
    std::unique_ptr<std::istream> dstream_;
    std::unique_ptr<DeserializingStream> deserializer_;
  };

  class CASADI_EXPORT StringSerializer : public SerializerBase {
  public:
    /** \brief Advanced serialization of CasADi objects
     * 
     * This class is intended for advanced users that want to circumvent the restrictions
     * of standard pickling/matlab save load, ie no raw SX/MX symbols allowed.
     * 
     * \example
     * x = SX.sym('x');
     * s = StringSerializer();
     * s.pack(x);
     * s.pack(sin(x));
     * 
     * data = s.encode();
     * 
     * s = StringDeserializer(data);
     * a = s.unpack();
     * b = s.unpack();
     * \endexample
     * 
     * Note:
     *  Saving SX/MX objects individually has a substantial overhead
     *  (both time and length of encoded string).
     *  You are encouraged to use the vector/list variants of 'save' for SX/MX to reduce
     *  the overhead.
     * 
     * 
     * \seealso Function::save, Function::serialize, StringDeserializer, FileSerializer
     * 
     */
    StringSerializer(const Dict& opts = Dict());
    ~StringSerializer();

    /** \brief Returns a string that holds the serialized objects
     * 
     * As a side effect, this method clears the internal buffer
    */
    std::string encode();
  };

  class CASADI_EXPORT FileSerializer : public SerializerBase {
  public:
    /** \brief Advanced serialization of CasADi objects
     * 
     * \seealso StringSerializer, FileDeserializer
     */
    FileSerializer(const std::string& fname, const Dict& opts = Dict());
    ~FileSerializer();
  };

  class CASADI_EXPORT StringDeserializer : public DeserializerBase {
  public:

    /** \brief Advanced deserialization of CasADi objects
     * 
     * \seealso StringDeserializer
     */
    StringDeserializer(const std::string& string);
    ~StringDeserializer();


    /** \brief Sets the string to deserialize objects from
    */
    void decode(const std::string& string);
  };

  class CASADI_EXPORT FileDeserializer : public DeserializerBase {
  public:
     /** \brief Advanced deserialization of CasADi objects
     * 
     * \seealso FileSerializer
     */
    FileDeserializer(const std::string& fname);
    ~FileDeserializer();
  };

} // namespace casadi

#endif // CASADI_SERIALIZER_HPP
