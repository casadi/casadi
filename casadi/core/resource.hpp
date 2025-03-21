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


#ifndef CASADI_RESOURCE_HPP
#define CASADI_RESOURCE_HPP

#include "shared_object.hpp"
#include "printable.hpp"

namespace casadi {
  // Forward declaration
  class ResourceInternal;
  class SerializingStream;
  class DeserializingStream;

  /** \brief RAII class for reading from a zip file */
  class CASADI_EXPORT Resource
    : public SharedObject,
      public SWIG_IF_ELSE(PrintableCommon, Printable<Resource>) {
  public:
    /** \brief Initialize with a path 
    *
    * If the path is a directory or empty, the path is passed through to the consumer.
    * Otherwise, the zip file is extracted to a temporary directory.
    *
    * Upon destruction, the temporary directory is removed.
    */
    Resource(const std::string& path);
    /// Default constructor
    explicit Resource();

    /// Readable name of the public class
    static std::string type_name() {return "Resource";}

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);

#ifndef SWIG
      /** \brief  Create from node */
      static Resource create(ResourceInternal *node);

      ResourceInternal* get() const;

      /// Access a member function or object
      const ResourceInternal* operator->() const;

      /// Reference to internal structure
      const ResourceInternal& operator*() const;
#endif

      /// Get path for a consumer
      const std::string& path() const;

    /** \brief Serialize an object */
    void serialize(SerializingStream &s) const;

    /** \brief Deserialize with type disambiguation */
    static Resource deserialize(DeserializingStream& s);
  };

} // namespace casadi

#endif // CASADI_RESOURCE_HPP
