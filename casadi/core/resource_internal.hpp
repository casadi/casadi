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


#ifndef CASADI_RESOURCE_INTERNAL_HPP
#define CASADI_RESOURCE_INTERNAL_HPP

#include "resource.hpp"
#include "shared_object.hpp"
#include "serializing_stream.hpp"

/// \cond INTERNAL
namespace casadi {


  /** \brief RAII class base for reading from resources */
  class CASADI_EXPORT ResourceInternal : public SharedObjectInternal {
    public:
      /** \brief Initialize with a path 
      *
      * If the path is a directory or empty, the path is passed through to the consumer.
      * Otherwise, the zip file is extracted to a temporary directory.
      *
      * Upon destruction, the temporary directory is removed.
      */
      explicit ResourceInternal();
      /// Get path for a consumer
      virtual const std::string& path() const = 0;
      void serialize(SerializingStream& s) const;

      virtual void serialize_type(SerializingStream& s) const;
      virtual void serialize_body(SerializingStream& s) const;

      static ResourceInternal* deserialize(DeserializingStream& s);

    protected:
      explicit ResourceInternal(DeserializingStream& s);
  };

  class CASADI_EXPORT DirResource : public ResourceInternal {
    public:
      /** \brief Initialize with a path */
      DirResource(const std::string& path);
      ~DirResource() override;
      /// Get path for a consumer
      const std::string& path() const override {return path_;}

      /** \brief Get type name */
      std::string class_name() const override {return "DirResource";}

      /// Print description
      void disp(std::ostream& stream, bool more) const override;

      void serialize_body(SerializingStream& s) const override;

      static ResourceInternal* deserialize(DeserializingStream& s);
    private:
      std::string path_;
    protected:
      explicit DirResource(DeserializingStream& s);
  };

  /** \brief RAII class for reading from a zip file

  \identifier{2c5} */
  class CASADI_EXPORT ZipResource : public ResourceInternal {
    public:
      /** \brief Initialize with a path 
      *
      * If the path is a directory or empty, the path is passed through to the consumer.
      * Otherwise, the zip file is extracted to a temporary directory.
      *
      * Upon destruction, the temporary directory is removed.

          \identifier{2c6} */
      ZipResource(const std::string& path);
      ~ZipResource() override;
      /// Get path for a consumer
      const std::string& path() const override {return dir_;}

      void unpack();

      /** \brief Get type name */
      std::string class_name() const override {return "ZipResource";}

      /// Print description
      void disp(std::ostream& stream, bool more) const override;

      void serialize_body(SerializingStream& s) const override;

      static ResourceInternal* deserialize(DeserializingStream& s);
    private:
      std::string lock_file_;
      std::string dir_;
      std::string path_;
    protected:
      explicit ZipResource(DeserializingStream& s);
  };

  /*
  ** \brief RAII class for reading from a zip held in memory *
  class CASADI_EXPORT ZipMemResource : public ResourceInternal {
    public:
    ZipMemResource(const std::stream& src);
      ~ZipMemResource() override;
      /// Get path for a consumer
      const std::string& path() const override {return dir_;}

      void unpack();
      ** \brief Get type name *
      std::string class_name() const override {return "ZipStreamResource";}

      /// Print description
      void disp(std::ostream& stream, bool more) const override;

      void serialize_body(SerializingStream& s) const override;

      static ResourceInternal* deserialize(DeserializingStream& s);
    private:
      std::string lock_file_;
      std::string dir_;
      std::stringstrean blob_;
    protected:
      explicit ZipMemResource(DeserializingStream& s);
  };*/

} // namespace casadi
/// \endcond
#endif // CASADI_RESOURCE_INTERNAL_HPP
