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


#include "resource_internal.hpp"
#include "casadi_misc.hpp"
#include "archiver_impl.hpp"

#include "filesystem_impl.hpp"


namespace casadi {


ResourceInternal::ResourceInternal() {

}

DirResource::DirResource(const std::string& path)
    : ResourceInternal(),
      path_(path) {
  path_ = path;
}

void DirResource::disp(std::ostream& stream, bool more) const {
  stream << "DirResource(\"" << path_ << "\")";
}

DirResource::~DirResource() {
}

// serialize: link|EMBED

void ZipResource::unpack() {
  casadi_assert(Filesystem::is_enabled(),
  "Unzipping '" + path_ + "' requires advanced filesystem access. Compile CasADi with WITH_GC=ON.\n"
  "Alternatively, manually unzip it into a direcory, "
  "and pass this directory name instead of the zip file name.");

  // Extract filename part of path
  std::string zip_file = Filesystem::filename(path_);

  lock_file_ = temporary_file(zip_file + ".", ".lock");
  dir_ = lock_file_.substr(0, lock_file_.size()-5) + ".unzipped";

  casadi_assert(Archiver::has_plugin("libzip"),
  "Unzipping '" + path_ + "' requires libzip. Compile CasADi with WITH_LIBZIP=ON.\n"
  "Alternatively, manually unzip it into a direcory, "
  "and pass this directory name instead of the zip file name.");


  Archiver::getPlugin("libzip").exposed.unpack(path_, dir_);
}

/*void ZipMemResource::unpack() {
  casadi_assert(Filesystem::is_enabled(),
  "Unzipping '" + path_ + "' requires advanced filesystem access. Compile CasADi with WITH_GC=ON.\n"
  "Alternatively, manually unzip it into a direcory, and pass this directory name instead of the zip file name.");

  // Extract filename part of path
  std::string zip_file = Filesystem::filename(path_);

  lock_file_ = temporary_file(zip_file + ".", ".lock");
  dir_ = lock_file_.substr(0, lock_file_.size()-5) + ".unzipped";

  casadi_assert(Archiver::has_plugin("libzip"),
  "Unzipping '" + path_ + "' requires libzip. Compile CasADi with WITH_LIBZIP=ON.\n"
  "Alternatively, manually unzip it into a direcory, and pass this directory name instead of the zip file name.");

  Archiver::getPlugin("libzip").exposed.unpack_from_stream(blob_, dir_);
}*/

ZipResource::ZipResource(const std::string& path)
    : ResourceInternal() {
  path_ = path;
  unpack();
}

/*ZipMemResource::ZipMemResource(const std::istream& src)
    : ResourceInternal() {
  blob_ << src.rdbuf();
  unpack();
}*/

void ZipResource::disp(std::ostream& stream, bool more) const {
  stream << "ZipResource(\"" << path_ << "\") -> \"" << dir_ << "\"";
}

ZipResource::~ZipResource() {
  uout() << "~ZipResource" << std::endl;
  try {
      Filesystem::remove_all(dir_);
  } catch (...) {
      casadi_warning("Error: Cannot remove temporary directory: " + dir_);
  }
  try {
      Filesystem::remove(lock_file_);
  } catch (...) {
      casadi_warning("Error: Cannot remove lock file: " + lock_file_);
  }
}

void ResourceInternal::serialize(SerializingStream& s) const {
  s.version("ResourceInternal", 1);
  serialize_type(s);
  serialize_body(s);
}

void ResourceInternal::serialize_type(SerializingStream& s) const {
  s.pack("ResourceInternal::type", class_name());
}

void ResourceInternal::serialize_body(SerializingStream& s) const {
  s.version("ResourceInternal", 1);
}

ResourceInternal::ResourceInternal(DeserializingStream& s) {
  s.version("ResourceInternal", 1);
}

ResourceInternal* ResourceInternal::deserialize(DeserializingStream& s) {
  std::string class_name;
  s.unpack("ResourceInternal::type", class_name);
  if (class_name=="DirResource") {
    return DirResource::deserialize(s);
  } else if (class_name=="ZipResource") {
    return ZipResource::deserialize(s);
  } else  {
     casadi_error("Cannot deserialize type '" + class_name + "'");
  }
}

ResourceInternal* ZipResource::deserialize(DeserializingStream& s) {
  ZipResource* ret = new ZipResource(s);
  return ret;
}

ResourceInternal* DirResource::deserialize(DeserializingStream& s) {
  DirResource* ret = new DirResource(s);
  return ret;
}

ZipResource::ZipResource(DeserializingStream& s) : ResourceInternal(s) {
  s.unpack("ZipResource::path", path_);
  unpack();
}

DirResource::DirResource(DeserializingStream& s) : ResourceInternal(s) {

}

void ZipResource::serialize_body(SerializingStream& s) const {
  ResourceInternal::serialize_body(s);
  s.pack("ZipResource::path", path_);

  std::ifstream binary(path_, std::ios_base::binary);
  casadi_assert(binary.good(),
    "Could not open zip file '" + path_ + "'.");
  s.pack("ZipResource::archive", binary);
}

void DirResource::serialize_body(SerializingStream& s) const {
  ResourceInternal::serialize_body(s);
  s.pack("DirResource::path", path_);
}

} // namespace casadi
