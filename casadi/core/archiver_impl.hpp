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


#ifndef CASADI_ARCHIVER_IMPL_HPP
#define CASADI_ARCHIVER_IMPL_HPP

#include "archiver.hpp"
#include "plugin_interface.hpp"
#include <iostream>

/// \cond INTERNAL
namespace casadi {

  /** \brief Archiver interface

    Zip and unzip archives without a core dependency

    @copydoc Archiver_doc
    \author Joris Gillis
    \date 2025

      \identifier{2d9} */
  class CASADI_EXPORT
  Archiver : public PluginInterface<Archiver> {
  public:
  typedef bool (* Unpack)(const std::string& src,
    const std::string& target_dir);
  typedef bool (* UnpackFromStringStream)(std::stringstream& src,
      const std::string& target_dir);
  typedef bool (* Pack)(const std::string& src_dir,
    const std::string& path);
  typedef bool (* PackToStream)(const std::string& src_dir,
    std::ostream& src);

    // Creator function for internal class
    typedef Archiver* (*Creator)();

    static const std::string meta_doc;

    // No static functions exposed
    struct Exposed{
      Unpack unpack;
      UnpackFromStringStream unpack_from_stringstream;
      Pack pack;
      PackToStream pack_to_stream;
    };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

    /// Infix
    static const std::string infix_;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_ARCHIVER_IMPL_HPP
