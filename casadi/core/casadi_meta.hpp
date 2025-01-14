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


#ifndef CASADI_CASADI_META_HPP
#define CASADI_CASADI_META_HPP

#include <casadi/core/casadi_export.h>
#include <string>

namespace casadi {
  /** \brief Collects global CasADi meta information
  *
  *  \author Joris Gillis
  *  \date 2012

      \identifier{23k} */
  class CASADI_EXPORT CasadiMeta {
    private:
      /// No instances are allowed
      CasadiMeta();
    public:
    /** \brief Obtain the version number of CasADi

    *  The format is 'x.y.z' or 'x.y.z+'
    *
    *  The variant without + indicates that the version is an official release
    *
    *  The variant with + indicates that the version is more recent than 'x.y.z',
    *     and might be more recent than 'x.y.w'  with w>z.
    *
    *  \see getGitRevision getGitDescribe

        \identifier{iz} */
    static const char* version();
    /** \brief Obtain the git hash of this build

    *      (only available if built from a git repo)

        \identifier{j0} */
    static const char* git_revision();
    /** \brief Obtain the git description of this build

    *      (only available if built from a git repo)

        \identifier{j1} */
    static const char* git_describe();
    /** \brief Obtain list of features that were compiled into this build

        \identifier{j2} */
    static const char* feature_list();
    /** \brief Obtain build type: RELEASE/Debug

        \identifier{j3} */
    static const char* build_type();
    /** \brief Obtain compiler identification

    * Provided by http://www.cmake.org/cmake/help/v2.8.10/cmake.html#variable:CMAKE_LANG_COMPILER_ID

        \identifier{j4} */
    static const char* compiler_id();
    /** \brief Obtain compiler

        \identifier{j5} */
    static const char* compiler();
    /** \brief Obtain compiler flags

        \identifier{j6} */
    static const char* compiler_flags();
    /** \brief Obtain modules list

        \identifier{j7} */
    static const char* modules();
    /** \brief Obtain plugins list

        \identifier{j8} */
    static const char* plugins();
    /** \brief Obtain install prefix

        \identifier{j9} */
    static const char* install_prefix();

    /** \brief Obtain shared library prefix

        \identifier{2be} */
    static const char* shared_library_prefix();

    /** \brief Obtain shared library suffix

        \identifier{2bf} */
    static const char* shared_library_suffix();

    /** \brief Obtain object file suffix

        \identifier{2bg} */
    static const char* object_file_suffix();
  };

}  // namespace casadi

#endif // CASADI_CASADI_META_HPP
