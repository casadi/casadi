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


#ifndef CASADI_CASADI_META_HPP
#define CASADI_CASADI_META_HPP

#include <string>
#include "casadi_common.hpp"

namespace casadi {
  /**
  * \brief Collects global CasADi meta information
  *
  *  \author Joris Gillis
  *  \date 2012
  */
  class CASADI_EXPORT CasadiMeta {
    private:
      /// No instances are allowed
      CasadiMeta();
    public:
#ifndef SWIG
      static const std::string version;
      static const std::string git_revision;
      static const std::string git_describe;
      static const std::string feature_list;
      static const std::string build_type;
      static const std::string compiler_id;
      static const std::string compiler;
      static const std::string compiler_flags;
      static const std::string modules;
      static const std::string plugins;
      static const std::string install_prefix;
#endif //SWIG
    /** \brief Obtain the version number of CasADi
    *  The format is 'x.y.z' or 'x.y.z+'
    *
    *  The variant without + indicates that the version is an official release
    *
    *  The variant with + indicates that the version is more recent than 'x.y.z',
    *     and might be more recent than 'x.y.w'  with w>z.
    *
    *  \see getGitRevision getGitDescribe
    */
    static std::string getVersion() { return version; }
    /** \brief Obtain the git hash of this build
    *      (only available if built from a git repo)
    */
    static std::string getGitRevision() { return git_revision; }
    /** \brief Obtain the git description of this build
    *      (only available if built from a git repo)
    */
    static std::string getGitDescribe() { return git_describe; }
    /** \brief Obtain list of features that were compiled into this build
    */
    static std::string getFeatureList() { return feature_list; }
    /** \brief Obtain build type: RELEASE/Debug
    */
    static std::string getBuildType() { return build_type; }
    /** \brief Obtain compiler identification
    * Provided by http://www.cmake.org/cmake/help/v2.8.10/cmake.html#variable:CMAKE_LANG_COMPILER_ID
    */
    static std::string getCompilerId() { return compiler_id; }
    /** \brief Obtain compiler
    */
    static std::string getCompiler() { return compiler; }
    /** \brief Obtain compiler flags
    */
    static std::string getCompilerFlags() { return compiler_flags; }
    /** \brief Obtain modules list
    */
    static std::string getModules() { return modules; }
    /** \brief Obtain plugins list
    */
    static std::string getPlugins() { return plugins; }
    /** \brief Obtain install prefix
    */
    static std::string getInstallPrefix() { return install_prefix; }
  };

}  // namespace casadi

#endif // CASADI_CASADI_META_HPP
