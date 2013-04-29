/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef CASADI_META_HPP
#define CASADI_META_HPP

#include <string>

namespace CasADi {
  /**
  * \brief Collects global CasADi meta information
  *
  *  \author Joris Gillis 
  *  \date 2012
  */
  class CasadiMeta {
    private:
      /// No instances are allowed
      CasadiMeta();
    public:
#ifndef SWIG
      static const std::string version;
      static const std::string git_revision;
      static const std::string git_describe;
#endif //SWIG
    /** \brief Obtain the version number of CasADi
    *  The format is x.y.z or x.y.z +
    *
    *  The variant without + indicates that the verion is an official release  
    *
    *  The variant with + indicates that the version is more recent than x.y.z,
    *     and might be more recent than x.y.w  with w>z.
    *     
    *  \see getGitRevision getGitDescribe
    */
    static std::string getVersion() { return version; }
    /** \brief Obtain the git hash of this build
    *      (only available if built from a git repo )
    */
    static std::string getGitRevision() { return git_revision; }
    /** \brief Obtain the git description of this build
    *      (only available if built from a git repo )
    */
    static std::string getGitDescribe() { return git_describe; }
  };

}

#endif //CASADI_META_HPP
