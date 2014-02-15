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

#ifndef CASADI_OPTIONS_HPP
#define CASADI_OPTIONS_HPP

#include <iostream>
#include <fstream>

namespace CasADi {
  /**
  * \brief Collects global CasADi options
  *
  *
  * Note to developers:  \n
  *  - use sparingly. Global options are - in general - a rather bad idea \n
  *  - this class must never be instantiated. Access its static members directly \n 
  *
  *  \author Joris Gillis 
  *  \date 2012
  */
  class CasadiOptions {
    private:
      /// No instances are allowed
      CasadiOptions();
    public:

#ifndef SWIG
      /** \brief Catch CasADi errors when they reach the python layer
      *  If set to true (which will always be the default), any CasADi errors are caught in the SWIG interface
      *  and converted to python exceptions. This allows the user to obtain python stacktraces, use try/except, etc...
      *
      *  If set to false, CasADi errors crash the python interface. This allows a CasADi developer
      *  to obtain a full C++ stacktrace with a debugger such as 'gdb'.
      *
      *  Default: true
      */
      static bool catch_errors_python;
      /** \brief Indicates wether simplifications should be made on the fly.
      * e.g   cos(-x) -> cos(x)
      * Default: true
      */
      static bool simplification_on_the_fly;
      
      /** \brief Stream on which profiling log should be written */
      static std::ofstream profilingLog;
      
      /** \brief flag to indicate if profiling is active */
      static bool profiling;
#endif //SWIG
      // Setter and getter for catch_errors_python
      static void setCatchErrorsPython(bool flag) { catch_errors_python = flag; }
      static bool getCatchErrorsPython() { return catch_errors_python; }

      // Setter and getter for simplification_on_the_fly
      static void setSimplificationOnTheFly(bool flag) { simplification_on_the_fly = flag; }
      static bool getSimplificationOnTheFly() { return simplification_on_the_fly; }
      
      /** \brief Start virtual machine profiling
      *
      *  When profiling is active, each primitive of an MX algorithm is profiling and dumped into the supplied file _filename_
      *  After the profiling is done, convert the supplied file to a viewable webpage with:
      * `python -mcasadi.tools.profilereport -o profiling.html _filename_`
      */
      static void startProfiling(const std::string &filename);
      static void stopProfiling();
  };

}

#endif //CASADI_OPTIONS_HPP
