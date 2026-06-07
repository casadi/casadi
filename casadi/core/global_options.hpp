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


#ifndef CASADI_GLOBAL_OPTIONS_HPP
#define CASADI_GLOBAL_OPTIONS_HPP

#include <fstream>
#include <iostream>

#include "casadi/core/casadi_common.hpp"
#include <casadi/core/casadi_export.h>

namespace casadi {

  /** \brief Collects global CasADi options
  *
  *
  * Note to developers:  \n
  *  - use sparingly. Global options are - in general - a rather bad idea \n
  *  - this class must never be instantiated. Access its static members directly \n
  *
  *  \author Joris Gillis
  *  \date 2012

      \identifier{23m} */
  class CASADI_EXPORT GlobalOptions {
    private:
      /// No instances are allowed
      GlobalOptions();
    public:

#ifndef SWIG
      /** \brief Indicates whether simplifications should be made on the fly.

      * e.g.   cos(-x) -> cos(x)
      * Default: true

          \identifier{17v} */
      static bool simplification_on_the_fly;

      static std::string casadipath;

      static std::string casadi_include_path;

      static bool hierarchical_sparsity;

      static casadi_int max_num_dir;

      static casadi_int start_index;

      static bool julia_initialized;

      static casadi_int copy_elision_min_size;

      static std::string temp_work_dir; // Temporary work directory

      /** \brief numpy interop mode (issue #2959).  Controls how an explicit
       *  `numpy.foo(M)` on a casadi value behaves in the Python bindings:
       *    0 (default): legacy casadi 3.7.2 behaviour + a Python
       *       FutureWarning -- a numeric value densifies to a numpy result,
       *       a symbolic value returns a casadi value.
       *    1: the casadi-aware numpy support -- `numpy.foo(M)` returns a
       *       python-only numpy-semantics array wrapper
       *       (casadi.ArrayInterface) following numpy's shape/axis contract
       *       for SX/MX/DM alike.
       *   -1: the same legacy behaviour as 0, but silent (no warning) --
       *       exact casadi 3.7.2 semantics.
       *  Operator arithmetic (`M + x`, `x * M`, ...) always stays casadi-
       *  typed and is unaffected.  Affects only Python bindings.
       *
       *  Temporary opt-in mechanism: exists to ease the transition and
       *  will be removed once numpy becomes the unconditional default.
       *  Process-global and not thread-safe -- it changes the *return type*
       *  of numpy dispatch, so set it once at startup rather than toggling
       *  it concurrently with running code.
       *  Probe with hasattr(casadi.GlobalOptions, "setNumpyMode").

          \identifier{2i0} */
      static int numpy_mode;

#endif //SWIG
      // Setter and getter for simplification_on_the_fly
      static void setSimplificationOnTheFly(bool flag) { simplification_on_the_fly = flag; }
      static bool getSimplificationOnTheFly() { return simplification_on_the_fly; }

      // Setter and getter for hierarchical_sparsity
      static void setHierarchicalSparsity(bool flag) { hierarchical_sparsity = flag; }
      static bool getHierarchicalSparsity() { return hierarchical_sparsity; }

      static void setCasadiPath(const std::string & path) { casadipath = path; }
      static std::string getCasadiPath() { return casadipath; }

      static void setCasadiIncludePath(const std::string & path) { casadi_include_path = path; }
      static std::string getCasadiIncludePath() { return casadi_include_path; }

      static void setMaxNumDir(casadi_int ndir) { max_num_dir=ndir; }
      static casadi_int getMaxNumDir() { return max_num_dir; }

      static void setCopyElisionMinSize(casadi_int sz) {
        copy_elision_min_size=sz;
        if (copy_elision_min_size!=-1 && copy_elision_min_size<2) {
          copy_elision_min_size = 2;
        }
      }
      static casadi_int getCopyElisionMinSize() { return copy_elision_min_size; }

      static void setTempWorkDir(const std::string& dir);
      static std::string getTempWorkDir() { return temp_work_dir; }

      /** \brief Set the numpy interop mode (issue #2959): 1 = casadi-aware
       *  numpy support, 0 (default) = legacy + FutureWarning, -1 = legacy
       *  but silent.  See the numpy_mode field.  Affects only Python
       *  bindings; ignored elsewhere.  Probe with
       *  hasattr(casadi.GlobalOptions, "setNumpyMode").

          \identifier{2i1} */
      static void setNumpyMode(int mode) { numpy_mode = mode; }

      /** \brief Get the current numpy interop mode.  See setNumpyMode.

          \identifier{2i2} */
      static int getNumpyMode() { return numpy_mode; }

  };

} // namespace casadi

#endif // CASADI_GLOBAL_OPTIONS_HPP
