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


#ifndef CASADI_CASADI_TYPES_HPP
#define CASADI_CASADI_TYPES_HPP

#include <climits>
#include <vector>
#include <utility>

#include "casadi_common.hpp"
#include "casadi_logger.hpp"

#ifdef SWIG
#define SWIG_OUTPUT(arg) OUTPUT
#define SWIG_INOUT(arg) INOUT
#define SWIG_CONSTREF(arg) const arg
#ifdef SWIGMATLAB
#define SWIG_IND1 true
#else // SWIGMATLAB
#define SWIG_IND1 false
#endif // SWIGMATLAB
#else // SWIG
#define SWIG_OUTPUT(arg) arg
#define SWIG_INOUT(arg) arg
#define SWIG_CONSTREF(arg) const arg &
#define SWIG_IND1 false
#endif // SWIG

namespace casadi {

  /// Forward declarations
  class SXElement;
  class MX;
  template<class T> class Matrix;
  class Function;
  class Sparsity;
  class CodeGenerator;
  class NlpSolver;
  class LinearSolver;
  class Integrator;
  class QpSolver;
  class StabilizedQpSolver;
  class QcqpSolver;
  class LpSolver;
  class SdpSolver;
  class SocpSolver;
  class SdqpSolver;
  class DleSolver;
  class CleSolver;
  class ImplicitFunction;

  class DerivativeGenerator;
  class Callback;
  class CustomEvaluate;
  class CustomFunction;
  class NlpBuilder;
  class Variable;
  class DaeBuilder;
  class XmlFile;

#ifndef SWIG

  // Type with a size corresponding to that of double (or smaller) that can be used to hold a set
  // of booleans. If the compiler supports C99 or has defined __SIZEOF_LONG_LONG__,
  // we shall use the long long datatype, which is 64 bits, otherwise long
  #if (defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L || defined(__SIZEOF_LONG_LONG__))
  typedef unsigned long long bvec_t;
  #else
  typedef unsigned long bvec_t;
  #endif

  // Number of directions we can deal with at a time
  // the size of bvec_t in bits (CHAR_BIT is the number of bits per byte, usually 8)
  const int bvec_size = CHAR_BIT*sizeof(bvec_t);

  // Make sure that the integer datatype is indeed smaller or equal to the double
  //assert(sizeof(bvec_t) <= sizeof(double)); // doesn't work - very strange
#endif // SWIG

} // namespace casadi

#endif // CASADI_CASADI_TYPES_HPP

