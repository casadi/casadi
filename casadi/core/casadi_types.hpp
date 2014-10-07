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

namespace casadi {

  /// Forward declarations
  class SXElement;
  class MX;
  template<class T> class Matrix;
  class Function;
  class Sparsity;
  /// \cond INTERNAL
  template<class T> class LPStructIOSchemeVector;
  template<class T> class QPStructIOSchemeVector;
  template<class T> class QCQPStructIOSchemeVector;
  template<class T> class SDPStructIOSchemeVector;
  template<class T> class SOCPStructIOSchemeVector;
  template<class T> class SDQPStructIOSchemeVector;
  template<class T> class DleStructIOSchemeVector;
  template<class T> class DpleVecStructIOSchemeVector;
  template<class T> class LrDleStructIOSchemeVector;
  template<class T> class LrDpleVecStructIOSchemeVector;
  template<class T> class CleStructIOSchemeVector;
  /// \endcond
  typedef LPStructIOSchemeVector<Sparsity> LPStructure;
  typedef QPStructIOSchemeVector<Sparsity> QPStructure;
  typedef QCQPStructIOSchemeVector<Sparsity> QCQPStructure;
  typedef SDPStructIOSchemeVector<Sparsity> SDPStructure;
  typedef SOCPStructIOSchemeVector<Sparsity> SOCPStructure;
  typedef SDQPStructIOSchemeVector<Sparsity> SDQPStructure;
  typedef DleStructIOSchemeVector<Sparsity> DleStructure;
  typedef LrDleStructIOSchemeVector<Sparsity> LrDleStructure;
  typedef CleStructIOSchemeVector<Sparsity> CleStructure;
  typedef DpleVecStructIOSchemeVector< std::vector<Sparsity> > DpleStructure;
  typedef LrDpleVecStructIOSchemeVector< std::vector<Sparsity> > LrDpleStructure;
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
  class SymbolicNLP;
  class Variable;
  class SymbolicOCP;
  class XmlFile;

#ifndef SWIG
  // The number of derivative directions for which the tool has been optimized
  const int optimized_num_dir = 64;

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

