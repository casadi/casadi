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

#ifndef CASADI_TYPES_HPP
#define CASADI_TYPES_HPP

#include <climits>
#include <cassert>
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
  /// \endcond
  typedef LPStructIOSchemeVector<Sparsity> LPStructure;
  typedef QPStructIOSchemeVector<Sparsity> QPStructure;
  typedef QCQPStructIOSchemeVector<Sparsity> QCQPStructure;
  typedef SDPStructIOSchemeVector<Sparsity> SDPStructure;
  typedef SOCPStructIOSchemeVector<Sparsity> SOCPStructure;
  typedef SDQPStructIOSchemeVector<Sparsity> SDQPStructure;
  class NLPSolver;
  class LinearSolver;
  class Integrator;
  class QPSolver;
  class StabilizedQPSolver;
  class QCQPSolver;
  class LPSolver;
  class SDPSolver;
  class SOCPSolver;
  class SDQPSolver;
  class ImplicitFunction;

  class DerivativeGenerator;
  class Callback;
  class CustomEvaluate;
  class CustomFunction;

  /// Function pointer to a nonlinear solver creator function
  typedef NLPSolver (*NLPSolverCreator)(const Function& nlp);

  /// Function pointer to a linear solver creator function
  typedef LinearSolver (*linearSolverCreator)(const Sparsity& sparsity, int nrhs);

  /// Function pointer to a LP solver creator function
  typedef LPSolver (*LPSolverCreator)(const LPStructure& st);

  /// Function pointer to an integrator creator function
  typedef Integrator (*integratorCreator)(const Function& f, const Function& g);

  /// Function pointer to a QP solver creator function
  typedef QPSolver (*QPSolverCreator)(const QPStructure& st);

  /// Function pointer to a Stabilized QP solver creator function
  typedef StabilizedQPSolver (*StabilizedQPSolverCreator)(const QPStructure& st);

  /// Function pointer to a QCQP solver creator function
  typedef QCQPSolver (*QCQPSolverCreator)(const QCQPStructure& st);

  /// Function pointer to an SDP solver creator function
  typedef SDPSolver (*SDPSolverCreator)(const SDPStructure& st);

  /// Function pointer to an SDQP solver creator function
  typedef SDQPSolver (*SDQPSolverCreator)(const SDQPStructure& st);

  /// Function pointer to an SOCP solver creator function
  typedef SOCPSolver (*SOCPSolverCreator)(const SOCPStructure& st);

  /// Function pointer to an implicit function creator
  typedef ImplicitFunction (*implicitFunctionCreator)(const Function& f, const Function& jac,
                                                      const LinearSolver& linsol);


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

#endif // CASADI_TYPES_HPP

