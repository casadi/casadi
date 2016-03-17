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
  class SXElem;
  class MX;
  template<class T> class Matrix;
  class Function;
  class Sparsity;
  class CodeGenerator;
  class NlpBuilder;
  class Variable;
  class DaeBuilder;
  class XmlFile;
  class Compiler;

#ifndef SWIG
  // Workarond for MinGW bug
#if defined(__MINGW32__) || defined(__MINGW64__)
  template<typename T>
  std::string to_string(const T& n) {
    std::stringstream s;
    s << n;
    return s.str();
  }
#else
  using std::to_string;
#endif

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

  ///@{
  /** \brief  Function pointer types for the C API */
  typedef void (*signal_t)(void);
  typedef int (*getint_t)(void);
  typedef const char* (*name_t)(int i);
  typedef const int* (*sparsity_t)(int i);
  typedef int (*work_t)(int* sz_arg, int* sz_res, int* sz_iw, int* sz_w);
  typedef int (*eval_t)(const double** arg, double** res, int* iw, double* w, int mem);
  typedef void (*simple_t)(const double* arg, double* res);
  ///@}

  /// Inputs of the symbolic representation of the DAE
  enum DeIn {
    DE_T,
    DE_X,
    DE_Z,
    DE_P,
    DE_RX,
    DE_RZ,
    DE_RP,
    DE_NUM_IN};

  /// Shortnames for DAE symbolic representation inputs
  const std::vector<std::string> DE_INPUTS = {"t", "x", "z", "p", "rx", "rz", "rp"};

  /// Inputs of the symbolic representation of the DAE
  enum DeOut {
    DE_ODE,
    DE_ALG,
    DE_QUAD,
    DE_RODE,
    DE_RALG,
    DE_RQUAD,
    DE_NUM_OUT};

  /// Shortnames for DAE symbolic representation outputs
  const std::vector<std::string> DE_OUTPUTS = {"ode", "alg", "quad", "rode", "ralg", "rquad"};

  /// Input arguments of an ODE/DAE function
  enum DAEInput {
    /// Differential state
    DAE_X,
    /// Algebraic state
    DAE_Z,
    /// Parameter
    DAE_P,
    /// Explicit time dependence
    DAE_T,
    /// Number of arguments
    DAE_NUM_IN
  };

  /// Output arguments of an DAE function
  enum DAEOutput {
    /// Right hand side of the implicit ODE
    DAE_ODE,
    /// Right hand side of algebraic equations
    DAE_ALG,
    /// Right hand side of quadratures equations
    DAE_QUAD,
    /// Number of arguments
    DAE_NUM_OUT
  };

  /// Input arguments of an ODE/DAE backward integration function
  enum RDAEInput {
    /// Backward differential state
    RDAE_RX,
    /// Backward algebraic state
    RDAE_RZ,
    /// Backward  parameter vector
    RDAE_RP,
    /// Forward differential state
    RDAE_X,
    /// Forward algebraic state
    RDAE_Z,
    /// Parameter vector
    RDAE_P,
    /// Explicit time dependence
    RDAE_T,
    /// Number of arguments
    RDAE_NUM_IN
  };

  /// Output arguments of an ODE/DAE backward integration function
  enum RDAEOutput {
    /// Right hand side of ODE
    RDAE_ODE,
    /// Right hand side of algebraic equations
    RDAE_ALG,
    /// Right hand side of quadratures
    RDAE_QUAD,
    /// Number of arguments
    RDAE_NUM_OUT
  };

  /// Input arguments of an integrator
  enum IntegratorInput {
    /// Differential state at the initial time
    INTEGRATOR_X0,
    /// Parameters
    INTEGRATOR_P,
    /// Initial guess for the algebraic variable
    INTEGRATOR_Z0,
    /// Backward differential state at the final time
    INTEGRATOR_RX0,
    /// Backward parameter vector
    INTEGRATOR_RP,
    /// Initial guess for the backwards algebraic variable
    INTEGRATOR_RZ0,
    /// Number of input arguments of an integrator
    INTEGRATOR_NUM_IN
  };

  /// Output arguments of an integrator
  enum IntegratorOutput {
    /// Differential state at the final time
    INTEGRATOR_XF,
    /// Quadrature state at the final time
    INTEGRATOR_QF,
    /// Algebraic variable at the final time
    INTEGRATOR_ZF,
    /// Backward differential state at the initial time
    INTEGRATOR_RXF,
    /// Backward quadrature state at the initial time
    INTEGRATOR_RQF,
    /// Backward algebraic variable at the initial time
    INTEGRATOR_RZF,
    /// Number of output arguments of an integrator
    INTEGRATOR_NUM_OUT
  };

  /// Input arguments of an NLP function
  enum NLPInput {
    /// Decision variable
    NL_X,
    /// Fixed parameter
    NL_P,
    /// Number of NLP inputs
    NL_NUM_IN
  };

  /// Shortname for onput arguments of an NLP function
  const std::vector<std::string> NL_INPUTS = {"x", "p"};

  /// Output arguments of an NLP function
  enum NLPOutput {
    /// Objective function
    NL_F,
    /// Constraint function
    NL_G,
    /// Number of NLP outputs
    NL_NUM_OUT
  };

  /// Shortname for output arguments of an NLP function
  const std::vector<std::string> NL_OUTPUTS = {"f", "g"};

  /// Input arguments of an NLP objective gradient function]
  enum GradFInput {
    /// Decision variable
    GRADF_X,
    /// Fixed parameter
    GRADF_P,
    /// Number of inputs
    GRADF_NUM_IN
  };

  /// Output arguments of an NLP objective gradient function
  enum GradFOutput {
    /// Jacobian of the constraints
    GRADF_GRAD,
    /// Objective function
    GRADF_F,
    /// Constraint function
    GRADF_G,
    /// Number of outputs
    GRADF_NUM_OUT
  };

  /// Input arguments of an NLP Jacobian function
  enum JacGInput {
    /// Decision variable
    JACG_X,
    /// Fixed parameter
    JACG_P,
    /// Number of inputs
    JACG_NUM_IN
  };

  /// Output arguments of an NLP Jacobian function
  enum JacGOutput {
    /// Jacobian of the constraints
    JACG_JAC,
    /// Objective function
    JACG_F,
    /// Constraint function
    JACG_G,
    /// Number of outputs
    JACG_NUM_OUT
  };

  /// Input arguments of an NLP Hessian function
  enum HessLagInput {
    /// Decision variable
    HESSLAG_X,
    /// Fixed parameter
    HESSLAG_P,
    /// Multiplier for f. Just a scalar factor for the objective that the
    /// NLP solver might use to scale the objective
    HESSLAG_LAM_F,
    /// Multiplier for g
    HESSLAG_LAM_G,
    /// Number of inputs
    HESSLAG_NUM_IN
  };

  /// Output arguments of an NLP Hessian function
  enum HessLagOutput {
    /// Hessian of the Lagrangian
    HESSLAG_HESS,
    /// Objective function
    HESSLAG_F,
    /// Constraint function
    HESSLAG_G,
    /// Gradient of the Lagrangian with respect to x
    HESSLAG_GRAD_X,
    /// Gradient of the Lagrangian with respect to p
    HESSLAG_GRAD_P,
    /// Number of outputs
    HESSLAG_NUM_OUT
  };

  /// Input arguments of an NLP Solver
  enum NlpsolInput {
    /// Decision variables, initial guess (nx x 1)
    NLPSOL_X0,
    /// Value of fixed parameters (np x 1)
    NLPSOL_P,
    /// Decision variables lower bound (nx x 1), default -inf
    NLPSOL_LBX,
    /// Decision variables upper bound (nx x 1), default +inf
    NLPSOL_UBX,
    /// Constraints lower bound (ng x 1), default -inf
    NLPSOL_LBG,
    /// Constraints upper bound (ng x 1), default +inf
    NLPSOL_UBG,
    /// Lagrange multipliers for bounds on X, initial guess (nx x 1)
    NLPSOL_LAM_X0,
    /// Lagrange multipliers for bounds on G, initial guess (ng x 1)
    NLPSOL_LAM_G0,
    NLPSOL_NUM_IN
  };

  /// Output arguments of an NLP Solver
  enum NlpsolOutput {
    /// Decision variables at the optimal solution (nx x 1)
    NLPSOL_X,
    /// Cost function value at the optimal solution (1 x 1)
    NLPSOL_F,
    /// Constraints function at the optimal solution (ng x 1)
    NLPSOL_G,
    /// Lagrange multipliers for bounds on X at the solution (nx x 1)
    NLPSOL_LAM_X,
    /// Lagrange multipliers for bounds on G at the solution (ng x 1)
    NLPSOL_LAM_G,
    /// Lagrange multipliers for bounds on P at the solution (np x 1)
    NLPSOL_LAM_P,
    NLPSOL_NUM_OUT
  };

  /// Input arguments of a QP problem
  enum QpsolInput {
    /// The square matrix H: sparse, (n x n). Only the lower triangular part is actually used.
    /// The matrix is assumed to be symmetrical.
    QPSOL_H,
    /// The vector g: dense,  (n x 1)
    QPSOL_G,
    /// The matrix A: sparse, (nc x n) - product with x must be dense.
    QPSOL_A,
    /// dense, (nc x 1)
    QPSOL_LBA,
    /// dense, (nc x 1)
    QPSOL_UBA,
    /// dense, (n x 1)
    QPSOL_LBX,
    /// dense, (n x 1)
    QPSOL_UBX,
    /// dense, (n x 1)
    QPSOL_X0,
    /// dense
    QPSOL_LAM_X0,
    QPSOL_NUM_IN};

  /// Output arguments of an QP Solver
  enum QpsolOutput {
    /// The primal solution
    QPSOL_X,
    /// The optimal cost
    QPSOL_COST,
    /// The dual solution corresponding to linear bounds
    QPSOL_LAM_A,
    /// The dual solution corresponding to simple bounds
    QPSOL_LAM_X,
    QPSOL_NUM_OUT};

  /// Input arguments of a linear solver
  enum LinsolInput {
    /// The square matrix A: sparse, (n x n)
    LINSOL_A,
    /// The right-hand-side matrix b: dense,  (n x m)
    LINSOL_B,
    LINSOL_NUM_IN};

  /// Output arguments of a linear solver
  enum LinsolOutput {
    /// Solution to the linear system of equations
    LINSOL_X,
    LINSOL_NUM_OUT};

#endif // SWIG

} // namespace casadi

#endif // CASADI_CASADI_TYPES_HPP

