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


#ifndef CASADI_SX_FUNCTION_HPP
#define CASADI_SX_FUNCTION_HPP

#include "x_function.hpp"

#ifdef WITH_OPENCL
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
#endif // WITH_OPENCL
/// \cond INTERNAL

namespace casadi {
  /** \brief  An atomic operation for the SXElem virtual machine */
  struct ScalarAtomic {
    int op;     /// Operator index
    int i0;
    union {
      double d;
      struct { int i1, i2; };
    };
  };

#ifdef WITH_OPENCL
  /** \brief Singleton for the sparsity propagation kernel
      TODO: Move to a separate file and make non sparsity pattern specific
      \author Joel Andersson
      \date 2013
  */
  class SparsityPropagationKernel {
  public:
    // Default constructor
    SparsityPropagationKernel();

    // Destructor
    ~SparsityPropagationKernel();

    // Copy constructor and equality operator
    // (not implemented, declared to prevent use of the default ones)
    SparsityPropagationKernel(const SparsityPropagationKernel& sparsityPropagationKernel);
    SparsityPropagationKernel&
      operator=(const SparsityPropagationKernel& sparsityPropagationKernel);

    // Data members (all public)
    cl_device_id device_id;
    cl_context context;
    cl_command_queue command_queue;
    cl_platform_id platform_id;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
  };
#endif // WITH_OPENCL

/** \brief  Internal node class for SXFunction
    Do not use any internal class directly - always use the public Function
    \author Joel Andersson
    \date 2010-2015
*/
class CASADI_EXPORT SXFunction :
        public XFunction<SXFunction, Matrix<SXElem>, SXNode>{
  public:
    /** \brief Constructor */
    SXFunction(const std::string& name,
                       const std::vector<Matrix<SXElem> >& inputv,
                       const std::vector<Matrix<SXElem> >& outputv);

  /** \brief  Destructor */
  virtual ~SXFunction();

  /** \brief  Evaluate numerically, work vectors given */
  virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

  /** \brief  evaluate symbolically while also propagating directional derivatives */
  virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

  /** \brief Calculate forward mode directional derivatives */
  virtual void evalFwd(const std::vector<std::vector<SX> >& fseed,
                       std::vector<std::vector<SX> >& fsens);

  /** \brief Calculate reverse mode directional derivatives */
  virtual void evalAdj(const std::vector<std::vector<SX> >& aseed,
                       std::vector<std::vector<SX> >& asens);

  /** \brief Create call to (cached) derivative function, forward mode  */
  virtual void forward_sx(const std::vector<SX>& arg, const std::vector<SX>& res,
                          const std::vector<std::vector<SX> >& fseed,
                          std::vector<std::vector<SX> >& fsens,
                          bool always_inline, bool never_inline);

  /** \brief Create call to (cached) derivative function, reverse mode  */
  virtual void reverse_sx(const std::vector<SX>& arg, const std::vector<SX>& res,
                          const std::vector<std::vector<SX> >& aseed,
                          std::vector<std::vector<SX> >& asens,
                          bool always_inline, bool never_inline);

  /** \brief  Check if smooth */
  bool is_smooth() const;

  /** \brief  Print the algorithm */
  virtual void print(std::ostream &stream) const;

  /** \brief Get type name */
  virtual std::string type_name() const;

  /** \brief Check if the function is of a particular type */
  virtual bool is_a(const std::string& type, bool recursive) const;

  /** \brief Gradient expression */
  virtual SX grad_sx(int iind=0, int oind=0);

  /** \brief Tangent expression */
  virtual SX tang_sx(int iind=0, int oind=0);

  /** \brief Jacobian expression */
  virtual SX jac_sx(int iind=0, int oind=0, bool compact=false, bool symmetric=false,
                    bool always_inline=true, bool never_inline=false);

  /** \brief Hessian expression */
  virtual SX hess_sx(int iind, int oind);

  ///@{
  /** \brief Get function input(s) and output(s)  */
  virtual const SX sx_in(int ind) const;
  virtual const std::vector<SX> sx_in() const;
  ///@}

  /// Get free variables (SX)
  virtual SX free_sx() const {return free_vars_;}

  /** \brief Does the function have free variables */
  virtual bool has_free() const { return !free_vars_.empty();}

  /** \brief Hessian (forward over adjoint) via source code transformation */
  SX hess(int iind=0, int oind=0);

  /** \brief Get the number of atomic operations */
  virtual int getAlgorithmSize() const { return algorithm_.size();}

  /** \brief Get the length of the work vector */
  virtual int getWorkSize() const { return sz_w();}

  /** \brief Get an atomic operation operator index */
  virtual int getAtomicOperation(int k) const { return algorithm_.at(k).op;}

  /** \brief Get the (integer) input arguments of an atomic operation */
  virtual std::pair<int, int> getAtomicInput(int k) const {
    const ScalarAtomic& atomic = algorithm_.at(k);
    return std::pair<int, int>(atomic.i1, atomic.i2);
  }

  /** \brief Get the floating point output argument of an atomic operation */
  virtual double getAtomicInputReal(int k) const {
    return algorithm_.at(k).d;
  }

  /** \brief Get the (integer) output argument of an atomic operation */
  virtual int getAtomicOutput(int k) const { return algorithm_.at(k).i0;}

  /** \brief Number of nodes in the algorithm */
  virtual int n_nodes() const { return algorithm_.size() - nnz_out();}

  /** \brief  DATA MEMBERS */

  /** \brief  An element of the algorithm, namely a binary operation */
  typedef ScalarAtomic AlgEl;

  /** \brief  An element of the tape */
  template<typename T>
  struct TapeEl {
    T d[2];
  };

  /** \brief  all binary nodes of the tree in the order of execution */
  std::vector<AlgEl> algorithm_;

  /// work vector for symbolic calculations (allocated first time)
  std::vector<SXElem> s_work_;
  std::vector<SXElem> free_vars_;

  /// The expressions corresponding to each binary operation
  std::vector<SXElem> operations_;

  /// The expressions corresponding to each constant
  std::vector<SXElem> constants_;

  /// Default input values
  std::vector<double> default_in_;

  ///@{
  /** \brief Options */
  static Options options_;
  virtual const Options& get_options() const { return options_;}
  ///@}

  /** \brief  Initialize */
  virtual void init(const Dict& opts);

  /** \brief Generate code for the declarations of the C function */
  virtual void generateDeclarations(CodeGenerator& g) const;

  /** \brief Generate code for the body of the C function */
  virtual void generateBody(CodeGenerator& g) const;

  /** \brief  Propagate sparsity forward */
  virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

  /** \brief  Propagate sparsity backwards */
  virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

  /// Is the class able to propagate seeds through the algorithm?
  virtual bool spCanEvaluate(bool fwd) { return true;}

  /** \brief Return Jacobian of all input elements with respect to all output elements */
  virtual Function getFullJacobian();

  /** \brief Get default input value */
  virtual double default_in(int ind) const { return default_in_.at(ind);}

  /// With just-in-time compilation using OpenCL
  bool just_in_time_opencl_;

  /// With just-in-time compilation for the sparsity propagation
  bool just_in_time_sparsity_;

#ifdef WITH_OPENCL
  // Initialize sparsity propagation using OpenCL
  void allocOpenCL();

  // Propagate sparsity using OpenCL
  void evaluateOpenCL();

  // Free memory for sparsity propagation using OpenCL
  void freeOpenCL();

  // Initialize sparsity propagation using OpenCL
  void spAllocOpenCL();

  // Propagate sparsity using OpenCL
  void spEvaluateOpenCL(bool fwd);

  // Free memory for sparsity propagation using OpenCL
  void spFreeOpenCL();

  // Compile OpenCL program
  static void compileProgram(cl_program program);

  // Execute OpenCL kernel
  static void executeKernel(cl_kernel kernel);

  // OpenCL memory object for the numerical evaluation
  cl_program program_;

  // OpenCL memory object for the sparsity propagation
  cl_program sp_program_;

  // Buffers and kernels for numerical evaluation
  std::vector<cl_mem> input_memobj_, output_memobj_;
  cl_kernel kernel_;

  // Buffers and kernels for sparsity propagation
  std::vector<cl_mem> sp_input_memobj_, sp_output_memobj_;
  cl_kernel sp_fwd_kernel_, sp_adj_kernel_;

  // OpenCL context. TODO: Nothing class specific in this class, move to a central location
  static SparsityPropagationKernel sparsity_propagation_kernel_;

#endif // WITH_OPENCL
};


} // namespace casadi

/// \endcond
#endif // CASADI_SX_FUNCTION_HPP
