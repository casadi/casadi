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

#ifndef SX_FUNCTION_INTERNAL_HPP
#define SX_FUNCTION_INTERNAL_HPP

#include "sx_function.hpp"
#include "x_function_internal.hpp"

#ifdef WITH_LLVM
// Some forward declarations
namespace llvm{
  class Module;
  class Function;
} // namespace llvm
#endif // WITH_LLVM

namespace CasADi{

  struct AlgElData{
    // Partial derivatives
    double d[2];
};

/** \brief  Internal node class for SXFunction
  A regular user should never work with any Node class. Use SXFunction directly.
  \author Joel Andersson 
  \date 2010
*/
class SXFunctionInternal : public XFunctionInternalCommon<SXFunctionInternal,Matrix<SX>,SXNode>{
  friend class SXFunction;
  
  protected:
    /** \brief  Constructor (only to be called from SXFunction, therefore protected) */
    SXFunctionInternal(const std::vector<Matrix<SX> >& inputv, const std::vector<Matrix<SX> >& outputv);

  public:

  /** \brief  Make a deep copy */
  virtual SXFunctionInternal* clone() const;
    
  /** \brief  Destructor */
  virtual ~SXFunctionInternal();

  /** \brief  Evaluate the function numerically */
  virtual void evaluate(int nfdir, int nadir);

  /** \brief  evaluate symbolically while also propagating directional derivatives */
  virtual void evalSX(const std::vector<SXMatrix>& input, std::vector<SXMatrix>& output, 
                      const std::vector<std::vector<SXMatrix> >& fwdSeed, std::vector<std::vector<SXMatrix> >& fwdSens, 
                      const std::vector<std::vector<SXMatrix> >& adjSeed, std::vector<std::vector<SXMatrix> >& adjSens,
                      bool output_given, bool eliminate_constants);
                          
  /** \brief  Check if smooth */
  bool isSmooth() const;

  /** \brief  Print the algorithm */
  virtual void print(std::ostream &stream) const;

  /** \brief Calculate the jacobian of a number of function outputs with respect to a number of function inputs, optionally include the function outputs */
  virtual FX jacobian(const std::vector<std::pair<int,int> >& jblocks);

  /** \brief Hessian of output oind with respect to input iind */
  virtual FX hessian(int iind=0, int oind=0);

  /** \brief Calculate the expression for the jacobian of a number of function outputs with respect to a number of function inputs, optionally include the function outputs */
  std::vector<Matrix<SX> > jac(const std::vector<std::pair<int,int> >& jblocks, bool compact=false, const std::vector<bool>& symmetric_block=std::vector<bool>());
  
  /** \brief  DATA MEMBERS */
  
  /** \brief  Indices of the nodes corresponding to the inputs */
  std::vector<std::vector<int> > input_ind_;
  
  /** \brief  Indices of the nodes corresponding the non-zeros of the outputs */
  std::vector<std::vector<int> > output_ind_;

  /** \brief  An elemenent of the algorithm, namely a binary operation */
  typedef SXAlgEl AlgEl;
  
  /** \brief  all binary nodes of the tree in the order of execution */
  std::vector<AlgEl> algorithm_;
  std::vector<AlgElData> pder_;

  /** \brief  Working vector for numeric calculation */
  std::vector<double> work_;
  std::vector<double> dwork_;
  int worksize_;

  /// work vector for symbolic calculations (allocated first time)
  std::vector<SX> swork_;
  std::vector<SX> free_vars_;
  std::vector<int> refcount_;
  
  /// The expressions corresponding to each binary operation
  std::vector<SX> binops_;
  
  /** \brief  Initialize */
  virtual void init();

  /** \brief  Update the number of sensitivity directions during or after initialization */
  virtual void updateNumSens(bool recursive);

  /** \brief  Print to a c file */
  static void printVector(std::ostream &cfile, const std::string& name, const std::vector<int>& v);

  /** \brief  Print operation i to a stream */
  void printOperation(std::ostream &stream, int i) const;
  
  /** \brief  Print to a c file */
  void generateCode(const std::string& filename);
      
  // Evaluate with inplace operations (experimental)
  bool evaluate_inplace_;
  
  /** \brief Clear the function from its symbolic representation, to free up memory, no symbolic evaluations are possible after this */
  void clearSymbolic();

  /// Reset the virtual machine for sparsity calculations
  void spReset(int iind, int oind);

  /// Propagate the sparsity seeds
  void spProp(bool fwd);
  
  /// Get the forward/adjoint sparsity seed
  inline bvec_t& spGet(bool get_input, int ind, int sdir){
    if(get_input){
      return iwork_[input_ind_[ind][sdir]];
    } else {
      return iwork_[output_ind_[ind][sdir]];
    }
  }

  /// Work vector for sparsity detection
  bvec_t *iwork_;
  
  /// With just-in-time compilation
  bool just_in_time_;
  
  #ifdef WITH_LLVM
  llvm::Module *jit_module_;
  llvm::Function *jit_function_;
  #endif // WITH_LLVM
};


} // namespace CasADi

#endif // SX_FUNCTION_INTERNAL_HPP
