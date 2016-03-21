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


#ifndef CASADI_KERNEL_SUM_HPP
#define CASADI_KERNEL_SUM_HPP

#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** KernelSum2D statement
      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT KernelSum : public FunctionInternal {
    friend class KernelSum2D;
  public:

    /** \brief Constructor (generic kernel_sum_2d) */
    KernelSum(const std::string& name, const Function& f,
                        const std::pair<int, int> & size,
                        double r,
                        int n);

    /** \brief  Destructor */
    virtual ~KernelSum();

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() { return f_.n_in()-1;}
    virtual size_t get_n_out() { return f_.n_out();}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int i) {
      if (i==0) {
        return Sparsity::dense(size_);
      } else if (i==1) {
        return Sparsity::dense(2, n_);
      } else {
        return f_.sparsity_in(i+1);
      }
    }
    virtual Sparsity get_sparsity_out(int i) {
      return f_.sparsity_out(i);
    }
    /// @}

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /** \brief  Evaluate numerically, work vectors given */
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Is the class able to propagate seeds through the algorithm? */
    virtual bool spCanEvaluate(bool fwd) { return fwd; }

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives */
    virtual Function get_forward(const std::string& name, int nfwd, Dict& opts);
    virtual int get_n_forward() const { return 64;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives */
    virtual Function get_reverse(const std::string& name, int nadj, Dict& opts);
    virtual int get_n_reverse() const { return 64;}
    ///@}

    /** \brief  Print description */
    virtual void print(std::ostream &stream) const;

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(CodeGenerator& g) const;

    /** \brief Generate code for the body of the C function */
    virtual void generateBody(CodeGenerator& g) const;

    // Default case;
    Function f_;

    std::pair<int, int> size_;

    double r_;

    int n_;

    int nnz_out_;

    /// Nonzero step for outputs
    std::vector<int> step_out_;

  };

} // namespace casadi
/// \endcond

#endif // CASADI_KERNEL_SUM_HPP
