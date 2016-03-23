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


#ifndef CASADI_MAPACCUM_HPP
#define CASADI_MAPACCUM_HPP

#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** MapAccum statement
      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT Mapaccum : public FunctionInternal {
    friend class MapAccum;
  public:

    // Create function (use instead of constructor)
    static Function create(const std::string& name,
                           Function& f, int n,
                           const std::vector<int>& accum_in, const std::vector<int>& accum_out,
                           const Dict& opts, bool reverse=false);

    /** \brief  Destructor */
    virtual ~Mapaccum();

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() { return f_.n_in();}
    virtual size_t get_n_out() { return f_.n_out();}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int i) {
      if (i<n_accum_) {
        return f_.sparsity_in(i);
      } else {
        return repmat(f_.sparsity_in(i), 1, n_);
      }
    }
    virtual Sparsity get_sparsity_out(int i) {
      return repmat(f_.sparsity_out(i), 1, n_);
    }
    /// @}

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /// Evaluate the function (template)
    template<typename T, typename R>
    void evalGen(const T** arg, T** res, int* iw, T* w, R reduction) const;

    /** \brief Binary or, helper function */
    static inline bvec_t orop(bvec_t x, bvec_t y) { return x | y; }

    /** \brief  Evaluate numerically, work vectors given */
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /** \brief Quickfix to avoid segfault, #1552 */
    virtual bool canEvalSX() const {return true;}

    /** \brief  Evaluate symbolically, SXElem type, possibly nonmatching sparsity patterns */
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Is the class able to propagate seeds through the algorithm? */
    virtual bool spCanEvaluate(bool fwd) { return true; }

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

    int n_;

    int nnz_out_;

    /// Nonzero step for inputs
    std::vector<int> step_in_;

    /// Nonzero step for outputs
    std::vector<int> step_out_;

    /// Number of accumulated inputs/outputs
    int n_accum_;

    /// Indicates the order of accumulation
    bool reverse_;

    /// Total number of accumulator nonzeros
    int nnz_accum_;

  protected:
    /** \brief Constructor (generic mapaccum) */
    Mapaccum(const std::string& name, const Function& f, int n, int n_accum, bool reverse);

  };

} // namespace casadi
/// \endcond

#endif // CASADI_MAPACCUM_HPP
