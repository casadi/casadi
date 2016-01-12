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


#ifndef CASADI_MAP_HPP
#define CASADI_MAP_HPP

#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** A map Base class for different map operations
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT MapBase : public FunctionInternal {
  public:
    // Create function (use instead of constructor)
    static MapBase* create(const std::string& name,
                           const Function& f, int n, const Dict& opts);

    /** \brief Constructor */
    MapBase(const std::string& name) : FunctionInternal(name) {}

    /** \brief Destructor */
    virtual ~MapBase();

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() const { return f_.n_in();}
    virtual size_t get_n_out() const { return f_.n_out();}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int ind) const {
      return repmat(f_.sparsity_in(ind), 1, n_);
    }
    virtual Sparsity get_sparsity_out(int ind) const {
      return repmat(f_.sparsity_out(ind), 1, n_);
    }
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    virtual std::vector<std::string> get_ischeme() const { return f_.name_in();}
    virtual std::vector<std::string> get_oscheme() const { return f_.name_out();}
    /// @}

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

  protected:
    // Constructor (protected, use create function above)
    MapBase(const std::string& name, const Function& f, int n);

    // The function which is to be evaluated in parallel
    Function f_;

    /// Number of Function inputs
    int n_in_;

    /// Number of Function outputs
    int n_out_;

    // Number of times to evaluate this function
    int n_;
  };

  /** A map Map for evaluating a function serially
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT MapSerial : public MapBase {
    friend class MapBase;
  protected:
    // Constructor (protected, use create function in MapBase)
    MapSerial(const std::string& name, const Function& f, int n)
      : MapBase(name, f, n) {}

    /** \brief  Destructor */
    virtual ~MapSerial();

    /** \brief  Evaluate or propagate sparsities */
    template<typename T>
    void evalGen(const T** arg, T** res, int* iw, T* w) const;

    /** \brief  Evaluate numerically, work vectors given */
    virtual void eval(Memory& mem, const double** arg, double** res, int* iw, double* w) const;

    /** \brief  evaluate symbolically while also propagating directional derivatives */
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Is the class able to propagate seeds through the algorithm? */
    virtual bool spCanEvaluate(bool fwd) { return true; }

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(CodeGenerator& g) const;

    /** \brief Generate code for the body of the C function */
    virtual void generateBody(CodeGenerator& g) const;

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

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

  };

#ifdef WITH_OPENMP
  /** A map Evaluate in parallel using OpenMP
      Inherits from MapSerial to allow fallback to serial methods
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT MapOmp : public MapSerial {
    friend class MapBase;
  protected:
    // Constructor (protected, use create function in MapBase)
    MapOmp(const std::string& name, const Function& f, int n) : MapSerial(name, f, n) {}

    /** \brief  Destructor */
    virtual ~MapOmp();

    /// Evaluate the function numerically
    virtual void eval(Memory& mem, const double** arg, double** res, int* iw, double* w) const;

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(CodeGenerator& g) const;

    /** \brief Generate code for the body of the C function */
    virtual void generateBody(CodeGenerator& g) const;

    /** \brief  Initialize */
    virtual void init(const Dict& opts);
  };
#endif // WITH_OPENMP

  /** A map operation that can also reduce certain arguments
      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT MapReduce : public MapBase {
  public:
    /** Types of parallelization supported */
    enum ParallelizationType {
      PARALLELIZATION_SERIAL,
      PARALLELIZATION_OMP
    };

    /** \brief Constructor (generic map) */
    MapReduce(const std::string& name, const Function& f, int n,
      const std::vector<bool> &repeat_in, const std::vector<bool> &repeat_out);

    /** \brief  Destructor */
    virtual ~MapReduce();

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /// Evaluate the function (template)
    template<typename T, typename R>
      void evalGen(const T** arg, T** res, int* iw, T* w, R reduction) const;

    /** \brief  Evaluate numerically, work vectors given */
    virtual void eval(Memory& mem, const double** arg, double** res, int* iw, double* w) const;

    /** \brief Quickfix to avoid segfault, #1552 */
    virtual bool canEvalSX() const {return true;}

    /** \brief  Evaluate symbolically, SXElem type, possibly nonmatching sparsity patterns */
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Is the class able to propagate seeds through the algorithm? */
    virtual bool spCanEvaluate(bool fwd) { return fwd; }

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int ind) const {
      return repeat_in_.at(ind) ? repmat(f_.sparsity_in(ind), 1, n_)
        : f_.sparsity_in(ind);
    }
    virtual Sparsity get_sparsity_out(int ind) const {
      return repeat_out_.at(ind) ? repmat(f_.sparsity_out(ind), 1, n_)
        : f_.sparsity_out(ind);
    }
    /// @}

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

    /// Indicate which inputs are repeated
    std::vector<bool> repeat_in_;

    /// Indicate which outputs are repeated
    std::vector<bool> repeat_out_;

    /// Nonzero step for inputs
    std::vector<int> step_in_;

    /// Nonzero step for outputs
    std::vector<int> step_out_;

    int nnz_out_;

    ParallelizationType parallelization_;

  };

} // namespace casadi
/// \endcond

#endif // CASADI_MAP_HPP
