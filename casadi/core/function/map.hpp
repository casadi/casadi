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
    static Function create(const std::string& name,
                          const std::string& parallelization, Function& f, int n,
                          const std::vector<int>& reduce_in, const std::vector<int>& reduce_out,
                          const Dict& opts);

    // Create function (use instead of constructor)
    static MapBase* create(const std::string& name,
                          const std::string& parallelization, const Function& f, int n,
                          const std::vector<int>& reduce_in, const std::vector<int>& reduce_out);

    /** \brief Constructor */
    MapBase(const std::string& name) : FunctionInternal(name) {}

    /** \brief Destructor */
    virtual ~MapBase();

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() { return f_.n_in();}
    virtual size_t get_n_out() { return f_.n_out();}
    ///@}

    ///@{
    /** \brief Names of function input and outputs */
    virtual std::string get_name_in(int i) { return f_.name_in(i);}
    virtual std::string get_name_out(int i) { return f_.name_out(i);}
    /// @}

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /// Type of parallellization
    virtual std::string parallelization() const=0;

  protected:
    // Constructor (protected, use create function above)
    MapBase(const std::string& name, const Function& f, int n);

    /// Propagate optiosn to derivatives
    void propagate_options(Dict& opts);

    // The function which is to be evaluated in parallel
    Function f_;

    /// Number of Function inputs
    int n_in_;

    /// Number of Function outputs
    int n_out_;

    // Number of times to evaluate this function
    int n_;

    // Number of threads
    int n_threads_;
  };

  /** A map Base class for pure maps (no reduced in/out)
      \author Joris Gillis
      \date 2016
  */
  class CASADI_EXPORT PureMap : public MapBase {
    friend class MapBase;
  protected:
    // Constructor (protected, use create function in MapBase)
    PureMap(const std::string& name, const Function& f, int n)
      : MapBase(name, f, n) {}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int i) {
      return repmat(f_.sparsity_in(i), 1, n_);
    }
    virtual Sparsity get_sparsity_out(int i) {
      return repmat(f_.sparsity_out(i), 1, n_);
    }
    /// @}

    /** \brief  Evaluate or propagate sparsities */
    template<typename T>
    void evalGen(const T** arg, T** res, int* iw, T* w) const;

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


  /** A map operation that can also reduce certain arguments
      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT MapSum : public MapBase {
    friend class MapBase;
  public:

    /** \brief Constructor (generic map) */
    MapSum(const std::string& name, const Function& f, int n,
      const std::vector<bool> &repeat_in, const std::vector<bool> &repeat_out);

    /** \brief  Destructor */
    virtual ~MapSum();

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /// Evaluate the function (template)
    template<typename T, typename R>
      void evalGen(const T** arg, T** res, int* iw, T* w, R reduction) const;

    /** \brief  Evaluate numerically, work vectors given */
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /** \brief Quickfix to avoid segfault, #1552 */
    virtual bool canEvalSX() const {return true;}

    /** \brief  Evaluate symbolically, SXElem type, possibly nonmatching sparsity patterns */
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Is the class able to propagate seeds through the algorithm? */
    virtual bool spCanEvaluate(bool fwd) { return true; }

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int i) {
      return repeat_in_.at(i) ? repmat(f_.sparsity_in(i), 1, n_)
        : f_.sparsity_in(i);
    }
    virtual Sparsity get_sparsity_out(int i) {
      return repeat_out_.at(i) ? repmat(f_.sparsity_out(i), 1, n_)
        : f_.sparsity_out(i);
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

    int nnz_out_;

    /// Nonzero step for inputs
    std::vector<int> step_in_;

    /// Nonzero step for outputs
    std::vector<int> step_out_;

  };

  /** A map Map for evaluating a function serially
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT MapSerial : public PureMap {
    friend class MapBase;
    friend class PureMap;
  protected:
    // Constructor (protected, use create function in MapBase)
    MapSerial(const std::string& name, const Function& f, int n)
      : PureMap(name, f, n) {}

    /** \brief  Destructor */
    virtual ~MapSerial();

    /// Evaluate the function numerically
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /// Type of parallellization
    virtual std::string parallelization() const { return "serial"; }

  };

  /** A mapsum evaluation in serial
      \author Joris Gillis
      \date 2016
  */
  class CASADI_EXPORT MapSumSerial : public MapSum {
    friend class MapSum;
    friend class MapBase;
  protected:
    // Constructor (protected, use create function in MapBase)
    MapSumSerial(const std::string& name, const Function& f, int n,
      const std::vector<bool> &repeat_in, const std::vector<bool> &repeat_out)
      : MapSum(name, f, n, repeat_in, repeat_out) {}

    /// Evaluate the function numerically
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /** \brief  Destructor */
    virtual ~MapSumSerial();

    /// Type of parallellization
    virtual std::string parallelization() const { return "serial"; }

  };

#ifdef WITH_OPENMP
  /** A map Evaluate in parallel using OpenMP
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT MapOmp : public PureMap {
    friend class PureMap;
    friend class MapBase;
  protected:
    // Constructor (protected, use create function in MapBase)
    MapOmp(const std::string& name, const Function& f, int n) : PureMap(name, f, n) {}

    /** \brief  Destructor */
    virtual ~MapOmp();

    /// Evaluate the function numerically
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /// Type of parallellization
    virtual std::string parallelization() const { return "openmp"; }

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(CodeGenerator& g) const;

    /** \brief Generate code for the body of the C function */
    virtual void generateBody(CodeGenerator& g) const;

  };

  /** A mapsum Evaluate in parallel using OpenMP
      \author Joris Gillis
      \date 2016
  */
  class CASADI_EXPORT MapSumOmp : public MapSum {
    friend class MapSum;
    friend class MapBase;
  protected:
    // Constructor (protected, use create function in MapBase)
    MapSumOmp(const std::string& name, const Function& f, int n,
      const std::vector<bool> &repeat_in, const std::vector<bool> &repeat_out)
      : MapSum(name, f, n, repeat_in, repeat_out) {}

    /** \brief  Destructor */
    virtual ~MapSumOmp();

    /// Evaluate the function numerically
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /// Type of parallellization
    virtual std::string parallelization() const { return "openmp"; }

  };
#endif // WITH_OPENMP


} // namespace casadi
/// \endcond

#endif // CASADI_MAP_HPP
