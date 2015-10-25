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


#ifndef CASADI_MAP_INTERNAL_HPP
#define CASADI_MAP_INTERNAL_HPP

#include "map.hpp"
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
    static MapBase* create(const Function& f, int n, const Dict& opts);

    /** \brief  Destructor */
    virtual ~MapBase();

    /** \brief  Initialize */
    virtual void init();

  protected:
    // Constructor (protected, use create function above)
    MapBase(const Function& f, int n);

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
    MapSerial(const Function& f, int n) : MapBase(f, n) {}

    /** \brief  Destructor */
    virtual ~MapSerial();

    /** \brief  clone function */
    virtual MapSerial* clone() const { return new MapSerial(*this);}

    /** \brief  Evaluate or propagate sparsities */
    template<typename T>
    void evalGen(const T** arg, T** res, int* iw, T* w);

    /** \brief  Evaluate numerically, work vectors given */
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /** \brief  evaluate symbolically while also propagating directional derivatives */
    virtual void evalSX(const SXElement** arg, SXElement** res,
                        int* iw, SXElement* w);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Initialize */
    virtual void init();
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
    MapOmp(const Function& f, int n) : MapSerial(f, n) {}

    /** \brief  clone function */
    virtual MapOmp* clone() const { return new MapOmp(*this);}

    /** \brief  Destructor */
    virtual ~MapOmp();

    /// Evaluate the function numerically
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /** \brief  Initialize */
    virtual void init();
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
    MapReduce(const Function& f, int n,
      const std::vector<bool> &repeat_in, const std::vector<bool> &repeat_out);

    /** \brief  clone function */
    virtual MapReduce* clone() const { return new MapReduce(*this);}

    /** \brief  Destructor */
    virtual ~MapReduce();

    /** \brief  Initialize */
    virtual void init();

    /// Evaluate the function (template)
    template<typename T, typename R>
      void evalGen(const T** arg, T** res, int* iw, T* w,
                   void (FunctionInternal::*ptrEval)(const T** arg, T** res, int* iw, T* w),
                   R reduction);

    /** \brief  Evaluate numerically, work vectors given */
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /** \brief Quickfix to avoid segfault, #1552 */
    virtual bool canEvalSX() const {return true;}

    /** \brief  Evaluate symbolically, SXElement type, possibly nonmatching sparsity patterns */
    virtual void evalSX(const SXElement** arg, SXElement** res,
                                int* iw, SXElement* w);

    /** \brief  Evaluate symbolically, MX type */
    //virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res);

    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Is the class able to propagate seeds through the algorithm? */
    virtual bool spCanEvaluate(bool fwd) { return true; }

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives */
    virtual Function getDerForward(const std::string& name, int nfwd, Dict& opts);
    virtual int numDerForward() const { return 64;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives */
    virtual Function getDerReverse(const std::string& name, int nadj, Dict& opts);
    virtual int numDerReverse() const { return 64;}
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

#endif // CASADI_MAP_INTERNAL_HPP
