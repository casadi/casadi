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

  /** Evaluate in parallel
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT Map : public FunctionInternal {
  public:
    // Create function (use instead of constructor)
    static Function create(const std::string& name,
                          const std::string& parallelization, Function& f, int n,
                          const Dict& opts);

    /** \brief Destructor */
    virtual ~Map();

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int i) {
      return repmat(f_.sparsity_in(i), 1, n_);
    }
    virtual Sparsity get_sparsity_out(int i) {
      return repmat(f_.sparsity_out(i), 1, n_);
    }
    /// @}

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

    /** \brief  Evaluate or propagate sparsities */
    template<typename T>
    void evalGen(const T** arg, T** res, int* iw, T* w) const;

    /// Evaluate the function numerically
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /// Type of parallellization
    virtual std::string parallelization() const { return "serial"; }

    /** \brief  evaluate symbolically while also propagating directional derivatives */
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    /** \brief  Propagate sparsity forward */
    virtual void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    ///@{
    /// Is the class able to propagate seeds through the algorithm?
    virtual bool has_spfwd() const { return true;}
    virtual bool has_sprev() const { return true;}
    ///@}

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(CodeGenerator& g) const;

    /** \brief Is codegen supported? */
    virtual bool has_codegen() const { return true;}

    /** \brief Generate code for the body of the C function */
    virtual void generateBody(CodeGenerator& g) const;

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives */
    virtual Function get_forward(const std::string& name, int nfwd,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names,
                                 const Dict& opts);
    virtual int get_n_forward() const { return 64;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives */
    virtual Function get_reverse(const std::string& name, int nadj,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names,
                                 const Dict& opts);
    virtual int get_n_reverse() const { return 64;}
    ///@}

  protected:
    // Constructor (protected, use create function)
    Map(const std::string& name, const Function& f, int n);

    // The function which is to be evaluated in parallel
    Function f_;

    // Number of times to evaluate this function
    int n_;
  };

  /** A map Evaluate in parallel using OpenMP
      Note: Do not use this class with much more than the intended number of
      threads for the parallel evaluation as it will cause excessive memory use.

      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT MapOmp : public Map {
    friend class Map;
  protected:
    // Constructor (protected, use create function in Map)
    MapOmp(const std::string& name, const Function& f, int n) : Map(name, f, n) {}

    /** \brief  Destructor */
    virtual ~MapOmp();

    /// Evaluate the function numerically
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /// Type of parallellization
    virtual std::string parallelization() const { return "openmp"; }

    /** \brief Generate code for the body of the C function */
    virtual void generateBody(CodeGenerator& g) const;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_MAP_HPP
