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

  /** Map statement
      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT MapInternal : public FunctionInternal {
    friend class Map;
  public:

    enum ParallelizationType {PARALLELIZATION_SERIAL, PARALLELIZATION_OMP};

    /** \brief Constructor (generic map) */
    MapInternal(const Function& f, int n,
      const std::vector<bool> &repeat_in, const std::vector<bool> &repeat_out);

    /** \brief  clone function */
    virtual MapInternal* clone() const { return new MapInternal(*this);}

    /** \brief  Destructor */
    virtual ~MapInternal();

    /** \brief  Initialize */
    virtual void init();

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T** arg, T** res, int* iw, T* w);

    /** \brief  Evaluate numerically, work vectors given */
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /** \brief  Evaluate symbolically, SXElement type, possibly nonmatching sparsity patterns */
    virtual void evalSX(const SXElement** arg, SXElement** res,
                                int* iw, SXElement* w);

    /** \brief  Evaluate symbolically, MX type */
    //virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res);

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives */
    virtual Function getDerForward(const std::string& name, int nfwd, const Dict& opts);
    virtual int numDerForward() const { return 64;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives */
    virtual Function getDerReverse(const std::string& name, int nadj, const Dict& opts);
    virtual int numDerReverse() const { return 64;}
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
