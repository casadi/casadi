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


#ifndef CASADI_MAPACCUM_INTERNAL_HPP
#define CASADI_MAPACCUM_INTERNAL_HPP

#include "mapaccum.hpp"
#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** MapAccum statement
      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT MapAccumInternal : public FunctionInternal {
    friend class MapAccum;
  public:

    /** \brief Constructor (generic mapaccum) */
    MapAccumInternal(const Function& f, int n,
      const std::vector<bool>& input_accum,
      const std::vector<int>& output_accum,
      bool reverse);

    /** \brief  clone function */
    virtual MapAccumInternal* clone() const { return new MapAccumInternal(*this);}

    /** \brief  Destructor */
    virtual ~MapAccumInternal();

    /** \brief  Initialize */
    virtual void init();

    /// Evaluate the function (template)
    template<typename T, typename R>
    void evalGen(const T** arg, T** res, int* iw, T* w, R reduction);

    /** \brief Binary or, helper function */
    static inline bvec_t orop(bvec_t x, bvec_t y) { return x | y; }

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

    // Default case;
    Function f_;

    int n_;

    int nnz_out_;

    /// Nonzero step for inputs
    std::vector<int> step_in_;

    /// Nonzero step for outputs
    std::vector<int> step_out_;

    std::vector<bool> input_accum_;
    std::vector<int> output_accum_;

    /// Indicates the order of accumulation
    bool reverse_;

    /// Total number of accumulator nonzeros
    int nnz_accum_;

  };

} // namespace casadi
/// \endcond

#endif // CASADI_MAPACCUM_INTERNAL_HPP
