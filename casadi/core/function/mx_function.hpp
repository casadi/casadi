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


#ifndef CASADI_MX_FUNCTION_HPP
#define CASADI_MX_FUNCTION_HPP

#include <set>
#include <map>
#include <vector>
#include <iostream>

#include "x_function.hpp"
#include "../mx/mx_node.hpp"

/// \cond INTERNAL

namespace casadi {

#ifndef SWIG
  /** \brief  An element of the algorithm, namely an MX node */
  struct MXAlgEl {
    /// Operator index
    int op;

    /// Data associated with the operation
    MX data;

    /// Work vector indices of the arguments
    std::vector<int> arg;

    /// Work vector indices of the results
    std::vector<int> res;
  };
#endif // SWIG

  /** \brief  Internal node class for MXFunction
      \author Joel Andersson
      \date 2010-2015
  */
  class CASADI_EXPORT MXFunction :
        public XFunction<MXFunction, MX, MXNode>{
  public:
    /** \brief  An element of the algorithm, namely an MX node */
    typedef MXAlgEl AlgEl;

    /** \brief  All the runtime elements in the order of evaluation */
    std::vector<AlgEl> algorithm_;

    /** \brief Offsets for elements in the w_ vector */
    std::vector<int> workloc_;

    /// Free variables
    std::vector<MX> free_vars_;

    /// Default input values
    std::vector<double> default_in_;

    /** \brief Constructor */
    MXFunction(const std::string& name, const std::vector<MX>& input,
                       const std::vector<MX>& output);

    /** \brief  Destructor */
    virtual ~MXFunction();

    /** \brief  Evaluate numerically, work vectors given */
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /** \brief  Print description */
    virtual void print(std::ostream &stream) const;

    /** \brief Get type name */
    virtual std::string type_name() const;

    /** \brief Check if the function is of a particular type */
    virtual bool is_a(const std::string& type, bool recursive) const;

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(CodeGenerator& g) const;

    /** \brief Codegen incref for dependencies */
    virtual void codegen_incref(CodeGenerator& g) const;

    /** \brief Codegen decref for dependencies */
    virtual void codegen_decref(CodeGenerator& g) const;

    /** \brief Generate code for the body of the C function */
    virtual void generateBody(CodeGenerator& g) const;

    /** \brief Extract the residual function G and the modified function Z out of an expression
     * (see Albersmeyer2010 paper) */
    void generate_lifted(Function& vdef_fcn, Function& vinit_fcn);

    /** \brief Generate a function that calculates a Jacobian function by operator overloading */
    virtual Function getNumericJacobian(const std::string& name, int iind, int oind,
                                        bool compact, bool symmetric, const Dict& opts);

    /** \brief Evaluate symbolically, SX type*/
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    /** \brief Evaluate symbolically, MX type */
    virtual void eval_mx(const MXVector& arg, MXVector& res, bool always_inline, bool never_inline);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fwdSeed,
                        std::vector<std::vector<MX> >& fwdSens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& adjSeed,
                        std::vector<std::vector<MX> >& adjSens);

    /** \brief Create call to (cached) derivative function, forward mode  */
    virtual void forward_mx(const std::vector<MX>& arg, const std::vector<MX>& res,
                            const std::vector<std::vector<MX> >& fseed,
                            std::vector<std::vector<MX> >& fsens,
                            bool always_inline, bool never_inline);

    /** \brief Create call to (cached) derivative function, reverse mode  */
    virtual void reverse_mx(const std::vector<MX>& arg, const std::vector<MX>& res,
                            const std::vector<std::vector<MX> >& aseed,
                            std::vector<std::vector<MX> >& asens,
                            bool always_inline, bool never_inline);

    /** \brief Expand the matrix valued graph into a scalar valued graph */
    Function expand(const std::vector<SX>& inputv);

    /// Get a vector of symbolic variables corresponding to the outputs
    virtual std::vector<MX> symbolicOutput(const std::vector<MX>& arg);

    /** \brief Create call */
    virtual std::vector<MX> create_call(const std::vector<MX>& arg);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /// Is the class able to propagate seeds through the algorithm?
    virtual bool spCanEvaluate(bool fwd) { return true;}

    // print an element of an algorithm
    void print(std::ostream &stream, const AlgEl& el) const;

    /** \brief Gradient expression */
    virtual MX grad_mx(int iind=0, int oind=0);

    /** \brief Tangent expression */
    virtual MX tang_mx(int iind=0, int oind=0);

    /** \brief Jacobian expression */
    virtual MX jac_mx(int iind=0, int oind=0, bool compact=false, bool symmetric=false,
                      bool always_inline=true, bool never_inline=false);

    ///@{
    /** \brief Get function input(s) and output(s)  */
    virtual const MX mx_in(int ind) const;
    virtual const std::vector<MX> mx_in() const;
    ///@}

    /// Get free variables (MX)
    virtual std::vector<MX> free_mx() const {return free_vars_;}

    /** \brief Does the function have free variables */
    virtual bool has_free() const { return !free_vars_.empty();}

    /** \brief Number of nodes in the algorithm */
    virtual int n_nodes() const { return algorithm_.size();}

    /** \brief Get default input value */
    virtual double default_in(int ind) const { return default_in_.at(ind);}
  };

} // namespace casadi
/// \endcond

#endif // CASADI_MX_FUNCTION_HPP

