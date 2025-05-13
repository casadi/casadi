/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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

#ifndef CASADI_CALLBACK_HPP
#define CASADI_CALLBACK_HPP

#include "function.hpp"

namespace casadi {
  /** Forward declaration of internal class */
  class CallbackInternal;

  /** \brief Callback function functionality

   This class provides a public API to the FunctionInternal class that
   can be subclassed by the user, who is then able to implement the different
   virtual method.
   Note that the Function class also provides a public API to FunctionInternal,
   but only allows calling, not being called.

   The user is responsible for not deleting this class for the lifetime
   of the internal function object.

   \author Joris Gillis, Joel Andersson
   \date 2015

      \identifier{o0} */
  class CASADI_EXPORT Callback : public Function {
  public:
    /** \brief Get type name

        \identifier{o1} */
    static std::string type_name() {return "Callback";}

    /** \brief Default constructor

        \identifier{o2} */
    Callback();

    /** \brief Copy constructor (throws an error)

        \identifier{o3} */
    Callback(const Callback& obj);

    /** \brief  Destructor

        \identifier{o4} */
    virtual ~Callback();

    /** \brief Construct internal object

     * This is the step that actually construct the internal object, as the
     * class constructor only creates a null pointer.
     * It should be called from the user constructor.

        \identifier{o5} */
    void construct(const std::string& name, const Dict& opts=Dict());

    /** \brief Initialize the object

     * This function is called after the object construction (for the whole class
     * hierarchy) is complete, but before the finalization step.
     * It is called recursively for the whole class hierarchy, starting with the
     * lowest level.

        \identifier{o6} */
    virtual void init() {}

    /** \brief Finalize the object

     * This function is called after the construction and init steps are completed,
     * but before user functions are called.
     * It is called recursively for the whole class hierarchy, starting with the
     * highest level.

        \identifier{o7} */
    virtual void finalize() {}

    /** \brief Evaluate numerically, using temporary matrices and work vectors
     *
     * This signature is not thread-safe.
     * For guaranteed thread-safety, use `eval_buffer`

        \identifier{o8} */
    virtual std::vector<DM> eval(const std::vector<DM>& arg) const;

    /** \brief A copy-free low level interface
     *
     * In Python, you will be passed two tuples of memoryview objects
     * Note that only the structural nonzeros are present in the memoryview objects/buffers.
     *
     * Make sure to override has_eval_buffer() to indicate support for this method.

        \identifier{o9} */
    virtual int eval_buffer(const double **arg, const std::vector<casadi_int>& sizes_arg,
                              double **res, const std::vector<casadi_int>& sizes_res) const;
    /** \brief Does the Callback class support a copy-free low level interface ?
     *

        \identifier{265} */
    virtual bool has_eval_buffer() const;

    /** \brief Get the number of inputs

     * This function is called during construction.

        \identifier{oa} */
    virtual casadi_int get_n_in();

    /** \brief Get the number of outputs

     * This function is called during construction.

        \identifier{ob} */
    virtual casadi_int get_n_out();

    /** \brief Get the sparsity of an input

     * This function is called during construction.

        \identifier{oc} */
    virtual Sparsity get_sparsity_in(casadi_int i);

    /** \brief Get the sparsity of an output

     * This function is called during construction.

        \identifier{od} */
    virtual Sparsity get_sparsity_out(casadi_int i);

    /** \brief Get the name of an input

     * This function is called during construction.

        \identifier{oe} */
    virtual std::string get_name_in(casadi_int i);

    /** \brief Get the name of an output

     * This function is called during construction.

        \identifier{of} */
    virtual std::string get_name_out(casadi_int i);

    /** \brief Do the derivative functions need nondifferentiated outputs?

        \identifier{og} */
    virtual bool uses_output() const;

    /** \brief Customize calls to the function factory

        \identifier{2de} */
    virtual Function get_factory(const std::string& name,
        const std::vector<std::string>& s_in,
        const std::vector<std::string>& s_out,
        const Function::AuxOut& aux,
        const Dict& opts) const;

    ///@{
    /** \brief Return Jacobian of all input elements with respect to all output elements

        \identifier{oh} */
    virtual bool has_jacobian() const;
    virtual Function get_jacobian(const std::string& name,
                                  const std::vector<std::string>& inames,
                                  const std::vector<std::string>& onames,
                                  const Dict& opts) const;
    ///@}

    ///@{
    /** \brief Return function that calculates forward derivatives

     *    forward(nfwd) returns a cached instance if available,
     *    and calls <tt>Function get_forward(casadi_int nfwd)</tt>
     *    if no cached version is available.

        \identifier{oi} */
    virtual bool has_forward(casadi_int nfwd) const;
    virtual Function get_forward(casadi_int nfwd, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const;
    ///@}

    ///@{
    /** \brief Return function that calculates adjoint derivatives

     *    reverse(nadj) returns a cached instance if available,
     *    and calls <tt>Function get_reverse(casadi_int nadj)</tt>
     *    if no cached version is available.

        \identifier{oj} */
    virtual bool has_reverse(casadi_int nadj) const;
    virtual Function get_reverse(casadi_int nadj, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const;
    ///@}

    ///@{
    /** \brief Return sparsity of Jacobian of all input elements

     * with respect to all output elements

        \identifier{ok} */
    virtual bool has_jac_sparsity(casadi_int oind, casadi_int iind) const { return false; }
    virtual Sparsity get_jac_sparsity(casadi_int oind, casadi_int iind, bool symmetric) const {
      return Sparsity(); }
    ///@}
  };

} // namespace casadi

#endif // CASADI_CALLBACK_HPP
