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
   */
  class CASADI_EXPORT Callback : public Function {
  public:
    /** \brief Get type name */
    static std::string type_name() {return "Callback";}

    /** \brief Default constructor */
    Callback();

    /** \brief Copy constructor (throws an error) */
    Callback(const Callback& obj);

    /** \brief  Destructor */
    virtual ~Callback();

    /** \brief Construct internal object
     * This is the step that actually construct the internal object, as the
     * class constructor only creates a null pointer.
     * It should be called from the user constructor.
     */
    void construct(const std::string& name, const Dict& opts=Dict());

    /** \brief Initialize the object
     * This function is called after the object construction (for the whole class
     * hierarchy) is complete, but before the finalization step.
     * It is called recursively for the whole class hierarchy, starting with the
     * lowest level.
     */
    virtual void init() {}

    /** \brief Finalize the object
     * This function is called after the construction and init steps are completed,
     * but before user functions are called.
     * It is called recursively for the whole class hierarchy, starting with the
     * highest level.
    */
    virtual void finalize() {}

    /** \brief Evaluate numerically, using temporary matrices and work vectors
     * 
     * This signature is not thread-safe.
     * For guaranteed thread-safety, use `eval_buffer`
     */
    virtual std::vector<DM> eval(const std::vector<DM>& arg) const;

    /** \brief A copy-free low level interface 
     * 
     * In Python, you will be passed two tuples of memoryview objects
    */
    virtual int eval_buffer(const double **arg, const std::vector<casadi_int>& sizes_arg,
                              double **res, const std::vector<casadi_int>& sizes_res) const;
    virtual bool has_eval_buffer() const;

    /** \brief Get the number of inputs
     * This function is called during construction.
     */
    virtual casadi_int get_n_in();

    /** \brief Get the number of outputs
     * This function is called during construction.
     */
    virtual casadi_int get_n_out();

    /** \brief Get the sparsity of an input
     * This function is called during construction.
     */
    virtual Sparsity get_sparsity_in(casadi_int i);

    /** \brief Get the sparsity of an output
     * This function is called during construction.
     */
    virtual Sparsity get_sparsity_out(casadi_int i);

    /** \brief Get the sparsity of an input
     * This function is called during construction.
     */
    virtual std::string get_name_in(casadi_int i);

    /** \brief Get the sparsity of an output
     * This function is called during construction.
     */
    virtual std::string get_name_out(casadi_int i);

    /** \brief Do the derivative functions need nondifferentiated outputs? */
    virtual bool uses_output() const;

    ///@{
    /** \brief Return Jacobian of all input elements with respect to all output elements */
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
     */
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
     */
    virtual bool has_reverse(casadi_int nadj) const;
    virtual Function get_reverse(casadi_int nadj, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const;
    ///@}

    ///@{
    /** \brief Return sparsity of Jacobian of all input elements
     * with respect to all output elements */
    virtual bool has_jacobian_sparsity() const { return false; }
    virtual Sparsity get_jacobian_sparsity() const { return Sparsity(); }
    ///@}
  };

} // namespace casadi

#endif // CASADI_CALLBACK_HPP
