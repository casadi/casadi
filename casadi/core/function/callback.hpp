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
    /** \brief Default constructor */
    Callback();

    /** \brief Copy constructor (throws an error) */
    Callback(const Callback& obj);

    /** \brief Create an owning reference, given a pointer to a derived object */
    static Function create(const std::string& name, Callback* n,
                           const Dict& opts=Dict());

    /** \brief Construct internal object
     * This is the step that actually construct the internal object, as the
     * class constructor only creates a null pointer.
     * It should be called from the user constructor.
     */
    void construct(const std::string& name, const Dict& opts=Dict());

#ifndef SWIG
    /** \brief Transfer ownership to the internal class
     * With a call to this function, the public class will be owned by the
     * internal class.
     * For this to work, the object must have been created with "new" (and must
     * not be deleted with "delete" as this is handled internally.
     * There also has to be at least one owning reference to the class.
     */
    void transfer_ownership();
#endif // SWIG

    /** \brief  Destructor */
    virtual ~Callback();

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

    /** \brief Evaluate numerically, temporary matrices and work vectors */
    virtual std::vector<DM> eval(const std::vector<DM>& arg);

#ifndef SWIG
    /** \brief Evaluate numerically, work vectors given */
    virtual void eval(const double** arg, double** res, int* iw, double* w, int mem);
#endif // SWIG

   /** \brief Get the number of inputs
     * This function is called during construction.
     */
    virtual int get_n_in();

    /** \brief Get the number of outputs
     * This function is called during construction.
     */
    virtual int get_n_out();

    /** \brief Get the sparsity of an input
     * This function is called during construction.
     */
    virtual Sparsity get_sparsity_in(int i);

    /** \brief Get the sparsity of an output
     * This function is called during construction.
     */
    virtual Sparsity get_sparsity_out(int i);

    /** \brief Get the sparsity of an input
     * This function is called during construction.
     */
    virtual std::string get_name_in(int i);

    /** \brief Get the sparsity of an output
     * This function is called during construction.
     */
    virtual std::string get_name_out(int i);

    ///@{
    /** \brief Return Jacobian of all input elements with respect to all output elements */
    virtual bool has_jacobian() const;
    virtual Function get_jacobian(const std::string& name, const Dict& opts);
    ///@}

    /** \brief [DEPRECATED] Overload get_forward_new instead */
    virtual Function get_forward(const std::string& name, int nfwd, Dict& opts);

    ///@{
    /** \brief Return function that calculates forward derivatives
     *    forward(nfwd) returns a cached instance if available,
     *    and calls <tt>Function get_forward(int nfwd)</tt>
     *    if no cached version is available.
     */
    virtual Function get_forward_new(const std::string& name, int nfwd,
                                     const std::vector<std::string>& i_names,
                                     const std::vector<std::string>& o_names,
                                     const Dict& opts);
    virtual int get_n_forward() const;
    ///@}

    /** \brief [DEPRECATED] Overload get_reverse_new instead */
    virtual Function get_reverse(const std::string& name, int nadj, Dict& opts);

    ///@{
    /** \brief Return function that calculates adjoint derivatives
     *    reverse(nadj) returns a cached instance if available,
     *    and calls <tt>Function get_reverse(int nadj)</tt>
     *    if no cached version is available.
     */
    virtual Function get_reverse_new(const std::string& name, int nadj,
                                     const std::vector<std::string>& i_names,
                                     const std::vector<std::string>& o_names,
                                     const Dict& opts);
    virtual int get_n_reverse() const;
    ///@}

#ifndef SWIG
    private:
    /** \brief  Access functions of the node */
    CallbackInternal* operator->();

    /** \brief  Const access functions of the node */
    const CallbackInternal* operator->() const;
#endif // SWIG
  };

} // namespace casadi

#endif // CASADI_CALLBACK_HPP
