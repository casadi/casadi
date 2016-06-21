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


#ifndef CASADI_SYMBOLIC_MX_HPP
#define CASADI_SYMBOLIC_MX_HPP

#include "mx_node.hpp"

/// \cond INTERNAL

namespace casadi {
  /** \brief Represents a symbolic MX
      \author Joel Andersson
      \date 2010
      A regular user is not supposed to work with this Node class.
      This user can call MX(name, n, m) directly.
  */
  class CASADI_EXPORT SymbolicMX : public MXNode {
  public:

    /** \brief  Constructors */
    explicit SymbolicMX(const std::string& name, int nrow=1, int ncol=1);

    /** \brief  Constructors */
    explicit SymbolicMX(const std::string& name, const Sparsity & sp);

    /// Destructor
    virtual ~SymbolicMX() {}

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /// Evaluate the function numerically
    virtual void eval(const double** arg, double** res, int* iw, double* w, int mem) const;

    /// Evaluate the function symbolically (SX)
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    /** \brief  Evaluate symbolically (MX) */
    virtual void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief  Propagate sparsity forward */
    virtual void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Get the name */
    virtual const std::string& name() const;

    /** \brief Get the operation */
    virtual int op() const { return OP_PARAMETER;}

    /** \brief  Check if valid function input */
    virtual bool is_valid_input() const { return true;}

    /** \brief Get the number of symbolic primitives */
    virtual int n_primitives() const { return 1;}

    /** \brief Get symbolic primitives */
    virtual void primitives(std::vector<MX>::iterator& it) const;

    /** \brief Split up an expression along symbolic primitives */
    virtual void split_primitives(const MX& x, std::vector<MX>::iterator& it) const;

    /** \brief Join an expression along symbolic primitives */
    virtual MX join_primitives(std::vector<MX>::const_iterator& it) const;

    /** \brief Detect duplicate symbolic expressions */
    virtual bool has_duplicates();

    /** \brief Reset the marker for an input expression */
    virtual void resetInput();

  protected:
    // Name of the variable
    std::string name_;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_SYMBOLIC_MX_HPP
