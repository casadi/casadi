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

    /** \brief  Clone function */
    virtual SymbolicMX* clone() const;

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /// Evaluate the function numerically
    virtual void evalD(cp_double* input, p_double* output, int* itmp, double* rtmp);

    /// Evaluate the function symbolically (SX)
    virtual void evalSX(cp_SXElement* input, p_SXElement* output,
                            int* itmp, SXElement* rtmp);

    /** \brief  Evaluate symbolically (MX) */
    virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(cp_bvec_t* arg,
                       p_bvec_t* res, int* itmp, bvec_t* rtmp);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(p_bvec_t* arg,
                       p_bvec_t* res, int* itmp, bvec_t* rtmp);

    /** \brief  Get the name */
    virtual const std::string& getName() const;

    /** \brief Get the operation */
    virtual int getOp() const { return OP_PARAMETER;}

    /** \brief  Check if valid function input */
    virtual bool isValidInput() const { return true;}

    /** \brief Get the number of symbolic primitives */
    virtual int numPrimitives() const { return 1;}

    /** \brief Get symbolic primitives */
    virtual void getPrimitives(std::vector<MX>::iterator& it) const;

    /** \brief Split up an expression along symbolic primitives */
    virtual void splitPrimitives(const MX& x, std::vector<MX>::iterator& it) const;

    /** \brief Join an expression along symbolic primitives */
    virtual MX joinPrimitives(std::vector<MX>::const_iterator& it) const;

    /** \brief Detect duplicate symbolic expressions */
    virtual bool hasDuplicates();

    /** \brief Reset the marker for an input expression */
    virtual void resetInput();

  protected:
    // Name of the variable
    std::string name_;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_SYMBOLIC_MX_HPP
