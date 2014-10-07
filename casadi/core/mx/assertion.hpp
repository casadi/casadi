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


#ifndef CASADI_ASSERTION_HPP
#define CASADI_ASSERTION_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL
namespace casadi {
  /** \brief Assertion
      \author Joris Gillis
      \date 2013
  */
  class CASADI_CORE_EXPORT Assertion : public MXNode {
  public:

    /// Constructor
    Assertion(const MX& x, const MX& y, const std::string & s);

    /// Clone function
    virtual Assertion* clone() const { return new Assertion(*this);}

    /// Destructor
    virtual ~Assertion() {}

    /// Evaluate the function symbolically (MX)
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                            MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                            bool output_given);

    /// Evaluate the function numerically
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                           std::vector<double>& rtmp);

    /// Evaluate the function symbolically (SX)
    virtual void evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                            std::vector<SXElement>& rtmp);

    /// Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief Get the operation */
    virtual int getOp() const { return OP_ASSERTION;}

    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);


  private:
    std::string fail_message_;
  };


} // namespace casadi

/// \endcond

#endif // CASADI_ASSERTION_HPP
