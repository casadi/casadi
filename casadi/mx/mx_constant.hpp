/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef MX_CONSTANT_HPP
#define MX_CONSTANT_HPP

#include "mx_node.hpp"

namespace CasADi{

/** \brief Represents an MX that is only composed of a constant.
	\author Joel Andersson 
	\date 2010

	A regular user is not supposed to work with this Node class.
	This user can call MX(double) directly, or even rely on implicit typecasting.
	\sa zeros , ones
*/
class MXConstant : public MXNode{
  public:

    /** \brief  Constructor */
    MXConstant(const Matrix<double> &x);

    /** \brief  Clone function */
    virtual MXConstant* clone() const;

    /** \brief  Print */
    virtual void print(std::ostream &stream, const std::vector<std::string>& args) const;

    /** \brief  Evaluate the function and store the result in the node */
    virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrV& fwdSens, const DMatrixPtrV& adjSeed, DMatrixPtrVV& adjSens, int nfwd, int nadj);

    /** \brief  Indicate that the node is constant */
    virtual bool isConstant() const;
    
    /// Symbolic evaluation (matrix graph)
    virtual MX eval(const std::vector<MX>& x){return MX::create(this);}
    
    /** \brief  Evaluate symbolically (SX) */
    virtual void evaluateSX(const std::vector<SXMatrix*> &input, SXMatrix& output);

    /** \brief  data member */
    Matrix<double> x_;

};

} // namespace CasADi


#endif // MX_CONSTANT_HPP
