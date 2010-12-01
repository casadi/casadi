/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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
  MXConstant(const double *x, int n, int m, char order ='R');
/** \brief Constructor that initialize to all zeros
     \sa zeros
*/
  MXConstant(int n, int m); // all zeros

/** \brief  Clone function */
  virtual MXConstant* clone() const;

/** \brief  Print */
  virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
/* virtual void evaluateAdj();*/
  
/** \brief  Evaluate the second order derivative (adjoint of the forward gradient) and store the result in the dependency nodes */
/** \brief    virtual void evaluateFoA(); */

/** \brief  Access an element matrix style */
   double operator()(int i, int j) const;
   double& operator()(int i, int j);

/** \brief  Access an element vector style */
   double operator[](int k) const;
   double& operator[](int k);

/** \brief  Print */
  friend std::ostream& operator<<(std::ostream &stream, const MXConstant &x);

  virtual bool isConstant() const;
  
/** \brief  data member */
  protected:
  std::vector<double> data;

};

} // namespace CasADi


#endif // MX_CONSTANT_HPP
