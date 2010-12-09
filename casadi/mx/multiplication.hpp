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

#ifndef MULTIPLICATION_HPP
#define MULTIPLICATION_HPP

#include "mx_node.hpp"

namespace CasADi{
/** \brief Represents matrix products
  \author Joel Andersson 
  \date 2010
*/
class Multiplication : public MXNode{
public:

/** \brief  Constructor */
Multiplication(const MX& x, const MX& y);

/** \brief  Clone function */
virtual Multiplication* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
//virtual void evaluateAdj();
  
protected:

/** \brief  Form the product of the ni-by-nk matrix t1 and the nk-by-nj matrix t2 and add the result to the ni-by-nj matrix t3 */
  void matrix_matrix_mult(
      const std::vector<double>& t1, 
      const std::vector<double>& t2, 
      std::vector<double>& t3);

/** \brief  Form the product of the ni-by-nj matrix t3 and the transpose of the nk-by-nj matrix t2 and add the result to the ni-by-nk matrix t1 */
  void matrix_matrix_mult1(
      std::vector<double>& t1, 
      const std::vector<double>& t2, 
      const std::vector<double>& t3);

/** \brief  Form the product of the transpose of the ni-by-nk matrix t1 and the ni-by-nj matrix t3 and add the result to the nk-by-nj matrix t2 */
  void matrix_matrix_mult2(
      const std::vector<double>& t1, 
      std::vector<double>& t2, 
      const std::vector<double>& t3);

/** \brief  Dimensions */
  int ni, nk, nj;

};

} // namespace CasADi


#endif // MULTIPLICATION_HPP
