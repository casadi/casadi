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

#ifndef MATRIX_ELEMENT_HPP
#define MATRIX_ELEMENT_HPP

#include "reordering.hpp"

namespace CasADi{
/** \brief  A submatrix (currently element)
    NOTE This class should inherit from reordering!
    \author Joel Andersson 
    \date 2010-2011
*/
class MatrixElement : public Reordering {
  public:

    /** \brief  Constructor */
    MatrixElement(const MX& x, int i, int j);

    /** \brief  Clone function */
    virtual MatrixElement* clone() const;

    /** \brief  Print */
    virtual void print(std::ostream &stream=std::cout) const;

  protected:
    int i_, j_;
};

} // namespace CasADi


#endif // MATRIX_ELEMENT_HPP
