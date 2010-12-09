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

#ifndef MATRIX_SIZE_HPP
#define MATRIX_SIZE_HPP

#include <iostream>

namespace CasADi{

/** \brief  Container class for size of a matrix
  \author Joel Andersson 
  \date 2010	
*/
class MatrixSize{
  public:
  MatrixSize() : nrow(0), ncol(0){ }  // default constructor 0-by-0 matrix
  MatrixSize(int nrow_, int ncol_) : nrow(nrow_), ncol(ncol_){ }  // size n-by-m

/*  bool operator==(const MatrixSize &x) const{
     return nrow == x.nrow && ncol == x.ncol;
  }

  bool operator!=(const MatrixSize &x) const{
      return nrow != x.nrow || ncol != x.ncol;
  }*/
 
 friend std::ostream& operator<<(std::ostream &stream, const MatrixSize &x) {
      return stream << "[" << x.nrow << "," << x.ncol << "]"; 
  }

  int nrow, ncol;
};

} // namespace CasADi

#endif // MATRIX_SIZE_HPP
