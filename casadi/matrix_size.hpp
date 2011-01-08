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

/*
 NOTE: THIS CLASS HAS BEEN MARKED FOR DELETION! Joel
*/
  
  
/** \brief  Container class for size of a matrix
  \author Joel Andersson 
  \date 2010	
*/
class MatrixSize{
  public:
  MatrixSize() : nrow(0), ncol(0){ }  // default constructor 0-by-0 matrix
  MatrixSize(int nrow_, int ncol_) : nrow(nrow_), ncol(ncol_){ }  // size n-by-m
 
#ifndef SWIG
 friend std::ostream& operator<<(std::ostream &stream, const MatrixSize &x) {
      return stream << "[" << x.nrow << "," << x.ncol << "]"; 
  }
#endif  
  
  int nrow, ncol;
};

#ifdef SWIG
%typemap(in) const MatrixSize & {
       $1 = new CasADi::MatrixSize(PyTupleToMatrixSize($input));
}

%typemap(freearg) const MatrixSize & {
    if ($1) {
        delete $1;
    }
}
#endif // SWIG


} // namespace CasADi

#ifdef SWIG
%inline %{
CasADi::MatrixSize PyTupleToMatrixSize(PyObject* tup) {
                  if (!PyTuple_Check(tup))  throw CasADi::CasadiException("__getitem__: expecting tuple");
      if(PyTuple_Size(tup)!=2) throw CasADi::CasadiException("__getitem__: not 2D");
      return CasADi::MatrixSize(PyInt_AsLong(PyTuple_GetItem(tup,0)),PyInt_AsLong(PyTuple_GetItem(tup,1)));
}
%}
#endif // SWIG



#endif // MATRIX_SIZE_HPP
