#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
#! CasADi tutorial 1
#! ==================
#! This tutorial file explains the use of CasADi's Matrix<T> in a python context.
#! Matrix<T> is a general class for sparse matrices. We inspect it with the help of Matrix<double>
#! Let's start with the import statements to load CasADi.
from casadi import *
from numpy import *
#! Contructors & printing
#! --------------------------------------
#! The python name for Matrix<double> is DMatrix
a = DMatrix.zeros(3,4)
print a
#! The string representation shows only the structural non-zero entries. In this case there are none.
#! Let's make a DMatrix with some structural non-zero entries.
w = DMatrix(Sparsity(4,3,[0,2,2,3],[1,2,1]),[3,2.3,8])
print w
#! Internally, the Matrix<> class uses a Compressed Column Format which containts the offset to the first nonzero on each column ...
print "column offsets: ", w.colind()
#! ... the row for each nonzero ...
print "row: ", w.row()
#! ... and the nonzero data entries:
print "nonzeros: ", w.data()
#! Conversion
#! --------------
#! DMatrix can easily be converted into other data formats
print w.data()
print w.toArray()
print array(w)
print w.toMatrix()
print matrix(w)
print w.toCsc_matrix()

