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
from casadi import *

x = MX.sym("x",5)
y = MX.sym("y")
x2 = SX.sym("x",5)
y2 = SX.sym("y")

fcn = MXFunction([x,y],[4*vertcat((x[2:5],x[0:2])) + y*x])
fcn.init()
js = IMatrix(fcn.jacSparsity(),1)
js.printDense()

fcn2 = SXFunction([x2,y2],[4*vertcat((x2[2:5],x2[0:2])) + y2*x2])
fcn2.init()
js2 = IMatrix(fcn2.jacSparsity(),1)
js2.printDense()

fcn3 = MXFunction([x,y],fcn2.call([x,y]))
fcn3.init()
js3 = IMatrix(fcn3.jacSparsity(),1)
js3.printDense()
