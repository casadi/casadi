#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class SDPtests(casadiTestCase):
  def test_example1(self):
    self.message("Example1")
    # Originates from http://sdpa.indsys.chuo-u.ac.jp/sdpa/files/sdpa-c.6.2.0.manual.pdf
    b = DMatrix([48,-8,20])

    A = vertcat([DMatrix([[10,4],[4,0]]),DMatrix([[0,0],[0,-8]]),DMatrix([[0,-8],[-8,-2]])])

    makeSparse(A)

    A.printMatrix()

    C = DMatrix([[-11,0],[0,23]])

    makeSparse(C)

    dsp = DSDPSolver(C.sparsity(),A.sparsity())

    dsp.init()

    dsp.input(SDP_C).set(C)
    dsp.input(SDP_B).set(b)
    dsp.input(SDP_A).set(A)

    dsp.evaluate()
    
    self.checkarray(dsp.output(SDP_PRIMAL_COST),DMatrix(-41.9),digits=5)
    self.checkarray(dsp.output(SDP_DUAL_COST),DMatrix(-41.9),digits=5)
    self.checkarray(dsp.output(SDP_PRIMAL),DMatrix([-1.1,-2.7375,-0.55]),digits=5)
    
    self.checkarray(dsp.output(SDP_DUAL),DMatrix([[5.9,-1.375],[-1.375,1]]),digits=5)
    self.checkarray(dsp.output(SDP_PRIMAL_P),DMatrix.zeros(2,2),digits=5)
    
  def test_example2(self):
    self.message("Example2")
    # Originates from http://sdpa.indsys.chuo-u.ac.jp/sdpa/files/sdpa-c.6.2.0.manual.pdf
    b = DMatrix([1.1, -10, 6.6 , 19 , 4.1])


    C = blockdiag([DMatrix([[-1.4,-3.2],[-3.2,-28]]),DMatrix([[15,-12,2.1],[-12,16,-3.8],[2.1,-3.8,15]]),1.8,-4.0]);
    
    sp = C.sparsity()
    
    flatdata = [[0.5,5.2,5.2,-5.3,7.8,-2.4,6.0,-2.4,4.2,6.5,6.0,6.5,2.1,-4.5,-3.5],
    [1.7,7.0,7.0,-9.3,-1.9,-0.9,-1.3,-0.9,-0.8,-2.1,-1.3,-2.1,4.0,-0.2,-3.7],
 [6.3,-7.5,-7.5,-3.3,0.2,8.8,5.4,8.8,3.4,-0.4,5.4,-0.4,7.5,-3.3,-4.0],
  [-2.4,-2.5,-2.5,-2.9,3.4,-3.2,-4.5,-3.2,3.0,-4.8,-4.5,-4.8,3.6,4.8,9.7],
  [-6.5,-5.4,-5.4,-6.6,6.7,-7.2,-3.6,-7.2,7.3,-3.0,-3.6,-3.0,-1.4,6.1,-1.5]]

    A = vertcat([DMatrix(sp,data) for data in flatdata])
    makeSparse(A)


    dsp = DSDPSolver(C.sparsity(),A.sparsity())

    dsp.init()

    dsp.input(SDP_C).set(C)
    dsp.input(SDP_B).set(b)
    dsp.input(SDP_A).set(A)

    dsp.evaluate()
    DMatrix.setPrecision(10)
    self.checkarray(dsp.output(SDP_PRIMAL_COST),DMatrix(3.20626934048e1),digits=5)
    self.checkarray(dsp.output(SDP_DUAL_COST),DMatrix(3.20626923535e1),digits=5)
    self.checkarray(dsp.output(SDP_PRIMAL),DMatrix([1.551644595,0.6709672545,0.9814916693,1.406569511,0.9421687787]),digits=5)
    
    self.checkarray(dsp.output(SDP_DUAL),DMatrix(sp,[2.640261206,0.5605636589,0.5605636589,3.717637107,0.7615505416,-1.513524657,1.139370202,-1.513524657,3.008016978,-2.264413045,1.139370202,-2.264413045,1.704633559,0,0]),digits=5)
    self.checkarray(dsp.output(SDP_PRIMAL_P),DMatrix(sp,[0,0,0,0,7.119155551,5.024671489,1.916294752,5.024671489,4.414745792,2.506021978,1.916294752,2.506021978,2.048124139,0.3432465654,4.391169489]),digits=5)
    
if __name__ == '__main__':
    unittest.main()

