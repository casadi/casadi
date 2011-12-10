from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *
from casadi.tools import *

class Toolstests(casadiTestCase):

  def test_variables(self):
      self.message("Variables")
      # -*- coding: utf-8 -*-
             
      p = Variables()

      p.x = ssym("x",2)
      p.z = ssym("z",2,4)
      p.y = ssym("y",3,2)

      p.xother = Variables()
      p.xother.a = SX("x")
      p.xother.b = diag(ssym("[a,b]"))

      self.assertEqual(p.o_x,0)
      self.assertEqual(p.o_xother,2)
      self.assertEqual(p.xother.o_a,0)
      self.assertEqual(p.xother.o_b,1)
      self.assertEqual(p.o_y,5)
      self.assertEqual(p.o_z,11)
      
      self.assertEqual(p.I_x,0,"Index")
      self.assertEqual(p.I_y,2,"Index")
      self.assertEqual(p.I_z,3,"Index")
      self.assertEqual(p.xother.I_a,0,"Index")
      self.assertEqual(p.xother.I_b,1,"Index")
      
      self.checkarray(array(p.i_x),DMatrix([[0],[1]]),"index")
      self.checkarray(array(p.i_y),DMatrix([[5,8],[6,9],[7,10]]),"index")
      self.checkarray(array(p.i_z),DMatrix([[11,13,15,17],[12,14,16,18]]),"index")
      self.checkarray(array(p.xother.i_a),DMatrix(0),"index")
      self.checkarray(array(p.xother.i_b),DMatrix([[1,0],[0,2]]),"index")

      p.x_ = [5,8]

      p.xother.b_.setAll(7)

      p.z_ = DMatrix([[1,2,3,4],[5,6,7,8]])


      A = p.vecNZcat_()
      
      self.checkarray(A,DMatrix([5,8,0,7,7,0,0,0,0,0,0,1,5,2,6,3,7,4,8]),"vecNZcat")
      
      A = p.veccat_()
      
      self.checkarray(A,DMatrix([5,8,0,7,0,0,7,0,0,0,0,0,0,1,5,2,6,3,7,4,8]),"veccat")
      
      self.checkarray(A[p.i_z],p.z_,"indexing round trip")
      self.checkarray(A[p.o_xother + p.xother.i_b],p.xother.b_,"indexing round trip 2")

      p = Variables()
      p.a = ssym("a",2)
      p.b = []
      p.b.append(ssym("b1",3))
      p.b.append(ssym("b2",3))
      p.c = ssym("c")
      
      self.checkarray(array(p.i_a),DMatrix([[0],[1]]),"index")
      self.checkarray(array(p.i_b[0]),DMatrix([[2],[3],[4]]),"index")
      self.checkarray(array(p.i_b[1]),DMatrix([[5],[6],[7]]),"index")
      self.checkarray(array(p.c),DMatrix(8),"index")

      self.assertEqual(p.o_a,0,"Offset")
      self.assertEqual(p.o_b[0],2,"Offset")
      self.assertEqual(p.o_b[1],5,"Offset")
      self.assertEqual(p.o_c,8,"Offset")
      
      self.assertEqual(p.I_a,0,"Index")
      self.assertEqual(p.I_b[0],1,"Index")
      self.assertEqual(p.I_b[1],2,"Index")
      self.assertEqual(p.I_c,3,"Index")
      
      
      p.b_[1].setAll(4)
      
      A = p.vecNZcat_()
 
      self.checkarray(A,DMatrix([0,0,0,0,0,4,4,4,0]),"vecNZcat")
     
      p.b_[0].setAll(3)
      A = p.veccat_()

      self.checkarray(A,DMatrix([0,0,3,3,3,4,4,4,0]),"vecNZcat")
      
      
      
      p = Variables()
      p.a = msym("a",2)
      #p.b = []
      #p.b.append(msym("b1",3))
      #p.b.append(msym("b2",3))
      #p.c = msym("c")
      
       

if __name__ == '__main__':
    unittest.main()

