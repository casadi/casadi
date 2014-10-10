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
from __future__ import division
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
import numpy as n
import unittest
from types import *
from helpers import *

scipy_available = True
try:
	from scipy.sparse import *
	from scipy import *
except:
	scipy_available = False
import warnings

class typemaptests(casadiTestCase):

  def setUp(self):
    pass
  
  def test_floordiv(self):
    self.message("make sure that floor_div raises errors")
    for x in [SXElement.sym("x"),MX.sym("x"),DMatrix([3]),SX.sym("x")]:
      for y in [2,2.0,DMatrix(3),numpy.matrix([2.0])]:
        print (x,y)
        self.assertRaises(Exception,lambda : x//y)
        self.assertRaises(Exception,lambda : y//x)
  
  def test_autoconversion(self):
    self.message("Auto conversion DMatrix")
    x=array([2.3])
    s = DMatrix([[1,2],[3,4]])
    n = array(s)
    
    self.checkarray(x[0]*s,s*x[0],"")
    self.checkarray(x[0]*s,n*x[0],"")
    
    self.checkarray(x[0]/s,1/(s/x[0]),"")
    self.checkarray(x[0]/s,x[0]/n,"")
    
    self.checkarray(x[0]-s,-(s-x[0]),"")
    self.checkarray(x[0]-s,x[0]-n,"")
    
    self.checkarray(x[0]+s,s+x[0],"")
    self.checkarray(x[0]+s,x[0]+n,"")
    
    w=array([2.3])[0]
    w+=s
    self.checkarray(w,2.3+n,"")

    w=array([2.3])[0]
    w-=s
    self.checkarray(w,2.3-n,"")
    
    w=array([2.3])[0]
    w*=s
    self.checkarray(w,2.3*n,"")
    
    w=array([2.3])[0]
    w/=s
    self.checkarray(w,2.3/n,"")
    
    x=[2.3]

    self.checkarray(x[0]*s,s*x[0],"")
    self.checkarray(x[0]*s,n*x[0],"")
    
    self.checkarray(x[0]/s,1/(s/x[0]),"")
    self.checkarray(x[0]/s,x[0]/n,"")
    
    self.checkarray(x[0]-s,-(s-x[0]),"")
    self.checkarray(x[0]-s,x[0]-n,"")
    
    self.checkarray(x[0]+s,s+x[0],"")
    self.checkarray(x[0]+s,x[0]+n,"")
    
    
    w=2.3
    w+=s
    self.checkarray(w,2.3+n,"")
    
    w=2.3
    w-=s
    self.checkarray(w,2.3-n,"")
    
    w=2.3
    w*=s
    self.checkarray(w,2.3*n,"")
    
    w=2.3
    w/=s
    self.checkarray(w,2.3/n,"")

  def test_autoconversionMX(self):
    self.message("Auto conversion MX")
    s = DMatrix([[1,2],[3,4]])
    x = SXElement(3)
    y = MX(3)
    
    def doit(z,s,fun):
      function = None
      
      if type(z) in [type(SXElement()),type(SX())]:
        ztype = [type(SXElement()),type(SX())]
        function = SXFunction
      
      if type(z) in [type(MX())]:
        ztype = [type(MX())]
        function = MXFunction
        
      r = fun(z,s)
            
      if type(z) is type(SXElement()) and type(s) is type(SXElement()):
        self.assertTrue(type(r) is type(SXElement()))
        

      self.assertTrue(type(r) in ztype)
      
      hasNum = True
      if type(s) in [type(SXElement()),type(MX()),type(SX())]:
        hasNum = False
      
      if hasNum:
        dummy = [1.3,2.7,9.4,1.0]

        f=function([z],[r])
        f.init()
        f.setInput(dummy[0:f.input().size()])
        f.evaluate()
        
        f_=function([z],[z])
        f_.init()
        f_.setInput(dummy[0:f.input().size()])
        f_.evaluate()
        

        self.checkarray(fun(f_.getOutput(),DMatrix(s)),f.getOutput(),"operation")
      else:
        dummy = [1.3,2.7,9.4,1.0]
        dummy2 = [0.3,2.4,1.4,1.7]
        
        f=function([z,s],[r])
        f.init()
        f.setInput(dummy[0:f.input(0).size()],0)
        f.setInput(dummy2[0:f.input(1).size()],1)
        f.evaluate()
        
        f_=function([z,s],[z,s])
        f_.init()
        f_.setInput(dummy[0:f.input(0).size()],0)
        f_.setInput(dummy2[0:f.input(1).size()],1)
        f_.evaluate()

        self.checkarray(fun(f_.getOutput(0),f_.getOutput(1)),f.getOutput(),"operation")
    
    
    def tests(z,s):
      doit(z,s,lambda z,s: -z)
      doit(z,s,lambda z,s: z+s)
      doit(z,s,lambda z,s: s+z)
      doit(z,s,lambda z,s: s*z)
      doit(z,s,lambda z,s: z*s)
      doit(z,s,lambda z,s: z-s)
      doit(z,s,lambda z,s: s-z)
      doit(z,s,lambda z,s: z/s)
      doit(z,s,lambda z,s: s/z)
      doit(z,s,lambda z,s: z**s)
      doit(z,s,lambda z,s: s**z)
      doit(z,s,lambda z,s: fmin(s,z))
      doit(z,s,lambda z,s: fmax(s,z))
      doit(z,s,lambda z,s: min(s,z))
      doit(z,s,lambda z,s: max(s,z))
      doit(z,s,lambda z,s: constpow(s,z))
      doit(z,s,lambda z,s: constpow(z,s))
      
    nums = [array([[1,2],[3,4]]),DMatrix([[1,2],[3,4]]), DMatrix(4), array(4),4.0,4]
        
    ## numeric & SX
    for s in nums:
      for z in [SXElement.sym("x"), SX.sym("x"), SX.sym("x",2,2)]:
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        tests(z,s)
       
    # numeric & MX
    for s in nums:
      for z in [MX.sym("x",2,2)]:
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        tests(z,s)
        
    # SXElement & SX
    for s in [SXElement.sym("x"), SX.sym("x"), SX.sym("x",2,2)]:
      for z in [SXElement.sym("x"),SX.sym("x"), SX.sym("x",2,2)]:
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        tests(z,s)
         
    ## MX & MX
    for s in [MX.sym("x"),MX.sym("x",2,2)]:
      for z in [MX.sym("x"),MX.sym("x",2,2)]:
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        tests(z,s)
        
    for (s,x,y) in [
                  (matrix([[1,2],[3,4]]),SX.sym("x",2,2),MX.sym("x",2,2))    
                  ]:
      for z,ztype in zip([x,y],[[type(SX()),type(SXElement())],[type(MX())]]):
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        doit(z,s,lambda z,s: -z)
        -s
        doit(z,s,lambda z,s: z+s)
        doit(z,s,lambda s,z: s+z)
        doit(z,s,lambda z,s: s*z)
        doit(z,s,lambda s,z: z*s)
        doit(z,s,lambda z,s: z-s)
        doit(z,s,lambda s,z: s-z)
        doit(z,s,lambda z,s: z/s)
        doit(z,s,lambda s,z: s/z)
        
if __name__ == '__main__':
    unittest.main()
