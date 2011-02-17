from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class ADtests(casadiTestCase):

  def setUp(self):
    x=SX("x")
    y=SX("y")
    z=SX("z")
    
    out=SXMatrix(5,1)
    out[0,0]=x
    out[2,0]=x+2*y**2
    out[4,0]=x+2*y**3+3*z**4
    

    inp=SXMatrix(5,1)
    inp[0,0]=x
    inp[2,0]=y
    inp[4,0]=z
    
    self.sxinputs = {
       "column" : {
            "dense": [[x,y,z]],
            "sparse": [inp] }
        , "row": {
            "dense":  [SXMatrix([x,y,z]).T],
            "sparse": [inp.T]
        }
    }
    
    xyz=MX("xyz",3,1)

    self.mxinputs = {
       "column" : {
            "dense": [xyz]
        }
    }
    
    self.mxoutputs = {
       "column": {
        "dense": [vertcat([xyz[0],xyz[0]+2*xyz[1]**2,xyz[0]+2*xyz[1]**3+3*xyz[2]**4])]
        }, "row": {
        "dense": [horzcat([xyz[0],xyz[0]+2*xyz[1]**2,xyz[0]+2*xyz[1]**3+3*xyz[2]**4])]
       }
    }
    
    
    self.sxoutputs = {
       "column": {
        "dense": [[x,x+2*y**2,x+2*y**3+3*z**4]],
        "sparse": [out]
        }, "row": {
          "dense":  [SXMatrix([x,x+2*y**2,x+2*y**3+3*z**4]).T],
          "sparse": [out.T]
       }
    }
    
    self.jacobians = {
      "dense" : {
        "dense" : lambda x,y,z: array([[1,0,0],[1,4*y,0],[1,6*y**2,12*z**3]]),
        "sparse" : lambda x,y,z: array([[1,0,0],[0,0,0],[1,4*y,0],[0,0,0],[1,6*y**2,12*z**3]])
        }
      ,
      "sparse" : {
        "dense" : lambda x,y,z: array([[1,0,0,0,0],[1,0,4*y,0,0],[1,0,6*y**2,0,12*z**3]]),
        "sparse" : lambda x,y,z:  array([[1,0,0,0,0],[0,0,0,0,0],[1,0,4*y,0,0],[0,0,0,0,0],[1,0,6*y**2,0,12*z**3]])
      }
    }
  
  
  def test_fwd(self):
    n=array([1.2,2.3,7])
    for inputshape in ["column","row"]:
      for outputshape in ["column","row"]:
        for inputtype in ["sparse","dense"]:
          for outputtype in ["sparse","dense"]:
            self.message("fwd AD on SX. Input %s %s, Output %s %s" % (inputtype,inputshape,outputtype,outputshape) )
            f=SXFunction(self.sxinputs[inputshape][inputtype],self.sxoutputs[outputshape][outputtype])
            f.init()
            f.input().set(n)
            self.assertEqual(f.fwdSeed().shape,f.input().shape,"fwdSeed shape")
            self.assertEqual(f.fwdSeed().size(),f.input().size(),"fwdSeed shape")
            J = self.jacobians[inputtype][outputtype](*n)
            for d in [array([1,0,0]),array([0,2,0]),array([1.2,4.8,7.9])]:
              f.fwdSeed().set(d)
              f.evaluate(1,0)
              seed = array(f.fwdSeed())
              if inputshape == "row":
                seed = seed.T
              sens = array(f.fwdSens())
              if outputshape == "row":
                sens = sens.T
              self.checkarray(sens,dot(J,seed),"AD")

  def test_adj(self):
    n=array([1.2,2.3,7])
    for inputshape in ["column","row"]:
      for outputshape in ["column","row"]:
        for inputtype in ["sparse","dense"]:
          for outputtype in ["sparse","dense"]:
            self.message("adj AD on SX. Input %s %s, Output %s %s" % (inputtype,inputshape,outputtype,outputshape) )
            f=SXFunction(self.sxinputs[inputshape][inputtype],self.sxoutputs[outputshape][outputtype])
            f.init()
            f.input().set(n)
            self.assertEqual(f.adjSeed().shape,f.output().shape,"adjSeed shape")
            self.assertEqual(f.adjSeed().size(),f.output().size(),"adjSeed shape")
            J = self.jacobians[inputtype][outputtype](*n)
            for d in [array([1,0,0]),array([0,2,0]),array([1.2,4.8,7.9])]:
              f.adjSeed().set(d)
              f.evaluate(0,1)
              seed = array(f.adjSeed())
              if outputshape == "row":
                seed = seed.T
              sens = array(f.adjSens())
              if inputshape == "row":
                sens = sens.T
              
              self.checkarray(sens,dot(J.T,seed),"AD")

  def test_fwdMX(self):
    n=array([1.2,2.3,7])
    for inputshape in ["column"]:
      for outputshape in ["column","row"]:
        for inputtype in ["dense"]:
          for outputtype in ["dense"]:
            self.message("fwd AD on MX. Input %s %s, Output %s %s" % (inputtype,inputshape,outputtype,outputshape) )
            f=MXFunction(self.mxinputs[inputshape][inputtype],self.mxoutputs[outputshape][outputtype])
            f.init()
            f.input().set(n)
            self.assertEqual(f.fwdSeed().shape,f.input().shape,"fwdSeed shape")
            self.assertEqual(f.fwdSeed().size(),f.input().size(),"fwdSeed shape")
            J = self.jacobians[inputtype][outputtype](*n)
            for d in [array([1,0,0]),array([0,2,0]),array([1.2,4.8,7.9])]:
              f.fwdSeed().set(d)
              f.evaluate(1,0)
              seed = array(f.fwdSeed())
              if inputshape == "row":
                seed = seed.T
              sens = array(f.fwdSens())
              if outputshape == "row":
                sens = sens.T
              self.checkarray(sens,dot(J,seed),"AD")    

  def test_adj(self):
    n=array([1.2,2.3,7])
    for inputshape in ["column"]:
      for outputshape in ["column","row"]:
        for inputtype in ["dense"]:
          for outputtype in ["dense"]:
            self.message("adj AD on MX. Input %s %s, Output %s %s" % (inputtype,inputshape,outputtype,outputshape) )
            f=MXFunction(self.mxinputs[inputshape][inputtype],self.mxoutputs[outputshape][outputtype])
            f.init()
            f.input().set(n)
            self.assertEqual(f.adjSeed().shape,f.output().shape,"adjSeed shape")
            self.assertEqual(f.adjSeed().size(),f.output().size(),"adjSeed shape")
            J = self.jacobians[inputtype][outputtype](*n)
            for d in [array([1,0,0]),array([0,2,0]),array([1.2,4.8,7.9])]:
              f.adjSeed().set(d)
              f.evaluate(0,1)
              seed = array(f.adjSeed())
              if outputshape == "row":
                seed = seed.T
              sens = array(f.adjSens())
              if inputshape == "row":
                sens = sens.T
              
              self.checkarray(sens,dot(J.T,seed),"AD")
if __name__ == '__main__':
    unittest.main()

