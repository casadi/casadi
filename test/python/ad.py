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
import casadi as c
import numpy
from numpy import random, array
import unittest
from types import *
from helpers import *
import itertools

warnings.filterwarnings("ignore",category=DeprecationWarning)

class ADtests(casadiTestCase):

  def setUp(self):
    x=SX.sym("x")
    y=SX.sym("y")
    z=SX.sym("z")
    w=SX.sym("w")

    out=SX(6,1)
    out[0,0]=x
    out[2,0]=x+2*y**2
    out[4,0]=x+2*y**3+3*z**4
    out[5,0]=w

    inp=SX(6,1)
    inp[0,0]=x
    inp[2,0]=y
    inp[4,0]=z
    inp[5,0]=w

    sp = Sparsity(1,6,[0, 1, 1, 2, 2, 3, 4],[0, 0, 0, 0]).T
    spT = Sparsity(6,1,[0, 4],[0, 2, 4, 5]).T

    self.sxinputs = {
       "column" : {
            "dense": [vertcat(*[x,y,z,w])],
            "sparse": [inp] }
        , "row": {
            "dense":  [vertcat(*[x,y,z,w]).T],
            "sparse": [inp.T]
       }, "matrix": {
          "dense": [c.reshape(vertcat(*[x,y,z,w]),2,2)],
          "sparse": [c.reshape(inp,3,2)]
        }
    }

    self.mxinputs = {
       "column" : {
            "dense": [MX.sym("xyzw",4,1)],
            "sparse": [MX.sym("xyzw",sp)]
        },
        "row" : {
            "dense": [MX.sym("xyzw",1,4)],
            "sparse": [MX.sym("xyzw",spT)]
        },
        "matrix": {
            "dense": [MX.sym("xyzw",2,2)],
            "sparse": [MX.sym("xyzw",c.reshape(inp,3,2).sparsity())]
        }
    }

    def temp1(xyz):
      X=MX(6,1)
      X[0,0]=xyz.nz[0]
      X[2,0]=xyz.nz[0]+2*xyz.nz[1]**2
      X[4,0]=xyz.nz[0]+2*xyz.nz[1]**3+3*xyz.nz[2]**4
      X[5,0]=xyz.nz[3]
      return [X]

    def temp2(xyz):
      X=MX(1,6)
      X[0,0]=xyz.nz[0]
      X[0,2]=xyz.nz[0]+2*xyz.nz[1]**2
      X[0,4]=xyz.nz[0]+2*xyz.nz[1]**3+3*xyz.nz[2]**4
      X[0,5]=xyz.nz[3]
      return [X]

    def testje(xyz):
      print(vertcat(*[xyz.nz[0],xyz.nz[0]+2*xyz.nz[1]**2,xyz.nz[0]+2*xyz.nz[1]**3+3*xyz.nz[2]**4,xyz.nz[3]]).shape)

    self.mxoutputs = {
       "column": {
        "dense":  lambda xyz: [vertcat(*[xyz.nz[0],xyz.nz[0]+2*xyz.nz[1]**2,xyz.nz[0]+2*xyz.nz[1]**3+3*xyz.nz[2]**4,xyz.nz[3]])],
        "sparse": temp1
        }, "row": {
        "dense": lambda xyz: [horzcat(*[xyz.nz[0],xyz.nz[0]+2*xyz.nz[1]**2,xyz.nz[0]+2*xyz.nz[1]**3+3*xyz.nz[2]**4,xyz.nz[3]])],
        "sparse": temp2
       },
       "matrix": {
          "dense": lambda xyz: [c.reshape(vertcat(*[xyz.nz[0],xyz.nz[0]+2*xyz.nz[1]**2,xyz.nz[0]+2*xyz.nz[1]**3+3*xyz.nz[2]**4,xyz.nz[3]]),(2,2))],
          "sparse": lambda xyz: [c.reshape(temp1(xyz)[0],(3,2))]
       }
    }


    self.sxoutputs = {
       "column": {
        "dense": [vertcat(*[x,x+2*y**2,x+2*y**3+3*z**4,w])],
        "sparse": [out]
        }, "row": {
          "dense":  [vertcat(*[x,x+2*y**2,x+2*y**3+3*z**4,w]).T],
          "sparse": [out.T]
      }, "matrix" : {
          "dense":  [c.reshape(vertcat(*[x,x+2*y**2,x+2*y**3+3*z**4,w]),2,2)],
          "sparse": [c.reshape(out,3,2)]
      }
    }

    self.jacobians = {
      "dense" : {
        "dense" : lambda x,y,z,w: array([[1,0,0,0],[1,4*y,0,0],[1,6*y**2,12*z**3,0],[0,0,0,1]]),
        "sparse" : lambda x,y,z,w: array([[1,0,0,0],[0,0,0,0],[1,4*y,0,0],[0,0,0,0],[1,6*y**2,12*z**3,0],[0,0,0,1]])
        }
      ,
      "sparse" : {
        "dense" : lambda x,y,z,w: array([[1,0,0,0,0,0],[1,0,4*y,0,0,0],[1,0,6*y**2,0,12*z**3,0],[0,0,0,0,0,1]]),
        "sparse" : lambda x,y,z,w:  array([[1,0,0,0,0,0],[0,0,0,0,0,0],[1,0,4*y,0,0,0],[0,0,0,0,0,0],[1,0,6*y**2,0,12*z**3,0],[0,0,0,0,0,1]])
      }
    }

  def test_SXeval_sx(self):
    n=array([1.2,2.3,7,1.4])
    for inputshape in ["column","row","matrix"]:
      for outputshape in ["column","row","matrix"]:
        for inputtype in ["dense","sparse"]:
          for outputtype in ["dense","sparse"]:
            self.message("eval_sx on SX. Input %s %s, Output %s %s" % (inputtype,inputshape,outputtype,outputshape) )
            f=Function("f", self.sxinputs[inputshape][inputtype],self.sxoutputs[outputshape][outputtype])
            f_in = DM(f.sparsity_in(0),n)
            r = f(f_in)
            J = self.jacobians[inputtype][outputtype](*n)

            seeds = [[1,0,0,0],[0,2,0,0],[1.2,4.8,7.9,4.6]]

            y = SX.sym("y",f.sparsity_in(0))

            fseeds = [DM(f.sparsity_in(0),x) for x in seeds]
            aseeds = [DM(f.sparsity_out(0),x) for x in seeds]
            res = f(y)
            fwdsens = forward([res], [y], [[x] for x in fseeds])
            adjsens = reverse([res], [y], [[x] for x in aseeds])
            fwdsens = [x[0] for x in fwdsens]
            adjsens = [x[0] for x in adjsens]

            fe = Function("fe", [y], [res])

            re = fe(f_in)

            self.checkarray(r,re)

            for sens,seed in zip(fwdsens,fseeds):
              fe = Function("fe", [y],[sens])
              re = fe(f_in)
              self.checkarray(c.vec(re),mtimes(J,c.vec(seed)),"AD")

            for sens,seed in zip(adjsens,aseeds):
              fe = Function("fe", [y],[sens])
              re = fe(f_in)
              self.checkarray(c.vec(re),mtimes(J.T,c.vec(seed)),"AD")

  def test_MXeval_mx(self):
    n=array([1.2,2.3,7,1.4])
    for inputshape in ["column","row","matrix"]:
      for outputshape in ["column","row","matrix"]:
        for inputtype in ["dense","sparse"]:
          for outputtype in ["dense","sparse"]:
            self.message("eval_mx on MX. Input %s %s, Output %s %s" % (inputtype,inputshape,outputtype,outputshape) )
            f=Function("f", self.mxinputs[inputshape][inputtype],self.mxoutputs[outputshape][outputtype](self.mxinputs[inputshape][inputtype][0]))
            f_in = DM(f.sparsity_in(0),n)
            r = f(f_in)
            J = self.jacobians[inputtype][outputtype](*n)

            seeds = [[1,0,0,0],[0,2,0,0],[1.2,4.8,7.9,4.6]]

            y = MX.sym("y",f.sparsity_in(0))

            fseeds = [DM(f.sparsity_in(0),x) for x in seeds]
            aseeds = [DM(f.sparsity_out(0),x) for x in seeds]
            res = f(y)
            fwdsens = forward([res],[y], [[x] for x in fseeds])
            adjsens = reverse([res],[y], [[x] for x in aseeds])
            fwdsens = [x[0] for x in fwdsens]
            adjsens = [x[0] for x in adjsens]

            fe = Function('fe', [y], [res])

            re = fe(f_in)

            self.checkarray(r,re)

            for sens,seed in zip(fwdsens,fseeds):
              fe = Function("fe", [y],[sens])
              re = fe(f_in)
              self.checkarray(c.vec(re),mtimes(J,c.vec(seed)),"AD")

            for sens,seed in zip(adjsens,aseeds):
              fe = Function("fe", [y],[sens])
              re = fe(f_in)
              self.checkarray(c.vec(re),mtimes(J.T,c.vec(seed)),"AD")

  @known_bug()  # Not implemented
  def test_MXeval_sx(self):
    n=array([1.2,2.3,7,1.4])
    for inputshape in ["column","row","matrix"]:
      for outputshape in ["column","row","matrix"]:
        for inputtype in ["dense","sparse"]:
          for outputtype in ["dense","sparse"]:
            self.message("eval_sx on MX. Input %s %s, Output %s %s" % (inputtype,inputshape,outputtype,outputshape) )
            f=Function("f", self.mxinputs[inputshape][inputtype],self.mxoutputs[outputshape][outputtype](self.mxinputs[inputshape][inputtype][0]))
            f_in = [0]*f.n_in();f_in[0]=n
            f_out = f(f_in)
            r = f_out[0]
            J = self.jacobians[inputtype][outputtype](*n)

            seeds = [[1,0,0,0],[0,2,0,0],[1.2,4.8,7.9,4.6]]

            y = SX.sym("y",f.sparsity_in(0))

            fseeds = [DM(f.sparsity_in(0),x) for x in seeds]
            aseeds = [DM(f.sparsity_out(0),x) for x in seeds]
            res = f(y)
            fwdsens = forward([res],[y],[[x] for x in fseeds])
            adjsens = reverse([res],[y],[[x] for x in aseeds])
            fwdsens = [x[0] for x in fwdsens]
            adjsens = [x[0] for x in adjsens]

            fe = Function("fe", [y], [res])

            fe_in = [0]*fe.n_in();fe_in[0]=n
            fe_out = fe(fe_in)

            self.checkarray(r,fe_out[0])

            for sens,seed in zip(fwdsens,fseeds):
              fe = Function("fe", [y],[sens])
              fe_in = [0]*fe.n_in();fe_in[0]=n
              fe_out = fe(fe_in)
              self.checkarray(c.vec(fe_out[0].T),mtimes(J,c.vec(seed.T)),"AD")

            for sens,seed in zip(adjsens,aseeds):
              fe = Function("fe", [y],[sens])
              fe_in = [0]*fe.n_in();fe_in[0]=n
              fe_out = fe(fe_in)
              self.checkarray(c.vec(fe_out[0].T),mtimes(J.T,c.vec(seed.T)),"AD")

  def test_MXeval_sx_reduced(self):
    n=array([1.2,2.3,7,1.4])
    for inputshape in ["column","row","matrix"]:
      for outputshape in ["column","row","matrix"]:
        for inputtype in ["dense","sparse"]:
          for outputtype in ["dense","sparse"]:
            self.message("eval_sx on MX. Input %s %s, Output %s %s" % (inputtype,inputshape,outputtype,outputshape) )
            f=Function("f", self.mxinputs[inputshape][inputtype],self.mxoutputs[outputshape][outputtype](self.mxinputs[inputshape][inputtype][0]))
            f_in = DM(f.sparsity_in(0),n)
            r = f(f_in)

            y = SX.sym("y",f.sparsity_in(0))

            res = f(y)
            fwdsens = forward([res],[y],[])
            adjsens = reverse([res],[y],[])

            fe = Function("fe", [y],[res])

            re = f(f_in)

            self.checkarray(r,re)

  def test_Jacobian(self):
    n=array([1.2,2.3,7,4.6])
    for inputshape in ["column","row","matrix"]:
      for outputshape in ["column","row","matrix"]:
        for inputtype in ["dense","sparse"]:
          for outputtype in ["dense","sparse"]:
            for mode in ["forward","reverse"]:
              self.message(" %s Jacobian on SX. Input %s %s, Output %s %s" % (mode,inputtype,inputshape,outputtype,outputshape) )
              opts = {}
              opts["ad_weight"] = 0 if mode=='forward' else 1
              opts["ad_weight_sp"] = 0 if mode=='forward' else 1
              f=Function("f", self.sxinputs[inputshape][inputtype],self.sxoutputs[outputshape][outputtype], opts)
              Jf=f.jacobian_old(0,0)
              J_in = DM(f.sparsity_in(0),n)
              Jout,_ = Jf(J_in)
              J = self.jacobians[inputtype][outputtype](*n)
              self.checkarray(array(Jout),J,"Jacobian\n Mode: %s\n Input: %s %s\n Output: %s %s"% (mode, inputshape, inputtype, outputshape, outputtype))

  def test_jacobianSX(self):
    n=array([1.2,2.3,7,4.6])
    for inputshape in ["column","row","matrix"]:
      for outputshape in ["column","row","matrix"]:
        for inputtype in ["dense","sparse"]:
          for outputtype in ["dense","sparse"]:
            self.message("jacobian on SX (SCT). Input %s %s, Output %s %s" % (inputtype,inputshape,outputtype,outputshape) )
            Jf=Function("Jf",
              self.sxinputs[inputshape][inputtype],
              [
                  jacobian(
                    SX(self.sxoutputs[outputshape][outputtype][0]),
                    SX(self.sxinputs[inputshape][inputtype][0])
                  )
              ]
            )
            J_in = DM(Jf.sparsity_in(0),n)
            J_out = Jf(J_in)
            J = self.jacobians[inputtype][outputtype](*n)
            self.checkarray(array(J_out),J,"jacobian")

  def test_jacsparsity(self):
    n=array([1.2,2.3,7,4.6])
    for inputshape in ["column","row","matrix"]:
      for outputshape in ["column","row","matrix"]:
        for inputtype in ["dense","sparse"]:
          for outputtype in ["dense","sparse"]:
            self.message("jacsparsity on SX. Input %s %s, Output %s %s" % (inputtype,inputshape,outputtype,outputshape) )
            f=Function("f", self.sxinputs[inputshape][inputtype],self.sxoutputs[outputshape][outputtype])
            J = self.jacobians[inputtype][outputtype](*n)
            self.checkarray(DM.ones(f.sparsity_jac(0, 0)),array(J!=0,int),"jacsparsity")

  def test_JacobianMX(self):
    n=array([1.2,2.3,7,4.6])
    for inputshape in ["column","row","matrix"]:
      for outputshape in ["column","row","matrix"]:
        for inputtype in ["dense","sparse"]:
          for outputtype in ["dense","sparse"]:
            for mode in ["forward","reverse"]:
              self.message("adj AD on MX. Input %s %s, Output %s %s" % (inputtype,inputshape,outputtype,outputshape) )
              opts = {}
              opts["ad_weight"] = 0 if mode=='forward' else 1
              opts["ad_weight_sp"] = 0 if mode=='forward' else 1
              f=Function("f", self.mxinputs[inputshape][inputtype],self.mxoutputs[outputshape][outputtype](self.mxinputs[inputshape][inputtype][0]), opts)
              Jf=f.jacobian_old(0,0)
              J_in = DM(f.sparsity_in(0),n)
              J_out,_ = Jf(J_in)
              J = self.jacobians[inputtype][outputtype](*n)
              self.checkarray(J_out,J,"Jacobian\n Mode: %s\n Input: %s %s\n Output: %s %s"% (mode, inputshape, inputtype, outputshape, outputtype))

  def test_jacsparsityMX(self):
    n=array([1.2,2.3,7,4.6])
    for inputshape in ["column","row","matrix"]:
      for outputshape in ["column","row","matrix"]:
        for inputtype in ["dense","sparse"]:
          for outputtype in ["dense","sparse"]:
            for mode in ["forward","reverse"]:
              self.message(" %s jacobian on MX (SCT). Input %s %s, Output %s %s" % (mode,inputtype,inputshape,outputtype,outputshape) )
              opts = {}
              opts["ad_weight"] = 0 if mode=='forward' else 1
              opts["ad_weight_sp"] = 0 if mode=='forward' else 1
              f=Function("f", self.mxinputs[inputshape][inputtype],self.mxoutputs[outputshape][outputtype](self.mxinputs[inputshape][inputtype][0]), opts)
              Jf=f.jacobian_old(0,0)
              J_in = DM(f.sparsity_in(0),n)
              J_out,_ = Jf(J_in)
              J = self.jacobians[inputtype][outputtype](*n)
              self.checkarray(array(J_out),J,"jacobian")
              self.checkarray(array(DM.ones(f.sparsity_jac(0, 0))),array(J!=0,int),"jacsparsity")



  def test_hessian(self):
    self.message("Jacobian chaining")
    x=SX.sym("x")
    y=SX.sym("y")
    z=SX.sym("z")
    n=array([1.2,2.3,7])
    f=Function("f", [vertcat(*[x,y,z])],[vertcat(*[x+2*y**3+3*z**4])])
    J=f.jacobian_old(0,0)
    m=MX.sym("m",3,1)
    JT,_ = J(m)
    JT = Function("JT", [m],[JT.T])
    JT(n)
    H = JT.jacobian_old(0,0)
    H(n)
    #H_out = H(H_in)

    #print array(JT_out[0])
    #print array(H_out[0])

  def test_bugshape(self):
    self.message("shape bug")
    x=SX.sym("x")
    y=SX.sym("y")

    inp=SX(5,1)
    inp[0,0]=x
    inp[3,0]=y

    f=Function("f", [inp],[vertcat(*[x+y,x,y])])
    J=f.jacobian_old(0,0)
    J(DM(f.sparsity_in(0),[2,7]))

    self.assertEqual(f.size1_out(0),3,"Jacobian shape bug")
    self.assertEqual(f.size2_out(0),1,"Jacobian shape bug")


  def test_bugglibc(self):
    self.message("Code that used to throw a glibc error")
    x=SX.sym("x")
    y=SX.sym("y")

    inp=SX(5,1)
    inp[0,0]=x
    inp[3,0]=y

    f=Function("f", [inp],[vertcat(*[x+y,x,y])])
    J=f.jacobian_old(0,0)
    J_in = DM(f.sparsity_in(0),[2,7])
    J_out,_ = J(J_in)

    f=Function("f", [inp],[vertcat(*[x+y,x,y])])
    J=f.jacobian_old(0,0)

  @memory_heavy()
  def test_MX(self):

    x = MX.sym("x",2)
    y = MX.sym("y",2,2)

    f1 = Function("f1", [x,y],[x+y[0,0],mtimes(y,x)])

    f2 = Function("f2", [x,y],[mtimes(MX.zeros(0,2),x)])

    f3 = Function("f3", [x,y],[MX.zeros(0,0),mtimes(y,x)])

    f4 = Function("f4", [x,y],[MX.zeros(0,2),mtimes(y,x)])

    ndir = 2

    in1 = [x,y]
    v1 = [DM([1.1,1.3]),DM([[0.7,1.5],[2.1,0.9]])]

    w=x[:]
    w[1]*=2

    w2=x[:]
    w2[1]*=x[0]

    ww=x[:]
    ww[[0,1]]*=x

    wwf=x[:]
    wwf[[1,0]]*=x

    wwr=x[:]
    wwr[[0,0,1,1]]*=2

    yy=y[:,:]

    yy[:,0] = x

    yy2=y[:,:]

    yy2[:,0] = x**2

    yyy=y[:,:]

    yyy[[1,0],0] = x

    yyy2=y[:,:]

    yyy2[[1,0],0] = x**2


    def remove_first(x):
      ret = DM(x)
      if ret.numel()>0:
        ret[0,0] = DM(1,1)
        return ret.sparsity()
      else:
        return ret.sparsity()

    def remove_last(x):
      ret = DM(x)
      if ret.nnz()>0:
        ret[ret.sparsity().row()[-1],ret.sparsity().get_col()[-1]] = DM(1,1)
        return ret.sparsity()
      else:
        return x

    spmods = [lambda x: x , remove_first, remove_last]

    # TODO: sparse seeding

    for inputs,values,out, jac in [
          (in1,v1,x,DM.eye(2)),
          (in1,v1,x.T,DM.eye(2)),
          (in1,v1,x**2,2*c.diag(x)),
          (in1,v1,(x**2).attachAssert(True),2*c.diag(x)),
          (in1,v1,(x**2).T,2*c.diag(x)),
          (in1,v1,c.reshape(x,(1,2)),DM.eye(2)),
          (in1,v1,c.reshape(x**2,(1,2)),2*c.diag(x)),
          (in1,v1,x+y.nz[0],DM.eye(2)),
          (in1,v1,x+y[0,0],DM.eye(2)),
          (in1,v1,x+x,2*DM.eye(2)),
          (in1,v1,x**2+x,2*c.diag(x)+DM.eye(2)),
          (in1,v1,x*x,2*c.diag(x)),
          (in1,v1,x*y.nz[0],DM.eye(2)*y.nz[0]),
          (in1,v1,x*y[0,0],DM.eye(2)*y[0,0]),
          (in1,v1,x[0],DM.eye(2)[0,:]),
          (in1,v1,(x**2)[0],horzcat(*[2*x[0],MX(1,1)])),
          (in1,v1,x[0]+x[1],DM.ones(1,2)),
          (in1,v1,sin(repmat(x**2,1,3)),repmat(cos(c.diag(x**2))*2*c.diag(x),3,1)),
          (in1,v1,sin(repsum((x**2).T,1,2)),cos(x[0]**2+x[1]**2)*2*x.T),
          (in1,v1,vertcat(*[x[1],x[0]]),sparsify(DM([[0,1],[1,0]]))),
          (in1,v1,vertsplit(x,[0,1,2])[1],sparsify(DM([[0,1]]))),
          (in1,v1,vertcat(*[x[1]**2,x[0]**2]),blockcat([[MX(1,1),2*x[1]],[2*x[0],MX(1,1)]])),
          (in1,v1,vertsplit(x**2,[0,1,2])[1],blockcat([[MX(1,1),2*x[1]]])),
          (in1,v1,vertsplit(x**2,[0,1,2])[1]**3,blockcat([[MX(1,1),6*x[1]**5]])),
          (in1,v1,horzcat(*[x[1],x[0]]).T,sparsify(DM([[0,1],[1,0]]))),
          (in1,v1,horzcat(*[x[1]**2,x[0]**2]).T,blockcat([[MX(1,1),2*x[1]],[2*x[0],MX(1,1)]])),
          (in1,v1,diagcat(*[x[1]**2,y,x[0]**2]),
            blockcat(  [[MX(1,1),2*x[1]]] + ([[MX(1,1),MX(1,1)]]*14)  + [[2*x[0],MX(1,1)]] )
          ),
          (in1,v1,horzcat(*[x[1]**2,x[0]**2]).T,blockcat([[MX(1,1),2*x[1]],[2*x[0],MX(1,1)]])),
          (in1,v1,x[[0,1]],sparsify(DM([[1,0],[0,1]]))),
          (in1,v1,(x**2)[[0,1]],2*c.diag(x)),
          (in1,v1,x[[0,0,1,1]],sparsify(DM([[1,0],[1,0],[0,1],[0,1]]))),
          (in1,v1,(x**2)[[0,0,1,1]],blockcat([[2*x[0],MX(1,1)],[2*x[0],MX(1,1)],[MX(1,1),2*x[1]],[MX(1,1),2*x[1]]])),
          (in1,v1,wwr,sparsify(DM([[2,0],[0,2]]))),
          (in1,v1,x[[1,0]],sparsify(DM([[0,1],[1,0]]))),
          (in1,v1,x[[1,0],0],sparsify(DM([[0,1],[1,0]]))),
          (in1,v1,w,sparsify(DM([[1,0],[0,2]]))),
          (in1,v1,w2,blockcat([[1,MX(1,1)],[x[1],x[0]]])),
          (in1,v1,ww,2*c.diag(x)),
          (in1,v1,wwf,vertcat(*[x[[1,0]].T,x[[1,0]].T])),
          (in1,v1,yy[:,0],DM.eye(2)),
          (in1,v1,yy2[:,0],2*c.diag(x)),
          (in1,v1,yyy[:,0],sparsify(DM([[0,1],[1,0]]))),
          (in1,v1,mtimes(y,x),y),
          (in1,v1,mtimes(x.T,y.T),y),
          (in1,v1,mac(y,x,DM.zeros(Sparsity.triplet(2,1,[1],[0]))),y[Sparsity.triplet(2,2,[1,1],[0,1])]),
          (in1,v1,mac(x.T,y.T,DM.zeros(Sparsity.triplet(2,1,[1],[0]).T)),y[Sparsity.triplet(2,2,[1,1],[0,1])]),
          (in1,v1,mtimes(y[Sparsity.triplet(2,2,[0,1,1],[0,0,1])],x),y[Sparsity.triplet(2,2,[0,1,1],[0,0,1])]),
          (in1,v1,mtimes(x.T,y[Sparsity.triplet(2,2,[0,1,1],[0,0,1])].T),y[Sparsity.triplet(2,2,[0,1,1],[0,0,1])]),
          (in1,v1,mtimes(y,x**2),y*2*vertcat(*[x.T,x.T])),
          (in1,v1,sin(x),c.diag(cos(x))),
          (in1,v1,sin(x**2),c.diag(cos(x**2)*2*x)),
          (in1,v1,x*y[:,0],c.diag(y[:,0])),
          (in1,v1,x*y.nz[[0,1]],c.diag(y.nz[[0,1]])),
          (in1,v1,x*y.nz[[1,0]],c.diag(y.nz[[1,0]])),
          (in1,v1,x*y[[0,1],0],c.diag(y[[0,1],0])),
          (in1,v1,x*y[[1,0],0],c.diag(y[[1,0],0])),
          (in1,v1,c.dot(x,x),(2*x).T),
          (in1,v1,c.dot(x**2,x),(3*x**2).T),
          #(in1,v1,c.det(horzcat(*[x,DM([1,2])])),DM([-1,2])), not implemented
          (in1,v1,f1.call(in1)[1],y),
          (in1,v1,f1.call([x**2,y])[1],y*2*vertcat(*[x.T,x.T])),
          (in1,v1,f2.call(in1)[0],DM.zeros(0,2)),
          (in1,v1,f2(x**2,y),DM.zeros(0,2)),
          (in1,v1,f3.call(in1)[0],DM.zeros(0,2)),
          (in1,v1,f3.call([x**2,y])[0],DM.zeros(0,2)),
          (in1,v1,f4.call(in1)[0],DM.zeros(0,2)),
          (in1,v1,f4.call([x**2,y])[0],DM.zeros(0,2)),
          #(in1,v1,f1([x**2,[]])[1],DM.zeros(2,2)),
          #(in1,v1,f1([[],y])[1],DM.zeros(2,2)),
          (in1,v1,vertcat(*[x,DM(0,1)]),DM.eye(2)),
          (in1,v1,project(x**2, sparsify(DM([0,1])).sparsity()),blockcat([[MX(1,1),MX(1,1)],[MX(1,1),2*x[1]]])),
          (in1,v1,c.dot(x,y[:,0]),y[:,0].T),
          (in1,v1,x.nz[IM([[1,0]])].T*y.nz[IM([[0,2]])],blockcat([[MX(1,1),y.nz[0]],[y.nz[2],MX(1,1)]])),
          (in1,v1,x.nz[c.diag([1,0])]*y.nz[c.diag([0,2])],blockcat([[MX(1,1),y.nz[0]],[MX(1,1),MX(1,1)],[MX(1,1),MX(1,1)],[y.nz[2],MX(1,1)]])),
     ]:
      print(out)
      fun = Function("fun", inputs,[out,jac])
      funsx = fun.expand("expand_fun")
      fun_ad = [Function("fun", inputs,[out,jac], {'ad_weight':w, 'ad_weight_sp':w}) for w in [0,1]]
      funsx_ad = [f.expand('expand_'+f.name()) for f in fun_ad]

      fun_out = fun.call(values)
      funsx_out = funsx.call(values)

      self.checkarray(fun_out[0],funsx_out[0])
      self.checkarray(fun_out[1],funsx_out[1])

      self.check_codegen(fun,inputs=values)

      J_ = fun_out[1]

      def vec(l):
        ret = []
        for i in l:
          ret.extend(i)
        return ret

      storage2 = {}
      storage = {}

      vf_mx = None

      for f in [fun, fun.expand('expand_'+fun.name())]:
        d1 = f.forward(ndir)
        d2 = f.reverse(ndir)

        num_in = f.n_in()
        num_out = f.n_out()

        # evalThings
        for sym in [MX.sym, SX.sym]:
          if f.is_a('MXFunction') and sym==SX.sym: continue
          if f.is_a('SXFunction') and sym==MX.sym: continue

          # dense
          for spmod,spmod2 in itertools.product(spmods,repeat=2):
            fseeds = [[sym("f",spmod(f.sparsity_in(i))) for i in range(f.n_in())]  for d in range(ndir)]
            aseeds = [[sym("a",spmod2(f.sparsity_out(i)))  for i in range(f.n_out())] for d in range(ndir)]
            inputss = [sym("i",f.sparsity_in(i)) for i in range(f.n_in())]

            res = f.call(inputss,True)
            fwdsens = forward(res,inputss,fseeds,dict(always_inline=True))
            adjsens = reverse(res,inputss,aseeds,dict(always_inline=True))

            fseed = [DM(fseeds[d][0].sparsity(),random.random(fseeds[d][0].nnz())) for d in range(ndir) ]
            aseed = [DM(aseeds[d][0].sparsity(),random.random(aseeds[d][0].nnz())) for d in range(ndir) ]
            vf = Function("vf", inputss+vec([fseeds[i]+aseeds[i] for i in range(ndir)]),list(res) + vec([list(fwdsens[i])+list(adjsens[i]) for i in range(ndir)]))

            vf_in = list(values)
            offset = len(inputss)

            for d in range(ndir):
              vf_in.append(fseed[d])
              for i in range(len(values)-1):
                vf_in.append(0)

              vf_in.append(aseed[d])
              vf_in.append(0)

            vf_out = vf.call(vf_in)
            self.check_codegen(vf,inputs=vf_in)

            offset = len(res)
            for d in range(ndir):
              seed = array(fseed[d].T).ravel()
              sens = array(vf_out[offset+0].T).ravel()
              offset+=len(inputss)
              self.checkarray(sens,mtimes(J_,seed),"eval Fwd %d %s" % (d,str(type(f))+str(sym)))

              seed = array(aseed[d].T).ravel()
              sens = array(vf_out[offset+0].T).ravel()
              offset+=len(inputss)

              self.checkarray(sens,mtimes(J_.T,seed),"eval Adj %d %s" % (d,str([vf_out[i] for i in range(vf.n_out())])))


            assert(offset==vf.n_out())

            # Complete random seeding
            random.seed(1)
            vf_in = []
            for i in range(vf.n_in()):
              vf_in.append(DM(vf.sparsity_in(i),random.random(vf.nnz_in(i))))

            vf_out = vf.call(vf_in)
            self.check_codegen(vf,inputs=vf_in)
            storagekey = (spmod,spmod2)
            if not(storagekey in storage):
              storage[storagekey] = []
            storage[storagekey].append(vf_out)

            # Added to make sure that the same seeds are used for SX and MX
            if sym is MX.sym:
              vf_mx = vf

          # Second order sensitivities
          for sym2 in [MX.sym, SX.sym]:

            if vf.is_a('MXFunction') and sym2==SX.sym: continue
            if vf.is_a('MXFunction') and sym2==MX.sym: continue

            for spmod_2,spmod2_2 in itertools.product(spmods,repeat=2):
              fseeds2 = [[sym2("f",vf_mx.sparsity_in(i)) for i in range(vf.n_in())] for d in range(ndir)]
              aseeds2 = [[sym2("a",vf_mx.sparsity_out(i))  for i in range(vf.n_out()) ] for d in range(ndir)]
              inputss2 = [sym2("i",vf_mx.sparsity_in(i)) for i in range(vf.n_in())]

              res2 = vf.call(inputss2,True)
              fwdsens2 = forward(res2,inputss2,fseeds2,dict(always_inline=True))
              adjsens2 = reverse(res2,inputss2,aseeds2,dict(always_inline=True))

              vf2 = Function("vf2", inputss2+vec([fseeds2[i]+aseeds2[i] for i in range(ndir)]),list(res2) + vec([list(fwdsens2[i])+list(adjsens2[i]) for i in range(ndir)]))

              random.seed(1)
              vf2_in = []
              for i in range(vf2.n_in()):
                vf2_in.append(DM(vf2.sparsity_in(i),random.random(vf2.nnz_in(i))))

              vf2_out = vf2.call(vf2_in)
              self.check_codegen(vf2,inputs=vf2_in)
              storagekey = (spmod,spmod2)
              if not(storagekey in storage2):
                storage2[storagekey] = []
              storage2[storagekey].append(vf2_out)

      # Remainder of eval testing
      for store,order in [(storage,"first-order"),(storage2,"second-order")]:
        for stk,st in list(store.items()):
          for i in range(len(st)-1):
            for k,(a,b) in enumerate(zip(st[0],st[i+1])):
              if b.numel()==0 and sparsify(a).nnz()==0: continue
              if a.numel()==0 and sparsify(b).nnz()==0: continue
              self.checkarray(sparsify(a),sparsify(b),("%s, output(%d)" % (order,k)))

      for expand in [False, True]:
        #  jacobian()
        for mode in ["forward","reverse"]:
          ind = 0 if mode=='forward' else 1
          f = fun_ad[ind] if expand  else funsx_ad[ind]

          Jf=f.jacobian_old(0,0)
          Jf_out = Jf.call(values)

          self.check_codegen(Jf,inputs=values)
          self.checkarray(Jf_out[0],J_)
          self.checkarray(DM.ones(Jf.sparsity_out(0)),DM.ones(J_.sparsity()),str(out)+str(mode))
          self.checkarray(DM.ones(f.sparsity_jac(0, 0)),DM.ones(J_.sparsity()))

      # Scalarized
      if out.is_empty(): continue
      s_i  = out.sparsity().row()[0]
      s_j  = out.sparsity().get_col()[0]
      s_k = s_i*out.size2()+s_j
      H_ = None

      for expand in [False, True]:
        for mode in ["forward","reverse"]:
          w = 0 if mode=='forward' else 1
          f = Function("fun", inputs,[out[s_i,s_j],jac[s_k,:].T], {'ad_weight':w, 'ad_weight_sp':w})
          if expand: f=f.expand('expand_'+f.name())
          f_out = f.call(values)
          J_ = f_out[1]

          Hf=f.hessian_old(0, 0)
          Hf_out = Hf.call(values)
          self.check_codegen(Hf,inputs=values)
          if H_ is None:
            H_ = Hf_out[0]
          self.checkarray(Hf_out[0],H_,failmessage=("mode: %s" % mode))

if __name__ == '__main__':
    unittest.main()
