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
from numpy import *
import unittest
from types import *
from helpers import *
from casadi.tools import *

class Toolstests(casadiTestCase):

  def test_structure(self):
  
  
    s = struct(['x','y','z'])
    
    print s
  
    with self.assertRaises(Exception):
      struct_symSX(['x','x','z'])
    with self.assertRaises(Exception): 
      struct_symSX(['x','y','z','y'])
      

    
    s = struct_symSX([('x','y'),'z'])
    self.assertEqual(s.size,3)
    
    init = s(0)
    
    init['x'] = 5
    init['y'] = 7
    init['z'] = 8
    
    x,y,z = s[...]
    
    self.checkarray(init.cat,DMatrix([5,7,8]))
    
    init[...] = [1,2,3]
    
    self.checkarray(init.cat,DMatrix([1,2,3]))
    
    init = s()
    init['x'] = 5
    init['y'] = 7
    init['z'] = 8
    
    self.checkarray(init.cat,DMatrix([5,7,8]))
    
    init[{}] = {'x': 12}
    self.checkarray(init.cat,DMatrix([12,7,8]))
    
    init = s(inf)
    self.checkarray(init.cat,DMatrix([inf,inf,inf]))
    
    f = DMatrix.zeros(3,1)
    init = s(f)
    f[0] = 5
    self.checkarray(init.cat,DMatrix([5,0,0]))
    
    s = struct_symSX(['x','y','z'],order=['x','y','z'])
    with self.assertRaises(Exception): 
      s = struct_symSX([('x','y'),'z'],order=['x','y','z'])
    
    s = struct_symSX(['x','y','z'])
    self.assertEqual(s.size,3)
    self.assertEqual(s[vertcat,['x','y']].shape[0],2)
    self.assertTrue(isinstance(s[list,vertcat,['x','y']],list))
    
    s = struct_symSX([entry('x',repeat=5),entry('y',repeat=6),entry('z')])
    self.assertEqual(s.size,12)
    self.assertEqual(len(s["x"]),5)
    self.assertEqual(len(s["y"]),6)
    self.assertTrue(s.cat.at(1).getName().startswith("x"))
    s = struct_symSX([(entry('x',repeat=5),entry('y',repeat=[6,5])),entry('z')])
    
        
    self.checkarray(s.f["x"],[0,6,12,18,24])
    self.checkarray(s.f["z"],[35])
    self.checkarray(s.f["y"],[1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,19,20,21,22,23,25,26,27,28,29,30,31,32,33,34])
    
    self.assertEqual(s.size,36)
    self.assertEqual(len(s["x"]),5)
    self.assertEqual(len(s["y"]),6)
    self.assertTrue(s.cat.at(1).getName().startswith("y"))
    s = struct_symSX([entry("x",shape=(3,2)),entry("y",shape=2),entry("z",shape=Sparsity.dense(3)),entry("w",shape=Sparsity.triu(5))])
    self.assertEqual(s.size,6+2+3+15)
    self.assertTrue(s["x"].sparsity()==Sparsity.dense(3,2))
    self.assertTrue(s["y"].sparsity()==Sparsity.dense(2,1))
    self.assertTrue(s["z"].sparsity()==Sparsity.dense(3,1))
    self.assertTrue(s["w"].sparsity()==Sparsity.triu(5))
    
    x  = SX.sym("x",2)
    x2 = SX.sym("x2",2)
    s = struct_symSX([entry('a',sym=x),'y','z'])
    self.assertTrue(s.cat.at(0).getName().startswith("x"))
    self.assertEqual(s.size,4)
    with self.assertRaises(Exception):
      struct_symSX([entry('a',sym=x+x),'y','z'])
    with self.assertRaises(Exception):
      struct_symSX([entry('a',sym=[x+x]),'y','z'])
    s = struct_symSX([entry('a',sym=vertcat([x,x2])),'y','z'])
    self.assertEqual(s.size,6)
    with self.assertRaises(Exception):
      s = struct_symSX([entry('a',repeat=6,sym=x),'y','z'])
    
    s = struct_symSX(['x','y','z'])
 
    S = struct_symSX([entry("X",sym=s)])
    self.assertTrue(S.cat.at(0).getName()=="x")
    S = struct_symSX([entry("X",struct=s)])
    self.assertTrue(S.cat.at(0).getName()=="X_x")
    S = struct_symSX([entry("X",repeat=[5],struct=s)])
    self.assertEqual(S.size,15)
    
    s0 = struct_symSX(['x','y','z'])
    s1 = struct_symSX(['x','y'])
    S0 = struct_symSX([entry("X",struct=s0),entry("Y",repeat=[5],struct=s0),entry("Z",struct=s0)])
    S1 = struct_symSX([entry("X",struct=s1),entry("Y",repeat=[5],struct=s1)])
    num0 = S0(0)
    num1 = S1(1)

    num0[nesteddict] = num1[nesteddict]
    
    S = struct_symSX([entry("P",shapestruct=(s0,s0))])
 
    num = S(0)
    num["P"] = DMatrix([[1,2,3],[4,5,6],[7,8,9]])
    self.checkarray(num["P",index["x"],index["y"]],DMatrix([2]))
    self.checkarray(num["P",index["x"],:],DMatrix([1,2,3]).T)
    self.checkarray(num["P",:,index["y"]],DMatrix([2,5,8]))
    self.checkarray(num["P",:,:],DMatrix([[1,2,3],[4,5,6],[7,8,9]]))
    
    self.checkarray(num["P",indexf[["x","y"]],indexf[["z","x"]]],DMatrix([[3,1],[6,4]]))
    
    self.checkarray(num["P",index[list,vertcat,["x","y"]],index[list,vertcat,["z","x"]]],DMatrix([[3,1],[6,4]]))
    self.checkarray(num["P",index[vertcat,["x","y"]],index[vertcat,["z","x"]]],DMatrix([3,4]))
    
    S = struct_symSX([entry("P",shapestruct=s0)])
    
    num = S(0)
    num["P"] = DMatrix([1,2,3])
    self.checkarray(num["P",index["x"]],DMatrix([1]))
    self.checkarray(num["P",:],DMatrix([1,2,3]))
    
    self.checkarray(num["P",index["x"],0],DMatrix([1]))
    self.checkarray(num["P",:,0],DMatrix([1,2,3]))
    
    with self.assertRaises(Exception):
      num["P",:,index["x"]]
    
    
    S = struct_symSX([entry("P",shapestruct=(3,s0))])
 
    num = S(0)
    num["P"] = DMatrix([[1,2,3],[4,5,6],[7,8,9]])
    with self.assertRaises(Exception):
      self.checkarray(num["P",index["x"],index["y"]],DMatrix([2]))
    with self.assertRaises(Exception):
      self.checkarray(num["P",index["x"],:],DMatrix([1,2,3]).T)
    self.checkarray(num["P",:,index["y"]],DMatrix([2,5,8]))
    self.checkarray(num["P",:,:],DMatrix([[1,2,3],[4,5,6],[7,8,9]]))
    
    self.checkarray(num["P",:,indexf[["z","x"]]],DMatrix([[3,1],[6,4],[9,7]]))
    
    self.checkarray(num["P",:,index[list,vertcat,["z","x"]]],DMatrix([[3,1],[6,4],[9,7]]))
    self.checkarray(num["P",:,index[vertcat,["z","x"]]],DMatrix([3,1,6,4,9,7]))
    
    S = struct_symSX([entry("P",shapestruct=s0)])
 
    
    
    s0 = struct_symSX(['x','y',entry('q',shape=4),'z'])

    S = struct_symSX([entry("P",shapestruct=(s0,s0))])
    
    num = S(0)
    num["P",indexf[["x","q"]],indexf[["x","q"]]] = 1
    
    self.checkarray(num["P"][[1,6],[1,6]],DMatrix.zeros(2,2))
    self.checkarray(num["P"][[0,2,3,4,5],[0,2,3,4,5]],DMatrix.ones(5,5))

    num = S(0)
    num["P",["x","q"],["x","q"]] = 1
    
    self.checkarray(num["P"][[1,6],[1,6]],DMatrix.zeros(2,2))
    self.checkarray(num["P"][[0,2,3,4,5],[0,2,3,4,5]],DMatrix.ones(5,5))
    
    num = S(0)
    num["P",["x","q"],["x","q"]] = DMatrix.ones(5,5)
    
    self.checkarray(num["P"][[1,6],[1,6]],DMatrix.zeros(2,2))
    self.checkarray(num["P"][[0,2,3,4,5],[0,2,3,4,5]],DMatrix.ones(5,5))
    
    with self.assertRaises(Exception):
      struct_symMX(['x','x','z'])
    with self.assertRaises(Exception): 
      struct_symMX(['x','y','z','y'])
    
    s = struct_symMX([('x','y'),'z'])
    s = struct_symMX(['x','y','z'],order=['x','y','z'])
    with self.assertRaises(Exception): 
      s = struct_symMX([('x','y'),'z'],order=['x','y','z'])
    
    s = struct_symMX(['x','y','z'])
    self.assertEqual(s.size,3)
    s = struct_symMX([entry('x',repeat=5),entry('y',repeat=6),entry('z')])
    self.assertEqual(s.size,12)
    self.assertEqual(len(s["x"]),5)
    self.assertEqual(len(s["y"]),6)
    s = struct_symMX([(entry('x',repeat=5),entry('y',repeat=[6,5])),entry('z')])
    self.assertEqual(s.size,36)
    self.assertEqual(len(s["x"]),5)
    self.assertEqual(len(s["y"]),6)
   
    
    s = struct_symMX([entry("x",shape=(3,2)),entry("y",shape=2),entry("z",shape=Sparsity.dense(3)),entry("w",shape=Sparsity.triu(5))])
    self.assertEqual(s.size,6+2+3+15)
    self.assertTrue(s["x"].sparsity()==Sparsity.dense(3,2))
    self.assertTrue(s["y"].sparsity()==Sparsity.dense(2,1))
    self.assertTrue(s["z"].sparsity()==Sparsity.dense(3,1))
    self.assertTrue(s["w"].sparsity()==Sparsity.triu(5))
    
    x  = MX.sym("x",2)
    x2 = MX.sym("x2",2)
    with self.assertRaises(Exception):
      s = struct_symMX([entry('x',sym=x),'y','z'])
    with self.assertRaises(Exception):
      struct_symMX([entry('x',sym=x+x),'y','z'])
    with self.assertRaises(Exception):
      struct_symMX([entry('x',sym=[x+x]),'y','z'])
    with self.assertRaises(Exception):
      struct_symMX([entry('x',sym=vertcat([x,x2])),'y','z'])
    with self.assertRaises(Exception):
      s = struct_,MX.sym([(2,'x',[x,x2]),'y','z'])
    
    s = struct_symMX(['x','y','z'])
    S = struct_symMX([entry("X",struct=s)])
    S = struct_symMX([entry("X",repeat=[5],struct=s)])
    
    x = SX.sym("x",2)
    y0 = sin(x) 
    y1 = cos(x)
    
    v = struct_SX([
      entry("x",expr=x),
      entry("y",expr=[y0,y1]),
      entry("w",expr=[[x,x**2],[x**3,x**4]])
    ])
    
    
    x = MX.sym("x",2)
    m = struct_symMX(['a','b'])
    y0 = sin(x) 
    y1 = cos(x)
    
    V = struct_MX([
      entry("x",expr=x),
      (entry("y",expr=[y0,y1],struct=m),
      entry("w",expr=[[x,x**2],[x**3,x**4]],struct=m)
      )
    ])
    self.assertEqual(V.size,14)

    self.assertTrue(isinstance(V.cat,MX)) 
    def isEqualV(a,b):
      ft = MXFunction([x,m],[a-b])
      ft.init()
      for i in range(ft.getNumInputs()):
        ft.setInput(numpy.random.rand(*ft.getInput(i).shape),i)
      ft.evaluate()
      self.checkarray(ft.getOutput(),DMatrix.zeros(*ft.getOutput().shape))
    
    isEqualV(V["x"],x)
    isEqualV(V["y",0],y0)
    isEqualV(V["y",1],y1)
    self.assertEqual(V["y",0,'a'].shape,(1,1))
    
    with self.assertRaises(Exception):
      V["y",0] = MX.sym("x",4) # shape mismatch
    abc = MX.sym("abc",2)
    V["y",0] = abc
    
    def isEqualV(a,b):
      ft = MXFunction([x,m,abc],[a-b])
      ft.init()
      for i in range(ft.getNumInputs()):
        ft.setInput(numpy.random.rand(*ft.getInput(i).shape),i)
      ft.evaluate()
      self.checkarray(ft.getOutput(),DMatrix.zeros(*ft.getOutput().shape))
      
    isEqualV(V["y",0],abc)

    states = struct_symSX([
                entry('x'),
                entry('y'),
                entry('z'),
                entry('u',shape=Sparsity.dense(4)),
                entry('v',repeat=[4,2]),
                entry('w',repeat=[6]),
                entry('p',repeat=[9],shape=Sparsity.dense(6))
             ],order=['x','y','z','u',('v','w'),'p'])
             
    shooting = struct_symSX([entry('X',struct=states,repeat=[4,5]),entry('U',repeat=[3])],order=[('X','U')])


    
    self.assertEqual(shooting.size,1503)
    s = shooting["X",:,:,{}]
    self.assertTrue(isinstance(s,list))
    self.assertEqual(len(s),4)
    self.assertTrue(isinstance(s[0],list))
    self.assertEqual(len(s[0]),5)
    self.assertTrue(isinstance(s[0][0],dict))
    self.assertTrue('x' in s[0][0])
    self.assertEqual(len(s[0][0]),7)
    if CasadiOptions.getSimplificationOnTheFly():
      self.assertTrue(isEqual(s[0][0]["x"],shooting["X",0,0,"x"]))
    
    
    init = shooting(nan)
    
    init['X',0,-1,'p',0] = 2
    self.checkarray(init.cat[shooting.i['X',0,-1,'p',0]],DMatrix.ones(6)*2)
    self.assertEqual(sum([i!=i for i in init.cat.data()]),1503-6)
    
    init['X',0,-1,'p',0] = [3]*6
    self.checkarray(init.cat[shooting.i['X',0,-1,'p',0]],DMatrix.ones(6)*3)
    self.assertEqual(sum([i!=i for i in init.cat.data()]),1503-6)
 
    init['X',0,-1,'p',0] = DMatrix([4]*6)
    self.checkarray(init.cat[shooting.i['X',0,-1,'p',0]],DMatrix.ones(6)*4)
    self.assertEqual(sum([i!=i for i in init.cat.data()]),1503-6)
    
    init['X',0,-1,'p',:] = 7
    
    self.checkarray(init.cat[shooting.i['X',0,-1,'p',1]],DMatrix.ones(6)*7)
    self.assertEqual(sum([i!=i for i in init.cat.data()]),1503-6*9)
    
    self.checkarray(init['X',0,-1,'p',horzcat,:],DMatrix.ones(6,9)*7)

    with self.assertRaises(Exception):
      init['X',0,-1,'p',:] = [1,2,3,4,5,6]
      
    init['X',0,-1,'p',:] = repeated([1,2,3,4,5,6])
    self.checkarray(init['X',0,-1,'p',vertcat,:,2],DMatrix.ones(9)*3)
    
    init = shooting(DMatrix(range(shooting.size)))

    self.checkarray(init['X',vertcat,:,horzcat,:],init['X',blockcat,:,:])

    init = shooting(nan)

    init['X',:,:,['x','y']] = repeated(repeated([6,5]))

    print init.cat

    init['X',:,:,{}] = repeated(repeated({'x': 9,'y': 3}))

    print init.cat
    
    V = struct_SX(shooting)
    
    V['X',:,:,['x','y']] = repeated(repeated([6,5]))
    
    print V
    
    V = struct_symMX(shooting)
    
    print V
    
    V = struct_MX(shooting)
    
    print V
    
    
    init = shooting(nan)
    

    
    init['X',0,-1,'p',horzcat,:] = DMatrix.ones(6,9)*2
    
    self.assertEqual(sum(init.cat==2),6*9)

    init['X',0,-1,'p',vertcat,:,2] = DMatrix.ones(9)*3
    
    self.assertEqual(sum(init.cat==2),6*9-9)
    self.assertEqual(sum(init.cat==3),9)

    self.checkarray(init['X',0,-1,'p',horzcat,:,[2,3]],vertcat([DMatrix.ones(1,9)*3,DMatrix.ones(1,9)*2]))
    
    init['X',:,-1,'p',horzcat,:] = repeated(DMatrix.ones(6,9)*2)
    
    self.assertEqual(sum(init.cat==2),4*6*9)
    
    init['X',0,0,'v',blockcat] = DMatrix.ones(4,2)*7
    
    self.assertEqual(sum(init.cat==7),4*2)

    init['X',:,0,'v',blockcat] = repeated(DMatrix.ones(4,2)*7)
    
    self.assertEqual(sum(init.cat==7),4*2*4)
    
    init['X',:,:,'v',blockcat] = repeated(DMatrix.ones(4,2)*7)
    self.assertEqual(sum(init.cat==7),4*2*4*5)
    
    init['X',0,0,'p',:] = range(9)
    
    
    print index["a"]

    init = shooting(range(shooting.size))
    for i in range(shooting.size):
      ci = shooting.getCanonicalIndex(i)
      self.assertEqual(i, init.__getitem__(ci))
      self.assertTrue("_".join(map(str,ci)).startswith(shooting.cat.at(i).getName()))
      
  def test_structure_prefix(self):
    self.message("structure prefix")
    s = struct(["x","y","z"])

    S = struct_symSX([entry("X",repeat=12,struct=s)])

    print S.__class__
    print S.prefix

    a = S.prefix["X"]

    num = S()

    init = num.prefix["X",-1]

    init["x"] = 12
    
    self.assertEqual(num["X",-1,"x"],init["x"])

    self.checkarray(num.cat,DMatrix([0]*(S.size-3) + [12,0,0]))

  def test_structure_repeated_dmatrix(self):
    self.message("repeated dmatrix")
    s = struct(["x","y","z"])
    d = DMatrix.zeros(s.size,12)
    a = s.repeated(d)
    
    a[:,"x"] = range(12)
    
    self.checkarray(a[4,"x"],DMatrix([4]))
    self.checkarray(d,vertcat([DMatrix(range(12)).T,DMatrix.zeros(1,12),DMatrix.zeros(1,12)]))

  def test_structure_squared_dmatrix(self):
    self.message("squared dmatrix")
    s = struct(["x","y","z"])
    d = DMatrix.zeros(3,3)
    a = s.squared(d)
    
    a["x","y"] = 2
    a["y","x"] = 1
    
    self.checkarray(a["x","y"],DMatrix([2]))
    self.checkarray(a["y","x"],DMatrix([1]))
    
    self.checkarray(d,DMatrix([[0,2,0],[1,0,0],[0,0,0]]))

  def test_structure_squared_repeated_dmatrix(self):
    self.message("squared repeated dmatrix")
    s = struct(["x","y","z"])
    
    d = DMatrix([[0,2,0,0,0,0,0,2,0],[1,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0]])
    a = s.squared_repeated(d)
    self.checkarray(a[0,"x","y"],DMatrix([2]))
    self.checkarray(a[0,"y","x"],DMatrix([1]))
    self.checkarray(a[-1,"x","y"],DMatrix([2]))
    self.checkarray(a[-1,"y","x"],DMatrix([1]))
    
    d = DMatrix.zeros(3,9)
    a = s.squared_repeated(d)
    
    a[0,"x","y"] = 2
    a[0,"y","x"] = 1
    a[-1,"x","y"] = 2
    a[-1,"y","x"] = 1
    self.checkarray(d,DMatrix([[0,2,0,0,0,0,0,2,0],[1,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0]]))

  def test_typemaps(self):
    self.message("typemaps")
    s = struct(["x","y","z"])
    d = DMatrix.zeros(3,3)
    a = s.squared(d)
    
    print type(a)
    print sin(a)
    print a+1
    
  def test_sparse(self):
    a = struct_symSX([entry("a",shape=Sparsity.diag(5))])
    b = struct_symMX([(entry("b",struct=a))])

    self.checkarray(b["b"].shape,(5,1))
    self.checkarray(b["b","a"].shape,(5,5))
    
  def test_symm(self):
    a = struct_symSX([entry("P",shape=(3,3),type='symm')])
    
    b = a()
    b["P"] = DMatrix([[0,3,6],[1,4,7],[2,5,8]])
    
    self.assertEqual(a.size,6)
    self.checkarray(a["P"].shape,(3,3))
    
    f = SXFunction([a],[a["P"]])
    f.init()
    f.setInput(range(6))
    f.evaluate()
    
    self.checkarray(f.getOutput(),DMatrix([[0,1,3],[1,2,4],[3,4,5]]))
    self.checkarray(b["P"],DMatrix([[0,3,6],[3,4,7],[6,7,8]]))
    self.checkarray(b.cat,DMatrix([0,3,4,6,7,8]))
    
    states = struct_symSX(["x","y","z"])

    a = struct_symSX([entry("P",shapestruct=(states,states),type='symm')])

    b = a()
    b["P"] = DMatrix([[0,3,6],[1,4,7],[2,5,8]])

    self.checkarray(b["P",["x","z"],["x","z"]] ,DMatrix([[0,6],[6,8]]))

    b["P","y","x"] = 66
    self.checkarray(b["P"],DMatrix([[0,66,6],[66,4,7],[6,7,8]]))

    
    b["P","x","y"] = 666
    
    self.checkarray(b["P"],DMatrix([[0,666,6],[666,4,7],[6,7,8]]))

    b["P",["x"],["y","z"]] = DMatrix([1,12]).T
    self.checkarray(b["P"],DMatrix([[0,1,12],[1,4,7],[12,7,8]]))
    
    b["P",["x","z"],["x","z"]] = DMatrix([[11,0],[0,11]])
    self.checkarray(b["P"],DMatrix([[11,1,0],[1,4,7],[0,7,11]]))
    
    a = struct_symSX([entry("P",shape=Sparsity.banded(3,1),type='symm')])

    b = a()
    with self.assertRaises(Exception):
      b["P"] = DMatrix([[11,1,2],[1,4,5],[2,5,8]])
      
    with self.assertRaises(Exception):
      b["P"] = sparse(DMatrix([[11,0,1],[1,4,5],[0,5,8]]))
    
    b["P"] = sparse(DMatrix([[11,1,0],[1,4,5],[0,5,8]]))
    
    self.checkarray(b["P"],DMatrix([[11,1,0],[1,4,5],[0,5,8]]))
    
    with self.assertRaises(Exception):
      b["P",:,0] = DMatrix([1,2,3])
    
    with self.assertRaises(Exception):
      b["P",0,:] = DMatrix([1,2,3]).T
      
    b["P",:,0] = sparse(DMatrix([1,2,0]))
    self.checkarray(b["P"],DMatrix([[1,2,0],[2,4,5],[0,5,8]]))

    b["P",0,:] = sparse(DMatrix([11,12,0])).T
    self.checkarray(b["P"],DMatrix([[11,12,0],[12,4,5],[0,5,8]]))
    
  def test_callableExtraIndex(self):
    a = struct_symSX([entry("a",shape=(5,3)),entry("b",shape=(4,3))])
    b = a()
    
    b["a",vec] = range(15)
    self.checkarray(b.cat,DMatrix(range(15)+[0]*12))
    
    self.checkarray(b["a",vec],DMatrix(range(15)))
    
    b["a",vec] = range(15)

    self.checkarray(b["a",vec],DMatrix(range(15)))
    
  def test_pickling(self):
    import pickle

    x = struct([entry("a",shape=(2,3)),entry("b",shape=(2,3))])

    s = pickle.dumps(x)

    w = pickle.loads(s)

    self.assertEqual(x.size,w.size)
    x = struct([entry("a",shape=(2,3)),entry("b",shape=(2,3))])

    y = x()

    y["b",:,2] = 12

    s = pickle.dumps(y)
    w = pickle.loads(s)

    self.checkarray(w["b",:,2],DMatrix([12,12]))
    
  def test_numpyint(self):
    state = struct_symSX(['x', 'y'])
    x = struct_symSX([entry('states', struct=state, repeat=10)])
    x_init = x()
    x_init['states', 0, 'x'] # OK
    a = [1,2,3]
    x_init['states', int32(0), 'x'] 
    x_init['states', int64(0), 'x']
    
  def test_numpyint(self):
    s = struct_symSX(map(entry, 'xyz')) # OK 
    print s['x']
    s = struct_symSX(map(entry, u'xyz')) # IndexError: list index out of range
    print s[u'x']
    
  def test_pickling_null(self):
    import pickle
    s = struct_symMX([
      entry("a",shape=(2,3)),
      entry("b",shape=(0,0))
    ])

    tt = pickle.dumps(s)

    print pickle.loads(tt)
    
  def test_bug_structSXMX(self):
    n= 2
    x_sx = struct_symSX([
        entry("x",shape=n),
        entry("S",shape=(n,n))
    ])
    
    x_mx = struct_symSX([
        entry("x",shape=n),
        entry("S",shape=(n,n))
    ])
    
    X_sx = struct_SX(x_sx)
    X_sx["x"] = DMatrix(range(n))
    X_sx["S"] = DMatrix(range(n,n+n*n),n,n)
   
    X_mx = struct_MX(x_sx)
    X_mx["x"] = DMatrix(range(n))
    X_mx["S"] = DMatrix(range(n,n+n*n),n,n)
    
    self.checkarray(x_sx.struct.map[("S",)],DMatrix(range(n,n+n*n),n,n))
    self.checkarray(x_mx.struct.map[("S",)],DMatrix(range(n,n+n*n),n,n))
    self.checkarray(X_sx.cat,DMatrix(range(n+n*n)))
    self.checkarray(X_mx.cat,DMatrix(range(n+n*n)))
    
    for s, S in [(x_sx,struct_symSX),(x_mx,struct_symMX)]:
      h = S([entry("w",struct=s)])
      hX = struct_SX(h)
      hX["w","x"] = DMatrix(range(n))
      hX["w","S"] = DMatrix(range(n,n+n*n),n,n) 
      
      self.checkarray(h.struct.map[("w","S",)],DMatrix(range(n,n+n*n),n,n))
      self.checkarray(hX.cat,DMatrix(range(n+n*n)))
      
      self.checkarray(h.struct.map[("w",)],DMatrix(range(n+n*n)))
      self.checkarray(hX["w"],DMatrix(range(n+n*n)))
      
      hX["w"] = DMatrix(range(n+n*n))
      self.checkarray(hX.cat,DMatrix(range(n+n*n)))
    
    n= 2
    m = 3
    x_sx = struct_symSX([
        entry("x",shape=n),
        entry("S",shape=(n,m))
    ])
    
    x_mx = struct_symSX([
        entry("x",shape=n),
        entry("S",shape=(n,m))
    ])
    
    X_sx = struct_SX(x_sx)
    X_sx["x"] = DMatrix(range(n))
    X_sx["S"] = DMatrix(range(n,n+n*m),n,m)
   
    X_mx = struct_MX(x_sx)
    X_mx["x"] = DMatrix(range(n))
    X_mx["S"] = DMatrix(range(n,n+n*m),n,m)
    
    self.checkarray(x_sx.struct.map[("S",)],DMatrix(range(n,n+n*m),n,m))
    self.checkarray(x_mx.struct.map[("S",)],DMatrix(range(n,n+n*m),n,m))
    self.checkarray(X_sx.cat,DMatrix(range(n+n*m)))
    self.checkarray(X_mx.cat,DMatrix(range(n+n*m)))
    
    n = 3
    x_sx = struct_symSX([
        entry("x",shape=n),
        entry("S",shape=Sparsity.triu(n))
    ])
    
    x_mx = struct_symSX([
        entry("x",shape=n),
        entry("S",shape=Sparsity.triu(n))
    ])
    
    X_sx = struct_SX(x_sx)
    X_sx["x"] = DMatrix(range(n))
    X_sx["S"] = DMatrix(Sparsity.triu(n),range(n,n+n*(n+1)/2))
   
    X_mx = struct_MX(x_sx)
    X_mx["x"] = DMatrix(range(n))
    X_mx["S"] = DMatrix(Sparsity.triu(n),range(n,n+n*(n+1)/2))
    
    self.checkarray(x_sx.struct.map[("S",)],DMatrix(Sparsity.triu(n),range(n,n+n*(n+1)/2)))
    self.checkarray(x_mx.struct.map[("S",)],DMatrix(Sparsity.triu(n),range(n,n+n*(n+1)/2)))
    self.checkarray(X_sx.cat,DMatrix(range(n+n*(n+1)/2)))
    self.checkarray(X_mx.cat,DMatrix(range(n+n*(n+1)/2)))
    
  def test_MX_result(self):
    s = struct_symMX(["a",entry("b",shape=2),entry("c",shape=(2,2))])

    V = s.cat

    d = s(V)

    print d["b"]

    f = MXFunction([V],[d["a"],d["b"],d["c"]])
    f.init()

    s_ = s()
    s_["a"] = 1
    s_["b"] = DMatrix([2,3])
    s_["c"] = DMatrix([[4,5],[6,7]])
    f.setInput(s_)
    f.evaluate()
    
    self.checkarray(f.getOutput(0),s_["a"])
    self.checkarray(f.getOutput(1),s_["b"])
    self.checkarray(f.getOutput(2),s_["c"])
   
  def test_issue1116(self):
    S = struct([entry("A",shape=(4,4),type="symm")])
    v = S()
    v["A"]  = np.eye(4)
    
    a = DMatrix(v["A"])
    v["A"]  = DMatrix(np.eye(4))
    b = DMatrix(v["A"])
    
    self.checkarray(a,b)

  def test_jacobian(self):
    states = struct_symSX(["x","y"])
    controls = struct_symSX(["u","v","w"])

    #   u v w
    # x 1 3 5
    # y 2 4 6

    for t,F in [(SX,SXFunction), (MX,MXFunction)]:


      J = t.sym("J",states.size,controls.size)

      J = states.product(controls,J)

      f = F([J],[J["x","v"], J["x",:] , J["y",["v","w"]],  J[:,"u"] ])
      f.init()
      f.setInput(range(1,7))

      f.evaluate()

      self.checkarray(f.output(0),DMatrix([3]))
      self.checkarray(f.output(1),DMatrix([1,3,5]).T)
      self.checkarray(f.output(2),DMatrix([4,6]).T)
      self.checkarray(f.output(3),DMatrix([1,2]))
   
  def test_empty_bug(self):
    
    Params = struct_symMX([ ])
    params = Params()
    
    self.checkarray(params.shape,(0,1))

  def test_empty_expr_bug(self):
    
    eq = MX.sym("X")

    g = struct_MX([ entry( 'equality', expr = eq),
                    entry( 'inequality', expr = [] )   ])
    
    self.checkarray(g.shape,(1,1))
    
    self.assertTrue(len(g["inequality"])==0)
    
if __name__ == '__main__':
    unittest.main()

