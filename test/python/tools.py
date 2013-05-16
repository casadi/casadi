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
from casadi.tools import *

class Toolstests(casadiTestCase):

  def test_collection(self):
      self.message("Collection")

      p = Collection()
      p.x = ssym("x",2)
      p.z = ssym("z",2,4)
      p.y = ssym("y",3,2)

      p.freeze()
      
      self.assertTrue(isinstance(p[...],SXMatrix))
      self.assertEqual(p.size,16)
      self.assertEqual(p[...].shape[0],16)
      self.assertEqual(p[...].shape[1],1)

      self.assertTrue(isEqual(p.x,p["x"]))
      self.assertTrue(isEqual(p.y,p["y"]))
      self.assertTrue(isEqual(p.z,p["z"]))

      self.checkarray(p.i_x,IMatrix([0,1]),"")
      self.checkarray(p.i_z,IMatrix([[2,4,6,8],[3,5,7,9]]),"")
      self.checkarray(p.i_y,IMatrix([[10,13],[11,14],[12,15]]),"")
      
      self.checkarray(p.iv_x,list(IMatrix([0,1])),"")
      
      p = Collection()
      p.x = [ssym("x")]
      p.z = [ssym("z%d" % i) for i in range(2)]
      p.y = [[ssym("y%d%d"% (i,j)) for i in range(2)] for j in range(3)]
      p.setOrder(["x","y","z"])
      p.freeze()
      
      self.assertTrue(isinstance(p[...],SXMatrix))
      self.assertEqual(p.size,9)

      self.assertTrue(all(map(isEqual,p.x,p["x"])))
      self.assertTrue(all(map(isEqual,p.z,p["z"])))
      self.assertTrue(all(map(lambda a,b: all(map(isEqual,a,b)),p.y,p["y"])))
      self.checkarray(p[...],vertcat([p.x[0],p.y[0][0],p.y[0][1],p.y[1][0],p.y[1][1],p.y[2][0],p.y[2][1],p.z[0],p.z[1]]),"")

      self.checkarray(p.i_x[0],0,"")
      self.checkarray(p.i_z[0],7,"")
      self.checkarray(p.i_z[1],8,"")
      self.checkarray(p.i_y[0],[1,2],"")
      self.checkarray(p.i_y[1],[3,4],"")
      
      self.checkarray(p.i_y.__getitem__(0),[1,2],"")
      #self.checkarray(p.i_y.__getitem__((0,)),[1,2],"")

      self.checkarray(p.i_y,[[1,2],[3,4],[5,6]],"")
      
      #self.checkarray(p.i_y[:,0],[1,3,5],"")
      #self.checkarray(p.i_y[:,1],[2,4,6],"")
      
      #self.checkarray(p.i_y[0,:],[1,2],"")
      #self.checkarray(p.i_y[1,:],[3,4],"")
      
      #self.checkarray(p._i["y",0],[1,2],"")
      #self.checkarray(p._i["y",:,0],[1,3,5],"")
      
      p = Collection()
      p.x = [ssym("x")]
      p.z = [ssym("z%d" % i) for i in range(2)]
      p.y = [[ssym("y%d%d"% (i,j)) for i in range(2)] for j in range(3)]
      p.setOrder(["x",("y","z")])
      p.freeze()
      
      self.assertTrue(isinstance(p[...],SXMatrix))
      self.assertEqual(p.size,9)

      self.assertTrue(all(map(isEqual,p.x,p["x"])))
      self.assertTrue(all(map(isEqual,p.z,p["z"])))
      self.assertTrue(all(map(lambda a,b: all(map(isEqual,a,b)),p.y,p["y"])))
      self.checkarray(p[...],vertcat([p.x[0],p.y[0][0],p.y[0][1],p.z[0],p.y[1][0],p.y[1][1],p.z[1],p.y[2][0],p.y[2][1]]),"")

      self.checkarray(p.z[...],vertcat([p.z[0],p.z[1]]),"")
      self.checkarray(p.y[...],vertcat([p.y[0][0],p.y[0][1],p.y[1][0],p.y[1][1],p.y[2][0],p.y[2][1]]),"")
      
      self.checkarray(p.i_x[0],0,"")
      self.checkarray(p.i_z[0],3,"")
      self.checkarray(p.i_z[1],6,"")
      self.checkarray(p.i_z[...],IMatrix([3,6]),"")
      self.checkarray(p.i_y[0],[1,2],"")
      self.checkarray(p.i_y[1],[4,5],"")
      self.checkarray(p.i_y[2],[7,8],"")
      self.checkarray(p.i_y[...],IMatrix([1,2,4,5,7,8]),"")
    
      p = Collection()
      p.a = ssym("a")
      p.b = ssym("b")
      p.freeze()

      g = Collection()
      g.c = ssym("c")
      g.d = p
      g.e = ssym("e")
      g.freeze()
       
      self.assertEqual(g.size,4)
      self.checkarray(g[...],vertcat([g.c,p.a,p.b,g.e]),"")

      self.assertEqual(p.size,2)
      self.checkarray(p[...],vertcat([p.a,p.b]),"")
      
      self.checkarray(p.i_a,0)
      self.checkarray(p.i_b,1)

      self.checkarray(g.i_c,0)
      self.checkarray(g.i_d["a"],1)
      self.checkarray(g.i_d["b"],2)
      self.checkarray(g.i_e,3)
      
      self.checkarray(g.d[...],vertcat([p.a,p.b]),"")
      

  def test_variables(self):
      self.message("Variables")
             
      p = Variables()

      p.x = ssym("x",2)
      p.z = ssym("z",2,4)
      p.y = ssym("y",3,2)

      xother = Variables()
      xother.a = SX("x")
      xother.b = diag(ssym("[a,b]"))
      xother.freeze()
      
      p.xother = xother
      
      p.freeze()

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

      self.assertEqual(p.veccat().numel(),21)
      self.assertEqual(p.veccat().size(),19)
      self.assertEqual(p.vecNZcat().numel(),19)
      self.assertEqual(p.vecNZcat().size(),19)
      self.assertEqual(p.getNumel(),21)
      self.assertEqual(p.getSize(),19)
      
      self.assertTrue(p.lookup(('x',)) is p.x)
      self.assertTrue(p.lookup(('x',(0,))).toScalar().isEqual(p.x[0].toScalar()))
      self.assertTrue(p.lookup(('x',(1,))).toScalar().isEqual(p.x[1].toScalar()))
      self.assertTrue(p.lookup(('y',)) is p.y)
      for i in range(6):
        self.assertTrue(p.lookup(('y',(i,))).toScalar().isEqual(p.y[i].toScalar()))
      for i in range(3):
        for j in range(2):
          self.assertTrue(p.lookup(('y',(i,j))).toScalar().isEqual(p.y[i,j].toScalar()))
      self.assertTrue(p.lookup(('z',)) is p.z)
      self.assertTrue(p.lookup(('xother',)) is p.xother)
      self.assertTrue(p.lookup(('xother','a')) is p.xother.a)
      self.assertTrue(p.lookup(('xother','a',(0,))).isEqual(p.xother.a))
      self.assertTrue(p.lookup(('xother','b')) is p.xother.b)
      self.assertTrue(p.lookup(('xother','b',(1,1))).toScalar().isEqual(p.xother.b[1].toScalar()))
      
      
      for k in range(p.getSize()):
        roundtrip = p.lookup(p.reverselookup(k))
        if not(isinstance(roundtrip,SX)):
          roundtrip = roundtrip.toScalar()
        self.assertTrue(roundtrip.isEqual(p.vecNZcat()[k].toScalar()))
      
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
      b = []
      b.append(ssym("b1",3))
      b.append(ssym("b2",3))
      p.b = b
      p.c = ssym("c")
      p.freeze()
      

      self.assertEqual(p.veccat().numel(),9)
      self.assertEqual(p.veccat().size(),9)
      self.assertEqual(p.vecNZcat().numel(),9)
      self.assertEqual(p.vecNZcat().size(),9)
      self.assertEqual(p.getNumel(),9)
      self.assertEqual(p.getSize(),9)
      
      self.checkarray(array(p.i_a),DMatrix([[0],[1]]),"index")
      self.checkarray(array(p.i_b[0]),DMatrix([[2],[3],[4]]),"index")
      self.checkarray(array(p.i_b[1]),DMatrix([[5],[6],[7]]),"index")
      self.checkarray(array(p.i_c),DMatrix(8),"index")

      self.assertEqual(p.o_a,0,"Offset")
      self.assertEqual(p.o_b[0],2,"Offset")
      self.assertEqual(p.o_b[1],5,"Offset")
      self.assertEqual(p.o_c,8,"Offset")
      
      self.assertEqual(p.I_a,0,"Index")
      self.assertEqual(p.I_b[0],1,"Index")
      self.assertEqual(p.I_b[1],2,"Index")
      self.assertEqual(p.I_c,3,"Index")
      
      
      self.assertTrue(p.lookup(('b',)) is p.b)
      self.assertTrue(p.lookup(('b',0)) is p.b[0])
      self.assertTrue(p.lookup(('b',1)) is p.b[1])
      
      for i in range(3):
        self.assertTrue(p.lookup(('b',1,(i,))).toScalar().isEqual(p.b[1][i].toScalar()))

            
      for k in range(p.getSize()):
        roundtrip = p.lookup(p.reverselookup(k))
        if not(isinstance(roundtrip,SX)):
          roundtrip = roundtrip.toScalar()
        self.assertTrue(roundtrip.isEqual(p.vecNZcat()[k].toScalar()))
        
      p.b_[1].setAll(4)
      
      A = p.vecNZcat_()
 
      self.checkarray(A,DMatrix([0,0,0,0,0,4,4,4,0]),"vecNZcat")
     
      p.b_[0].setAll(3)
      A = p.veccat_()

      self.checkarray(A,DMatrix([0,0,3,3,3,4,4,4,0]),"vecNZcat")
      
      
      p = Variables()
      p.a = ssym("a",2)
      b = []
      b.append(ssym("b1",3,2))
      b.append(ssym("b2",3))
      p.b = b
      p.c = ssym("c")
      p.freeze()
      

      
      self.checkarray(array(p.i_a),DMatrix([[0],[1]]),"index")
      self.checkarray(array(p.i_b[0]),DMatrix([[2,5],[3,6],[4,7]]),"index")
      self.checkarray(array(p.i_b[1]),DMatrix([[8],[9],[10]]),"index")
      self.checkarray(array(p.i_c),DMatrix(11),"index")

      self.assertEqual(p.o_a,0,"Offset")
      self.assertEqual(p.o_b[0],2,"Offset")
      self.assertEqual(p.o_b[1],8,"Offset")
      self.assertEqual(p.o_c,11,"Offset")
      
      self.assertEqual(p.I_a,0,"Index")
      self.assertEqual(p.I_b[0],1,"Index")
      self.assertEqual(p.I_b[1],2,"Index")
      self.assertEqual(p.I_c,3,"Index")
      
      p.b_[1].setAll(4)
      
      A = p.vecNZcat_()
 
      self.checkarray(A,DMatrix([0,0,0,0,0,0,0,0,4,4,4,0]),"vecNZcat")
     
      p.b_[0].setAll(3)
      p.b_[0][2,0] = 9
      A = p.veccat_()

      self.checkarray(A,DMatrix([0,0,3,3,9,3,3,3,4,4,4,0]),"vecNZcat")

      p = Variables()
      p.a = msym("a",2)
      b = []
      b.append(msym("b1",3))
      b.append(msym("b2",3))
      p.b = b
      p.c = msym("c")
      p.freeze()
      
      f = MXFunction([p.veccat()],[p.a,p.b[0],p.b[1],p.c])
      f.init()
      
      f.input()[p.i_a]=[4,5]
      
      f.evaluate()
      
      p = Variables()
      p.a = ssym("a",2)
      b = []
      b.append(ssym("b1",3))
      b.append(ssym("b2",3))
      b.append([ssym("b3",3),ssym("b4",3)])
      p.b = b
      p.c = ssym("c")
      p.freeze()
      
      self.checkarray(array(p.i_a),DMatrix([[0],[1]]),"index")
      self.checkarray(array(p.i_b[0]),DMatrix([[2],[3],[4]]),"index")
      self.checkarray(array(p.i_b[1]),DMatrix([[5],[6],[7]]),"index")
      self.checkarray(array(p.i_b[2][0]),DMatrix([[8],[9],[10]]),"index")
      self.checkarray(array(p.i_b[2][1]),DMatrix([[11],[12],[13]]),"index")
      self.checkarray(array(p.i_c),DMatrix(14),"index")

      self.assertEqual(p.o_a,0,"Offset")
      self.assertEqual(p.o_b[0],2,"Offset")
      self.assertEqual(p.o_b[1],5,"Offset")
      self.assertEqual(p.o_b[2][0],8,"Offset")
      self.assertEqual(p.o_b[2][1],11,"Offset")
      self.assertEqual(p.o_c,14,"Offset")
      
      self.assertEqual(p.I_a,0,"Index")
      self.assertEqual(p.I_b[0],1,"Index")
      self.assertEqual(p.I_b[1],2,"Index")
      self.assertEqual(p.I_b[2][0],3,"Index")
      self.assertEqual(p.I_b[2][1],4,"Index")
      self.assertEqual(p.I_c,5,"Index")

      for k in range(p.getSize()):
        roundtrip = p.lookup(p.reverselookup(k))
        if not(isinstance(roundtrip,SX)):
          roundtrip = roundtrip.toScalar()
        self.assertTrue(roundtrip.isEqual(p.vecNZcat()[k].toScalar()))
        
        
      x = msym("a",2,3)
      
      p = Variables()
      p.a = x**2.2 # Was x**2 which simplifies to the _unary_ operation sq(x)
      p.b = [ x**3, sqrt(x)]
      p.c = [1/x]
      p.freeze(False)
      
      self.assertTrue(p.a.isBinary())
      
  def test_Numbers(self):
      p = Variables()

      p.x = ssym("x",2)
      p.z = ssym("z",2,4)
      p.y = ssym("y",3,2)

      xother = Variables()
      xother.a = SX("x")
      xother.b = diag(ssym("[a,b]"))
      xother.freeze()
      
      p.xother = xother
      
      p.freeze()
      
      p_ = Numbers(p)
      
      p_.y.setAll(3)
      
      p_.x = [4,5]
      
      p_.xother.a=12

      self.checkarray(p_.vecNZcat(),DMatrix([4,5,12,0,0,3,3,3,3,3,3,0,0,0,0,0,0,0,0]),"vecNZcat")
      
  def test_structure(self):
  
  
    s = struct(['x','y','z'])
    
    print s
  
    with self.assertRaises(Exception):
      struct_ssym(['x','x','z'])
    with self.assertRaises(Exception): 
      struct_ssym(['x','y','z','y'])
      

    
    s = struct_ssym([('x','y'),'z'])
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
    
    s = struct_ssym(['x','y','z'],order=['x','y','z'])
    with self.assertRaises(Exception): 
      s = struct_ssym([('x','y'),'z'],order=['x','y','z'])
    
    s = struct_ssym(['x','y','z'])
    self.assertEqual(s.size,3)
    self.assertEqual(s[vertcat,['x','y']].shape[0],2)
    self.assertTrue(isinstance(s[list,vertcat,['x','y']],list))
    
    s = struct_ssym([entry('x',repeat=5),entry('y',repeat=6),entry('z')])
    self.assertEqual(s.size,12)
    self.assertEqual(len(s["x"]),5)
    self.assertEqual(len(s["y"]),6)
    self.assertTrue(s.cat.at(1).getName().startswith("x"))
    s = struct_ssym([(entry('x',repeat=5),entry('y',repeat=[6,5])),entry('z')])
    
        
    self.checkarray(s.f["x"],[0,6,12,18,24])
    self.checkarray(s.f["z"],[35])
    self.checkarray(s.f["y"],[1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,19,20,21,22,23,25,26,27,28,29,30,31,32,33,34])
    
    self.assertEqual(s.size,36)
    self.assertEqual(len(s["x"]),5)
    self.assertEqual(len(s["y"]),6)
    self.assertTrue(s.cat.at(1).getName().startswith("y"))
    s = struct_ssym([entry("x",shape=(3,2)),entry("y",shape=2),entry("z",shape=sp_dense(3)),entry("w",shape=sp_tril(5))])
    self.assertEqual(s.size,6+2+3+15)
    self.assertTrue(s["x"].sparsity()==sp_dense(3,2))
    self.assertTrue(s["y"].sparsity()==sp_dense(2,1))
    self.assertTrue(s["z"].sparsity()==sp_dense(3,1))
    self.assertTrue(s["w"].sparsity()==sp_tril(5))
    
    x  = ssym("x",2)
    x2 = ssym("x2",2)
    s = struct_ssym([entry('a',sym=x),'y','z'])
    self.assertTrue(s.cat.at(0).getName().startswith("x"))
    self.assertEqual(s.size,4)
    with self.assertRaises(Exception):
      struct_ssym([entry('a',sym=x+x),'y','z'])
    with self.assertRaises(Exception):
      struct_ssym([entry('a',sym=[x+x]),'y','z'])
    s = struct_ssym([entry('a',sym=vertcat([x,x2])),'y','z'])
    self.assertEqual(s.size,6)
    with self.assertRaises(Exception):
      s = struct_ssym([entry('a',repeat=6,sym=x),'y','z'])
    
    s = struct_ssym(['x','y','z'])
 
    S = struct_ssym([entry("X",sym=s)])
    self.assertTrue(S.cat.at(0).getName()=="x")
    S = struct_ssym([entry("X",struct=s)])
    self.assertTrue(S.cat.at(0).getName()=="X_x")
    S = struct_ssym([entry("X",repeat=[5],struct=s)])
    self.assertEqual(S.size,15)
    
    s0 = struct_ssym(['x','y','z'])
    s1 = struct_ssym(['x','y'])
    S0 = struct_ssym([entry("X",struct=s0),entry("Y",repeat=[5],struct=s0),entry("Z",struct=s0)])
    S1 = struct_ssym([entry("X",struct=s1),entry("Y",repeat=[5],struct=s1)])
    num0 = S0(0)
    num1 = S1(1)

    num0[nesteddict] = num1[nesteddict]
    
    S = struct_ssym([entry("P",shapestruct=(s0,s0))])
 
    num = S(0)
    num["P"] = DMatrix([[1,2,3],[4,5,6],[7,8,9]])
    self.checkarray(num["P",index["x"],index["y"]],DMatrix([2]))
    self.checkarray(num["P",index["x"],:],DMatrix([1,2,3]).T)
    self.checkarray(num["P",:,index["y"]],DMatrix([2,5,8]))
    self.checkarray(num["P",:,:],DMatrix([[1,2,3],[4,5,6],[7,8,9]]))
    
    self.checkarray(num["P",indexf[["x","y"]],indexf[["z","x"]]],DMatrix([[3,1],[6,4]]))
    
    self.checkarray(num["P",index[list,vertcat,["x","y"]],index[list,vertcat,["z","x"]]],DMatrix([[3,1],[6,4]]))
    self.checkarray(num["P",index[vertcat,["x","y"]],index[vertcat,["z","x"]]],DMatrix([3,4]))
    
    S = struct_ssym([entry("P",shapestruct=s0)])
    
    num = S(0)
    num["P"] = DMatrix([1,2,3])
    self.checkarray(num["P",index["x"]],DMatrix([1]))
    self.checkarray(num["P",:],DMatrix([1,2,3]))
    
    self.checkarray(num["P",index["x"],0],DMatrix([1]))
    self.checkarray(num["P",:,0],DMatrix([1,2,3]))
    
    with self.assertRaises(Exception):
      num["P",:,index["x"]]
    
    
    S = struct_ssym([entry("P",shapestruct=(3,s0))])
 
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
    
    S = struct_ssym([entry("P",shapestruct=s0)])
 
    
    
    s0 = struct_ssym(['x','y',entry('q',shape=4),'z'])

    S = struct_ssym([entry("P",shapestruct=(s0,s0))])
    
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
      struct_msym(['x','x','z'])
    with self.assertRaises(Exception): 
      struct_msym(['x','y','z','y'])
    
    s = struct_msym([('x','y'),'z'])
    s = struct_msym(['x','y','z'],order=['x','y','z'])
    with self.assertRaises(Exception): 
      s = struct_msym([('x','y'),'z'],order=['x','y','z'])
    
    s = struct_msym(['x','y','z'])
    self.assertEqual(s.size,3)
    s = struct_msym([entry('x',repeat=5),entry('y',repeat=6),entry('z')])
    self.assertEqual(s.size,12)
    self.assertEqual(len(s["x"]),5)
    self.assertEqual(len(s["y"]),6)
    s = struct_msym([(entry('x',repeat=5),entry('y',repeat=[6,5])),entry('z')])
    self.assertEqual(s.size,36)
    self.assertEqual(len(s["x"]),5)
    self.assertEqual(len(s["y"]),6)
   
    
    s = struct_msym([entry("x",shape=(3,2)),entry("y",shape=2),entry("z",shape=sp_dense(3)),entry("w",shape=sp_tril(5))])
    self.assertEqual(s.size,6+2+3+15)
    self.assertTrue(s["x"].sparsity()==sp_dense(3,2))
    self.assertTrue(s["y"].sparsity()==sp_dense(2,1))
    self.assertTrue(s["z"].sparsity()==sp_dense(3,1))
    self.assertTrue(s["w"].sparsity()==sp_tril(5))
    
    x  = msym("x",2)
    x2 = msym("x2",2)
    with self.assertRaises(Exception):
      s = struct_msym([entry('x',sym=x),'y','z'])
    with self.assertRaises(Exception):
      struct_msym([entry('x',sym=x+x),'y','z'])
    with self.assertRaises(Exception):
      struct_msym([entry('x',sym=[x+x]),'y','z'])
    with self.assertRaises(Exception):
      struct_msym([entry('x',sym=vertcat([x,x2])),'y','z'])
    with self.assertRaises(Exception):
      s = struct_,msym([(2,'x',[x,x2]),'y','z'])
    
    s = struct_msym(['x','y','z'])
    S = struct_msym([entry("X",struct=s)])
    S = struct_msym([entry("X",repeat=[5],struct=s)])
    
    x = ssym("x",2)
    y0 = sin(x) 
    y1 = cos(x)
    
    v = struct_SX([
      entry("x",expr=x),
      entry("y",expr=[y0,y1]),
      entry("w",expr=[[x,x**2],[x**3,x**4]])
    ])
    
    
    x = msym("x",2)
    m = struct_msym(['a','b'])
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
    self.assertTrue(isEqual(V["x"],x))
    self.assertTrue(isEqual(V["y",0],y0))
    self.assertTrue(isEqual(V["y",1],y1))
    self.assertEqual(V["y",0,'a'].shape,(1,1))
    
    with self.assertRaises(Exception):
      V["y",0] = msym("x",4) # shape mismatch
    abc = msym("abc",2)
    V["y",0] = abc
    self.assertTrue(isEqual(V["y",0],abc))

    states = struct_ssym([
                entry('x'),
                entry('y'),
                entry('z'),
                entry('u',shape=sp_dense(4)),
                entry('v',repeat=[4,2]),
                entry('w',repeat=[6]),
                entry('p',repeat=[9],shape=sp_dense(6))
             ],order=['x','y','z','u',('v','w'),'p'])
             
    shooting = struct_ssym([entry('X',struct=states,repeat=[4,5]),entry('U',repeat=[3])],order=[('X','U')])


    
    self.assertEqual(shooting.size,1503)
    s = shooting["X",:,:,{}]
    self.assertTrue(isinstance(s,list))
    self.assertEqual(len(s),4)
    self.assertTrue(isinstance(s[0],list))
    self.assertEqual(len(s[0]),5)
    self.assertTrue(isinstance(s[0][0],dict))
    self.assertTrue('x' in s[0][0])
    self.assertEqual(len(s[0][0]),7)
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
    
    V = struct_msym(shooting)
    
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

    S = struct_ssym([entry("X",repeat=12,struct=s)])

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
    d = DMatrix.zeros(12,s.size)
    a = s.repeated(d)
    
    a[:,"x"] = range(12)
    
    self.checkarray(a[4,"x"],DMatrix([4]))
    self.checkarray(d,horzcat([range(12),DMatrix.zeros(12),DMatrix.zeros(12)]))

  def test_structure_squared_dmatrix(self):
    self.message("squared dmatrix")
    s = struct(["x","y","z"])
    d = DMatrix.zeros(3,3)
    a = s.squared(d)
    
    a["x","y"] = 2
    a["y","x"] = 1
    
    self.checkarray(d,DMatrix([[0,2,0],[1,0,0],[0,0,0]]))

  def test_structure_squared_repeated_dmatrix(self):
    self.message("squared repeated dmatrix")
    s = struct(["x","y","z"])
    d = DMatrix.zeros(9,3)
    a = s.squared_repeated(d)
    
    a[0,"x","y"] = 2
    a[0,"y","x"] = 1
    a[-1,"x","y"] = 2
    a[-1,"y","x"] = 1
    self.checkarray(d,DMatrix([[0,2,0],[1,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,2,0],[1,0,0],[0,0,0]]))
    
  def test_typemaps(self):
    self.message("typemaps")
    s = struct(["x","y","z"])
    d = DMatrix.zeros(3,3)
    a = s.squared(d)
    
    print sin(a)
    print a+1
    
  def test_sparse(self):
    a = struct_ssym([entry("a",shape=sp_diag(5))])
    b = struct_msym([(entry("b",struct=a))])

    self.checkarray(b["b"].shape,(5,1))
    self.checkarray(b["b","a"].shape,(5,5))
    
if __name__ == '__main__':
    unittest.main()

