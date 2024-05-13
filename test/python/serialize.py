#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
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
#LH9JCUD52HXU
from casadi import *
import casadi as c
import numpy
import unittest
from types import *
from helpers import *
import random
from collections import defaultdict
import sys


class SerializeTests(casadiTestCase):


  def test_compat(self):
  
    dir = "serialize_3.5.5"
    if not os.path.isdir(dir):
        return
    
    errors = {}
    
    inputs = defaultdict(list)
    outputs = defaultdict(list)
    names = os.listdir(dir)
    functions = [n for n in names if n.endswith(".casadi")]
    for n in names:
        data = inputs if "_in" in n else outputs
        data[n.split("_")[0]].append(n)
    for fun in functions:
      try:
          f = Function.load(os.path.join(dir,fun))
      except Exception as e:
          if "DllLibrary::init_handle" in str(e):
            continue
          if "DeSerialization of Integrator failed" in str(e):
            continue
          else:
             raise Exception(str(e))

      key = fun.split(".")[0]
      for in_name,out_name in zip(sorted(inputs[key]),sorted(outputs[key])):
          inp = f.generate_in(os.path.join(dir,in_name))
          outp_ref = f.generate_out(os.path.join(dir,out_name))
          outp = f.call(inp)
          self.assertEqual(len(outp),len(outp_ref))
          for o,o_ref in zip(outp,outp_ref):
            digits = 15
            if "SuperscsInterface" in str(f):
                digits = 7
            if sys.platform=="darwin":
                digits = 7 # Bug?
            self.checkarray(o,o_ref,digits=digits,message=fun+"/"+str(f))
            
  def test_identity(self):
    obj = ["foo",{"foo":"bar"},[{"a":3},{"b":9}],["a",5],{"foo": ["a",5]},{"foo": [["a",5],["b",2]]},[["a",5],["b",2]],[[["a",5],["b",2]]]]
    
    def check_equal(a,b):
        if isinstance(a,dict):
            assert list(sorted(a.keys()))==list(sorted(b.keys()))
            for k in a.keys():
                check_equal(a[k],b[k])
        elif isinstance(a,list):
            assert(len(a)==len(b))
            for i in range(len(a)):
                check_equal(a[i],b[i])
        else:
            assert a==b
        
    for e in obj:
        ss = StringSerializer()
        ss.pack(e)
        
        ds = StringDeserializer(ss.encode())
        r = ds.unpack()
        
        print(e,r)
        check_equal(e,r)
      
if __name__ == '__main__':
    unittest.main()
