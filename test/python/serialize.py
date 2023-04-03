#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
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
#LH9JCUD52HXU
from casadi import *
import casadi as c
import numpy
import unittest
from types import *
from helpers import *
import random
from collections import defaultdict


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
          if "CommonExternal" in str(e):
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
            self.checkarray(o,o_ref,digits=digits)
  
      
if __name__ == '__main__':
    unittest.main()
