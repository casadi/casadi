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
#
from casadi import *
import casadi as c
import numpy
from numpy import random, array
import unittest
from types import *
from helpers import *
import itertools


class DocStrings(casadiTestCase):

  def test_1(self):
    doc = Function.jacobian.__doc__
    print(doc)
    self.assertTrue("Calculate all Jacobian blocks Generates a function that take" not in doc)

  def test_2(self):
    doc = cse.__doc__
    self.assertTrue("Common subexpression elimination" in doc)

  def test_3(self):
    doc = nlpsol.__doc__
    self.assertTrue("Create an NLP solver" in doc)
    self.assertTrue("WORHP" in doc)
    
  def test_4(self):
    doc = conic.__doc__
    self.assertTrue("Second-order cone constraints" in doc)
    self.assertTrue("olves QPs using a Mehrotra predictor-corrector interior poi" in doc)
    self.assertTrue("min_lam" in doc)
    self.assertTrue("dump_out" in doc)
    self.assertTrue("CONIC_A" in doc)
    
  def test_5(self):
    doc = collocation_interpolators.__doc__
    self.assertTrue("A collocation method" in doc)

  def test_6(self):
    doc = adj.__doc__
    self.assertTrue("Matrix adjoint" in doc)

  def test_7(self):
    doc = taylor.__doc__
    self.assertTrue("sin(x)" in doc)
    
  def test_8(self):
    doc = erfinv.__doc__
    self.assertTrue("Inverse error function" in doc)


if __name__ == '__main__':
    unittest.main()
