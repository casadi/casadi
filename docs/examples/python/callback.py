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
# -*- coding: utf-8 -*-
from casadi import *
#
# How to use Callback
# Joel Andersson
#

class MyCallback(Callback):
  def __init__(self, name, d, opts={}):
    Callback.__init__(self)
    self.d = d
    self.construct(name, opts)

  # Number of inputs and outputs
  def get_n_in(self): return 1
  def get_n_out(self): return 1

  # Initialize the object
  def init(self):
     print('initializing object')

  # Evaluate numerically
  def eval(self, arg):
    x = arg[0]
    f = sin(self.d*x)
    return [f]

# Use the function
f = MyCallback('f', 0.5)
res = f(2)
print(res)

# You may call the Callback symbolically
x = MX.sym("x")
print(f(x))

# Derivates OPTION 1: finite-differences
eps = 1e-5
print((f(2+eps)-f(2))/eps)

f = MyCallback('f', 0.5, {"enable_fd":True})
J = Function('J',[x],[jacobian(f(x),x)])
print(J(2))

# Derivates OPTION 2: Supply forward mode
# Example from https://www.youtube.com/watch?v=mYOkLkS5yqc&t=4s

class Example4To3(Callback):
  def __init__(self, name, opts={}):
    Callback.__init__(self)
    self.construct(name, opts)

  def get_n_in(self): return 1
  def get_n_out(self): return 1

  def get_sparsity_in(self,i):
    return Sparsity.dense(4,1)

  def get_sparsity_out(self,i):
    return Sparsity.dense(3,1)

  # Evaluate numerically
  def eval(self, arg):
    a,b,c,d = vertsplit(arg[0])
    ret = vertcat(sin(c)*d+d**2,2*a+c,b**2+5*c)
    return [ret]


class Example4To3_Fwd(Example4To3):
  def has_forward(self,nfwd):
    # This example is written to work with a single forward seed vector
    # For efficiency, you may allow more seeds at once
    return nfwd==1
  def get_forward(self,nfwd,name,inames,onames,opts):
    
    class ForwardFun(Callback):
      def __init__(self, opts={}):
        Callback.__init__(self)
        self.construct(name, opts)

      def get_n_in(self): return 3
      def get_n_out(self): return 1

      def get_sparsity_in(self,i):
        if i==0: # nominal input
          return Sparsity.dense(4,1)
        elif i==1: # nominal output
          return Sparsity(3,1)
        else: # Forward seed
          return Sparsity.dense(4,1)

      def get_sparsity_out(self,i):
        # Forward sensitivity
        return Sparsity.dense(3,1)

      # Evaluate numerically
      def eval(self, arg):
        a,b,c,d = vertsplit(arg[0])
        a_dot,b_dot,c_dot,d_dot = vertsplit(arg[2])
        print("Forward sweep with", a_dot,b_dot,c_dot,d_dot)
        w0 = sin(c)
        w0_dot = cos(c)*c_dot
        w1 = w0*d
        w1_dot = w0_dot*d+w0*d_dot
        w2 = d**2
        w2_dot = 2*d_dot*d
        r0 = w1+w2
        r0_dot = w1_dot + w2_dot
        w3 = 2*a
        w3_dot = 2*a_dot
        r1 = w3+c
        r1_dot = w3_dot+c_dot
        w4 = b**2
        w4_dot = 2*b_dot*b
        w5 = 5*w0
        w5_dot = 5*w0_dot
        r2 = w4+w5
        r2_dot = w4_dot + w5_dot
        ret = vertcat(r0_dot,r1_dot,r2_dot)
        return [ret]
    # You are required to keep a reference alive to the returned Callback object
    self.fwd_callback = ForwardFun()
    return self.fwd_callback

    
f = Example4To3_Fwd('f')
x = MX.sym("x",4)
J = Function('J',[x],[jacobian(f(x),x)])
print(J(vertcat(1,2,0,3)))

# Derivates OPTION 3: Supply reverse mode

class Example4To3_Rev(Example4To3):
  def has_reverse(self,nadj):
    # This example is written to work with a single forward seed vector
    # For efficiency, you may allow more seeds at once
    return nadj==1
  def get_reverse(self,nfwd,name,inames,onames,opts):
    
    class ReverseFun(Callback):
      def __init__(self, opts={}):
        Callback.__init__(self)
        self.construct(name, opts)

      def get_n_in(self): return 3
      def get_n_out(self): return 1

      def get_sparsity_in(self,i):
        if i==0: # nominal input
          return Sparsity.dense(4,1)
        elif i==1: # nominal output
          return Sparsity(3,1)
        else: # Reverse seed
          return Sparsity.dense(3,1)

      def get_sparsity_out(self,i):
        # Reverse sensitivity
        return Sparsity.dense(4,1)

      # Evaluate numerically
      def eval(self, arg):
        a,b,c,d = vertsplit(arg[0])
        r0_bar,r1_bar,r2_bar = vertsplit(arg[2])
        print("Reverse sweep with", r0_bar, r1_bar, r2_bar)
        w0 = sin(c)
        w1 = w0*d
        w2 = d**2
        r0 = w1+w2
        w3 = 2*a
        r1 = w3+c
        w4 = b**2
        w5 = 5*w0
        r2 = w4+w5
        w4_bar = r2_bar
        w5_bar = r2_bar
        w0_bar = 5*w5_bar
        b_bar = 2*b*w4_bar
        w3_bar = r1_bar
        c_bar = r1_bar
        a_bar = 2*w3_bar
        w1_bar = r0_bar
        w2_bar = r0_bar
        d_bar = 2*d*w2_bar
        w0_bar = w0_bar + w1_bar*d
        d_bar = d_bar + w0*w1_bar
        c_bar = c_bar + cos(c)*w0_bar
        ret = vertcat(a_bar,b_bar,c_bar,d_bar)
        return [ret]
    # You are required to keep a reference alive to the returned Callback object
    self.rev_callback = ReverseFun()
    return self.rev_callback

    
f = Example4To3_Rev('f')
x = MX.sym("x",4)
J = Function('J',[x],[jacobian(f(x),x)])
print(J(vertcat(1,2,0,3)))

# Derivates OPTION 4: Supply full Jacobian

class Example4To3_Jac(Example4To3):
  def has_jacobian(self): return True
  def get_jacobian(self,name,inames,onames,opts):
    class JacFun(Callback):
      def __init__(self, opts={}):
        Callback.__init__(self)
        self.construct(name, opts)

      def get_n_in(self): return 2
      def get_n_out(self): return 1

      def get_sparsity_in(self,i):
        if i==0: # nominal input
          return Sparsity.dense(4,1)
        elif i==1: # nominal output
          return Sparsity(3,1)

      def get_sparsity_out(self,i):
        return sparsify(DM([[0,0,1,1],[1,0,1,0],[0,1,1,0]])).sparsity()

      # Evaluate numerically
      def eval(self, arg):
        a,b,c,d = vertsplit(arg[0])
        ret = DM(3,4)
        ret[0,2] = d*cos(c)
        ret[0,3] = sin(c)+2*d
        ret[1,0] = 2
        ret[1,2] = 1
        ret[2,1] = 2*b
        ret[2,2] = 5
        return [ret]

    # You are required to keep a reference alive to the returned Callback object
    self.jac_callback = JacFun()
    return self.jac_callback

f = Example4To3_Jac('f')
x = MX.sym("x",4)
J = Function('J',[x],[jacobian(f(x),x)])
print(J(vertcat(1,2,0,3)))

    
