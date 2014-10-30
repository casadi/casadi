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
#! CasADi tutorial
#! ==================
#! This tutorial file explains the use of structures
#! Structures are a python-only feature from the tools library:

from casadi.tools import *
#! The struct tools offer a way to structure your symbols and data
#! It allows you to make abstraction of ordering and indices.
#! 
#! Put simply, the goal is to eleminate code with 'magic numbers' numbers such as:
#!
#!    f = MXFunction([V],[ V[214] ])  # time
#!    ...
#!    x_opt = solver.getOutput()[::5]    # Obtain all optimized x's
#!
#!  and replace it with
#!    f = MXFunction([V],[ V["T"] ])
#!    ...
#!    shooting(solver.getOutput())["x",:]
#!
#! Introduction
#! --------------------

#! Create a structured SX.sym
states = struct_symSX(["x","y","z"])

print states

#! Superficially, states behaves like a dictionary
print states["y"]

#! To obtain aliases, use the Ellipsis index:
x,y,z = states[...]

#! The cat attribute will return the concatenated version of the struct. This will always be a column vector
print states.cat

#! This structure is of size:
print states.size, "=", states.cat.shape
  
  
f = SXFunction([states.cat],[x*y*z])

#! In many cases, states will be auto-cast to SX:
f = SXFunction([states],[x*y*z])

#! Expanded structure syntax and ordering 
#! --------------------------------------

#! The structure defnition above can also be written in expanded syntax:
simplestates = struct_symSX([
    entry("x"),
    entry("y"),
    entry("z")
  ])

#! More information can be attached to the entries
#!   shape argument  : specify sparsity/shape
states = struct_symSX([
    entry("x",shape=3),
    entry("y",shape=(2,2)),
    entry("z",shape=Sparsity.tril(2))
  ])
  
print states["x"]
print states["y"]
print states["z"]

#! Note that the cat version of this structure does only contain the nonzeros
print states.cat

#!   repeat argument  : specify nested lists
states = struct_symSX([
    entry("w",repeat=2),
    entry("v",repeat=[2,3]),
  ])

print states["w"]
print states["v"]

#! Notice that all v variables come before the v entries:
for i,s in enumerate(states.labels()):
  print i, s
  
#! We can influency this order by introducing a grouping bracket:

states = struct_symSX([
    "a",
    ( entry("w",repeat=2),
      entry("v",repeat=[2,3])
    ),
    "b"
  ])
  
#! Notice how the w and v variables are now interleaved:
for i,s in enumerate(states.labels()):
  print i, s
  
#! Nesting, Values and PowerIndex
#! -------------------------------

#! Structures can be nested. For example consider a statespace of two cartesian coordinates and a quaternion
states = struct_symSX(["x","y",entry("q",shape=4)])

shooting = struct_symSX([
  entry("X",repeat=[5,3],struct=states),
  entry("U",repeat=4,shape=1),
])

print shooting.size

#! The canonicalIndex is the combination of strings and numbers that uniquely defines the entries of a structure:
print shooting["X",0,0,"x"]

#! If we use more exoctic indices, we call this a powerIndex
print shooting["X",:,0,"x"]

#! Having structured symbolics is one thing. The numeric structures can be derived:
#! The following line allocates a DMatrix of correct size, initialised with zeros
init = shooting(0)

print init.cat

# We can use the powerIndex in the context of indexed assignent, too:
init["X",0,-1,"y"] = 12

#! The corresponding numerical value has changed now:
print init.cat

#! The entry that changed is in fact this one:
print init.f["X",0,-1,"y"]
print init.cat[13]

#! One can lookup the meaning of the 13th entry in the cat version as such:
#! Note that the canonicalIndex does not contain negative numbers
print shooting.getCanonicalIndex(13)

print shooting.labels()[13]

#! Other datatypes
#! ----------------

#! A symbolic structure is immutable

try:
  states["x"] = states["x"]**2
except Exception as e:
  print "Oops:", e
  
#! If you want to have a mutable variant, for example to contian the right hand side of an ode, use struct_SX:
rhs = struct_SX(states)

rhs["x"] = states["x"]**2
rhs["y"] = states["y"]*states["x"]
rhs["q"] = -states["q"]

print rhs.cat

#! Alternatively, you can supply the expressions at defintion time:
x,y,q = states[...]
rhs = struct_SX([
    entry("x",expr=x**2),
    entry("y",expr=x*y),
    entry("q",expr=-q)
  ])
  
print rhs.cat


#! One can also construct symbolic MX structures
V = struct_symMX(shooting)

print V

#! The catted version is one single MX from which all entries are derived:
print V.cat
print V.shape
print V["X",0,-1,"y"]

#! Similar to struct_SX, we have struct_MX:
V = struct_MX([
    (
    entry("X",expr=[[ MX.sym("x",6)**2 for j in range(3)] for i in range(5)]),
    entry("U",expr=[ -MX.sym("u") for i in range(4)])
    )
  ])

#! By default SX.sym structure constructor will create new SX.syms.
#! To recycle one that is already available, use the 'sym' argument: 
qsym = SX.sym("quaternion",4)
states = struct_symSX(["x","y",entry("q",sym=qsym)])
print states.cat

#! The 'sym' feature is not available for struct_MX, since it will construct one parent MX.

#! More powerIndex
#! ----------------

#! As illustrated before, powerIndex allows slicing
print init["X",:,:,"x"]

#! The repeated method duplicates its argument a number of times such that it matches the length that is needed at the lhs
init["X",:,:,"x"] = repeated(range(3))

print init["X",:,:,"x"]

#! Callables/functions can be thrown in in the powerIndex at any location.
#! They operate on subresults obtain from resolving the remainder of the powerIndex

print init["X",:,horzcat,:,"x"]
print init["X",vertcat,:,horzcat,:,"x"]
print init["X",blockcat,:,:,"x"]

#! Set all quaternions to 1,0,0,0
init["X",:,:,"q"] = repeated(repeated(DMatrix([1,0,0,0])))

#! {} can be used in the powerIndex to expand into a dictionary once
init["X",:,0,{}] = repeated({"y": 9})

print init["X",:,0,{}]

#! lists can be used in powerIndex in both list context or dict context:
print shooting["X",[0,1],[0,1],"x"]
print shooting["X",[0,1],0,["x","y"]]

#! nesteddict can be used to expand into a dictionary recursively
print init[nesteddict]

#! ... will expand entries as an ordered list
print init["X",:,0,...]

#! If the powerIndex ends at the boundary of a structure, it's catted version is returned:
print init["X",0,0]

#! If the powerIndex is longer than what could be resolved as structure, the remainder, extraIndex, is passed onto the resulting Casadi-matrix-type
print init["X",blockcat,:,:,"q",0]

print init["X",blockcat,:,:,"q",0,0]

#! shapeStruct and delegated indexing
#! -----------------------------------

#! When working with covariance matrices, both the rows and columns relate to states

states = struct(["x","y",entry("q",repeat=2)])
V = struct_symSX([
      entry("X",repeat=5,struct=states),
      entry("P",repeat=5,shapestruct=(states,states))
    ])

#! P has a 4x4 shape
print V["P",0]

#! Now we can use powerIndex-style in the extraIndex:
print V["P",0,["x","y"],["x","y"]]

#! There is a problem when we wich to use the full potential of powerIndex in these extraIndices:
#! The following is in fact invalid python syntax:
#   V["P",0,["q",:],["q",:]] 

#! We resolve this by using delegater objects index/indexf:

print V["P",0,indexf["q",:],indexf["q",:]] 

#! Of course, in this basic example, also the following would be allowed
print V["P",0,"q","q"] 

#! Prefixing
#! -----------

#! The prefix attribute allows you to create shorthands for long powerIndices

states = struct(["x","y","z"])
V = struct_symSX([
      entry("X",repeat=[4,5],struct=states)
])

num = V()

#! Consider the following statements:

num["X",0,0,"x"] = 1
num["X",0,0,"y"] = 2
num["X",0,0,"z"] = 3

#! Note the common part ["X",0,0].
#! We can pull this apart with prefix:

initial = num.prefix["X",0,0]

initial["x"] = 1
initial["y"] = 2
initial["z"] = 3

#! This is equivalent to the longer statements above

#! Helper constructors
#! -------------------

#! If you work with Simulator, ControlSimulator, you typically end up
#! with wanting to index a DMatrix that is n x N
#! with n the size of a statespace and N an arbitrary integer

states = struct(["x","y","z"])

#! We artificially construct here a DMatrix that could be a Simulator output.
output = DMatrix.zeros(states.size,8)

#! The helper construct is 'repeated' here. Instead of "states(output)", we have
outputs = states.repeated(output)

#! Now we have an object that supports powerIndexing:
outputs[-1] = DMatrix([1,2,3])
outputs[:,"x"] = range(8) 

print output
print outputs[5,{}]

#! Next we represent the 'squared' helper construct
#! Imagine we somehow obtain a matrix that represents covariance
P0 = DMatrix.zeros(states.size,states.size)

#! We can conveniently access it as follows:
P = states.squared(P0)
P["x","y"] = 2
P["y","x"] = 3

print P0

#! P itself is a rather queer object
print P

#! You can access its concents with a call:
print P()

#! But often, it will behave like a DMatrix transparantly:
P0.set(P)

#! Next we represent the 'squared_repeated' helper construct
#! Imagine we somehow obtain a matrix that represents a horizontal concatenation of covariance
P0 = horzcat([DMatrix.zeros(states.size,states.size),DMatrix.ones(states.size,states.size)])

#! We can conveniently access it as follows:
P = states.squared_repeated(P0)
P[0,"x","y"] = 2
P[:,"y","x"] = 3

print P0

#! Finally, we present the 'product' helper construct
controls = struct(["u","v"])

J0 = DMatrix.zeros(states.size,controls.size)

J = states.product(controls,J0)

J[:,"u"] = 3
J[["x","z"],:] = 2

print J()

#! Saving and loading
#! -------------------

#! It is possible to save and load some types of structures.
#! Supported types are pure structures (the ones created with 'struct') and numeric structures.
#! 
#!  Saving:
#!     mystructure.save("myfilename")
#!  Loading:
#!     struct_load("myfilename")

