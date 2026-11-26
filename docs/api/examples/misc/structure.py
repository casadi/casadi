# %%
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
# %% [markdown]
# # CasADi tutorial: Structures
#
# This tutorial file explains the use of **structures**.
# Structures are a Python-only feature from the `casadi.tools` library:

# %%
from casadi import DM, MX, SX, Function, Sparsity, blockcat, horzcat, vertcat
from casadi.tools.structure import (
    entry,
    indexf,
    nesteddict,
    repeated,
    struct,
    struct_MX,
    struct_SX,
    struct_symMX,
    struct_symSX,
)

# %% [markdown]
# The struct tools offer a way to structure your symbols and data.
# It allows you to make abstraction of ordering and indices.
#
# Put simply, the goal is to eliminate code with 'magic numbers' such as:
#
# ```python
# f = Function('f', [V], [V[214]])          # time
# ...
# x_opt = solver.getOutput()[::5]           # Obtain all optimized x's
# ```
#
# and replace it with:
#
# ```python
# f = Function('f', [V], [V["T"]])
# ...
# shooting(solver.getOutput())["x", :]
# ```

# %% [markdown]
# ## Introduction
#
# Create a structured SX.sym
# %%
states = struct_symSX(["x", "y", "z"])
print(states)

# %% [markdown]
# Superficially, `states` behaves like a dictionary:

# %%
print(states["y"])

# %% [markdown]
# To obtain aliases, use the Ellipsis index:

# %%
x, y, z = states[...]

# %% [markdown]
# The `cat` attribute will return the concatenated version of the struct.
# This will always be a column vector:

# %%
print(states.cat)

# %% [markdown]
# This structure is of size:

# %%
print(states.shape, "=", states.cat.shape)

# %%
f = Function("f", [states.cat], [x * y * z])

# %% [markdown]
# In many cases, `states` will be auto-cast to SX:

# %%
f = Function("f", [states], [x * y * z])

# %% [markdown]
# ## Expanded structure syntax and ordering
#
# The structure definition above can also be written in expanded syntax:

# %%
simplestates = struct_symSX([entry("x"), entry("y"), entry("z")])

# %% [markdown]
# More information can be attached to the entries:
#
# * `shape` argument: specify sparsity/shape

# %%
states = struct_symSX(
    [entry("x", shape=3), entry("y", shape=(2, 2)), entry("z", shape=Sparsity.lower(2))]
)

print(states["x"])
print(states["y"])
print(states["z"])

# %% [markdown]
# Note that the `cat` version of this structure only contains the nonzeros:

# %%
print(states.cat)

# %% [markdown]
# * `repeat` argument: specify nested lists

# %%
states = struct_symSX([entry("w", repeat=2), entry("v", repeat=[2, 3])])
print(states["w"])
print(states["v"])

# %% [markdown]
# Notice that all `w` variables come before the `v` entries:

# %%
for i, s in enumerate(states.labels()):
    print(i, s)

# %% [markdown]
# We can influence this order by introducing a grouping bracket:

# %%
states = struct_symSX(["a", (entry("w", repeat=2), entry("v", repeat=[2, 3])), "b"])

# %% [markdown]
# Notice how the `w` and `v` variables are now interleaved:

# %%
for i, s in enumerate(states.labels()):
    print(i, s)

# %% [markdown]
# ## Nesting, Values and PowerIndex
#
# Structures can be nested. For example, consider a statespace of two cartesian
# coordinates and a quaternion:

# %%
states = struct_symSX(["x", "y", entry("q", shape=4)])
shooting = struct_symSX(
    [
        entry("X", repeat=[5, 3], struct=states),
        entry("U", repeat=4, shape=1),
    ]
)
print(shooting.shape)

# %% [markdown]
# The *canonicalIndex* is the combination of strings and numbers that uniquely
# defines the entries of a structure:

# %%
print(shooting["X", 0, 0, "x"])

# %% [markdown]
# If we use more exotic indices, we call this a *powerIndex*:

# %%
print(shooting["X", :, 0, "x"])

# %% [markdown]
# Having structured symbolics is one thing. The numeric structures can be derived:
#
# The following line allocates a DM of correct size, initialised with zeros:

# %%
init = shooting(0)
print(init.cat)

# %% [markdown]
# We can use the powerIndex in the context of indexed assignment, too:

# %%
init["X", 0, -1, "y"] = 12

# %% [markdown]
# The corresponding numerical value has changed now:

# %%
print(init.cat)

# %% [markdown]
# The entry that changed is in fact this one:

# %%
print(init["X", 0, -1, "y"])
print(init.cat[13])

# %% [markdown]
# One can look up the meaning of the 13th entry in the `cat` version as such:
# Note that the canonicalIndex does not contain negative numbers:

# %%
print(shooting.getCanonicalIndex(13))
print(shooting.labels()[13])

# %% [markdown]
# ## Other datatypes
#
# A symbolic structure is immutable:

# %%
try:
    states["x"] = states["x"] ** 2
except Exception as e:
    print("Oops:", e)

# %% [markdown]
# If you want to have a mutable variant, for example to contain the right hand
# side of an ODE, use `struct_SX`:

# %%
rhs = struct_SX(states)
rhs["x"] = states["x"] ** 2
rhs["y"] = states["y"] * states["x"]
rhs["q"] = -states["q"]
print(rhs.cat)

# %% [markdown]
# Alternatively, you can supply the expressions at definition time:

# %%
x, y, q = states[...]
rhs = struct_SX([entry("x", expr=x**2), entry("y", expr=x * y), entry("q", expr=-q)])

print(rhs.cat)

# %% [markdown]
# One can also construct symbolic MX structures:

# %%
V = struct_symMX(shooting)
print(V)

# %% [markdown]
# The catted version is one single MX from which all entries are derived:

# %%
print(V.cat)
print(V.shape)
print(V["X", 0, -1, "y"])

# %% [markdown]
# Similar to `struct_SX`, we have `struct_MX`:

# %%
V = struct_MX(
    [
        (
            entry(
                "X", expr=[[MX.sym("x", 6) ** 2 for j in range(3)] for i in range(5)]
            ),
            entry("U", expr=[-MX.sym("u") for i in range(4)]),
        )
    ]
)

# %% [markdown]
# By default, the `struct_symSX` structure constructor will create new `SX.sym`s.
# To recycle one that is already available, use the `sym` argument:

# %%
qsym = SX.sym("quaternion", 4)
states = struct_symSX(["x", "y", entry("q", sym=qsym)])
print(states.cat)

# %% [markdown]
# The `sym` feature is not available for `struct_MX`, since it will construct
# one parent MX.

# %% [markdown]
# ## More powerIndex
#
# As illustrated before, powerIndex allows slicing
# %%
print(init["X", :, :, "x"])

# %% [markdown]
# The `repeated` method duplicates its argument a number of times such that it
# matches the length that is needed at the lhs:

# %%
init["X", :, :, "x"] = repeated(list(range(3)))
print(init["X", :, :, "x"])

# %% [markdown]
# Callables/functions can be thrown into the powerIndex at any location.
# They operate on subresults obtained from resolving the remainder of the
# powerIndex:

# %%
print(init["X", :, lambda v: horzcat(*v), :, "x"])
print(init["X", lambda v: vertcat(*v), :, lambda v: horzcat(*v), :, "x"])
print(init["X", blockcat, :, :, "x"])

# %% [markdown]
# Set all quaternions to 1,0,0,0:

# %%
init["X", :, :, "q"] = repeated(repeated(DM([1, 0, 0, 0])))

# %% [markdown]
# `{}` can be used in the powerIndex to expand into a dictionary once:

# %%
init["X", :, 0, {}] = repeated({"y": 9})
print(init["X", :, 0, {}])

# %% [markdown]
# Lists can be used in powerIndex in both list context or dict context:

# %%
print(shooting["X", [0, 1], [0, 1], "x"])
print(shooting["X", [0, 1], 0, ["x", "y"]])

# %% [markdown]
# `nesteddict` can be used to expand into a dictionary recursively:

# %%
print(init[nesteddict])

# %% [markdown]
# `...` will expand entries as an ordered list:

# %%
print(init["X", :, 0, ...])

# %% [markdown]
# If the powerIndex ends at the boundary of a structure, its catted version is
# returned:

# %%
print(init["X", 0, 0])

# %% [markdown]
# If the powerIndex is longer than what could be resolved as a structure, the
# remainder, *extraIndex*, is passed onto the resulting CasADi-matrix-type:

# %%
print(init["X", blockcat, :, :, "q", 0])
print(init["X", blockcat, :, :, "q", 0, 0])

# %% [markdown]
# ## shapeStruct and delegated indexing
#
# When working with covariance matrices, both the rows and columns relate to
# states:

# %%
states = struct(["x", "y", entry("q", repeat=2)])
V = struct_symSX(
    [
        entry("X", repeat=5, struct=states),
        entry("P", repeat=5, shapestruct=(states, states)),
    ]
)

# %% [markdown]
# `P` has a 4x4 shape:

# %%
print(V["P", 0])

# %% [markdown]
# Now we can use powerIndex-style in the extraIndex:

# %%
print(V["P", 0, ["x", "y"], ["x", "y"]])

# %% [markdown]
# There is a problem when we wish to use the full potential of powerIndex in
# these extraIndices. The following is in fact invalid Python syntax:
#
# ```python
# V["P", 0, ["q", :], ["q", :]]
# ```
#
# We resolve this by using delegater objects `index`/`indexf`:

# %%
print(V["P", 0, indexf["q", :], indexf["q", :]])

# %% [markdown]
# Of course, in this basic example, also the following would be allowed:

# %%
print(V["P", 0, "q", "q"])

# %% [markdown]
# ## Prefixing
#
# The `prefix` attribute allows you to create shorthands for long powerIndices:

# %%
states = struct(["x", "y", "z"])
V = struct_symSX([entry("X", repeat=[4, 5], struct=states)])
num = V()

# %% [markdown]
# Consider the following statements:

# %%
num["X", 0, 0, "x"] = 1
num["X", 0, 0, "y"] = 2
num["X", 0, 0, "z"] = 3

# %% [markdown]
# Note the common part `["X", 0, 0]`.
# We can pull this apart with `prefix`:

# %%
initial = num.prefix["X", 0, 0]
initial["x"] = 1
initial["y"] = 2
initial["z"] = 3

# %% [markdown]
# This is equivalent to the longer statements above.

# %% [markdown]
# ## Helper constructors
#
# If you work with `Simulator`, `ControlSimulator`, you typically end up
# with wanting to index a DM that is n x N,
# with n the size of a statespace and N an arbitrary integer:

# %%
states = struct(["x", "y", "z"])

# We artificially construct here a DM that could be a Simulator output.
output = DM.zeros(states.shape[0], 8)

# %% [markdown]
# The helper construct is `repeated` here. Instead of `states(output)`, we have:

# %%
outputs = states.repeated(output)

# %% [markdown]
# Now we have an object that supports powerIndexing:

# %%
outputs[-1] = DM([1, 2, 3])
outputs[:, "x"] = list(range(8))
print(output)
print(outputs[5, {}])

# %% [markdown]
# Next we represent the `squared` helper construct.
# Imagine we somehow obtain a matrix that represents covariance:

# %%
P0 = DM.zeros(states.shape[0], states.shape[0])

# %% [markdown]
# We can conveniently access it as follows:

# %%
P = states.squared(P0)
P["x", "y"] = 2
P["y", "x"] = 3
print(P0)

# %% [markdown]
# `P` itself is a rather queer object:

# %%
print(P)

# %% [markdown]
# You can access its contents with a call:

# %%
print(P())

# %% [markdown]
# But often, it will behave like a DM transparently:

# %%
print(P0 + P)

# %% [markdown]
# Next we represent the `squared_repeated` helper construct.
# Imagine we somehow obtain a matrix that represents a horizontal concatenation
# of covariances:

# %%
P0 = horzcat(
    DM.zeros(states.shape[0], states.shape[0]),
    DM.ones(states.shape[0], states.shape[0]),
)

# %% [markdown]
# We can conveniently access it as follows:

# %%
P = states.squared_repeated(P0)
P[0, "x", "y"] = 2
P[:, "y", "x"] = 3
print(P0)

# %% [markdown]
# Finally, we present the `product` helper construct:

# %%
controls = struct(["u", "v"])
J0 = DM.zeros(states.shape[0], controls.shape[0])
J = states.product(controls, J0)
J[:, "u"] = 3
J[["x", "z"], :] = 2
print(J())

# %% [markdown]
# ## Saving and loading
#
# It is possible to save and load some types of structures.
# Supported types are pure structures (the ones created with `struct`) and
# numeric structures.
#
# Saving:
#
# ```python
# mystructure.save("myfilename")
# ```
#
# Loading:
#
# ```python
# struct_load("myfilename")
# ```
#
# or:
#
# ```python
# import pickle
# with open("myfilename", "rb") as f:
#     mystructure = pickle.load(f)
# ```
