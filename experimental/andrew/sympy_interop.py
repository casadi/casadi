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
#!/usr/bin/python
#
# Functions for interoperability between sympy and casadi
#
# Notes:
# All sympy Symbols that have the same string have the same address.
# In casadi, the string is only for printing, each time you create an
# expression it is a new object.

import sympy
import casadi

def traverse(node, casadi_syms, rootnode):
	#print node
	#print node.args
	#print len(node.args)

	if len(node.args)==0:
		# Handle symbols
		if(node.is_Symbol):
			return casadi_syms[node.name]

		# Handle numbers and constants
		if node.is_Zero:
			return 0
		if node.is_Number:
			return float(node)

	trig = sympy.functions.elementary.trigonometric
	if len(node.args)==1:
		# Handle unary operators
		child = traverse(node.args[0], casadi_syms, rootnode) # Recursion!
		if type(node) == trig.cos:
			return casadi.cos(child)
		if type(node) == trig.sin:
			return casadi.sin(child)
		if type(node) == trig.tan:
			return casadi.tan(child)

		if type(node) == trig.cosh:
			return casadi.cosh(child)
		if type(node) == trig.sinh:
			return casadi.sinh(child)
		if type(node) == trig.tanh:
			return casadi.tanh(child)

		if type(node) == trig.cot:
			return 1/casadi.tan(child)
		if type(node) == trig.acos:
			return casadi.arccos(child)
		if type(node) == trig.asin:
			return casadi.arcsin(child)
		if type(node) == trig.atan:
			return casadi.arctan(child)

	if len(node.args)==2:
		# Handle binary operators
		left = traverse(node.args[0], casadi_syms, rootnode) # Recursion!
		right = traverse(node.args[1], casadi_syms, rootnode) # Recursion!
		if node.is_Pow:
			return left**right
		if type(node) == trig.atan2:
			return casadi.arctan2(left,right)

	if len(node.args)>=2:
		# Handle n-ary operators
		child_generator = ( traverse(arg,casadi_syms,rootnode) 
								for arg in node.args )
		if node.is_Add:
			return reduce(lambda x, y: x+y, child_generator)
		if node.is_Mul:
			return reduce(lambda x, y: x*y, child_generator)

	if node!=rootnode:
		raise Exception("No mapping to casadi for node of type " 
				+ str(type(node)))

def sympy2casadi(node):
	sympy_syms = node.free_symbols  # a python set of sympy Symbols
	# A dictionary matching symbol names to instances of casadi Matrix<SXElement>
	casadi_syms = {sym.name:casadi.ssym(sym.name) for sym in sympy_syms}
	return traverse(node, casadi_syms, node)

if __name__ == "__main__":
	import sympy
	import casadi

	# Generate relatively easy expression
	a,b,c = sympy.symbols('a b c')
	y = (sympy.sin(a)+b)/c

	# Generate a relatively complex sympy expression
	y=a
	for i in range(4):
		y = sympy.sin(y)*y

	# Convert to casadi!
	y_sympy = y
	y_casadi = sympy2casadi(y_sympy)
	print y
	print y_casadi

