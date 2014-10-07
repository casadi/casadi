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
#!/usr/bin/env python
# Parse a MathML (xml) file for content

import xml.etree.ElementTree as xml
import casadi

operatorMappings = [('cos',casadi.cos),
					('sin', casadi.sin),
					('tan', casadi.tan),
					('plus', lambda x,y: x+y),
					('minus', lambda x,y: x-y),
					('times', lambda x,y: x*y) ]

def parse_mathml(filename):
	pass

if __name__=='__main__':
	filename = 'example1.xml'
	#parse_mathml(filename)

	tree = xml.parse(filename)
	root = tree.getroot()
	semanticscontents = root.getchildren()[0].getchildren()
	expressionroots = []
	for el in semanticscontents:
		if el.attrib.has_key('encoding') and el.attrib['encoding']=='MathML-Content':
			expressionroots.append(el)
	assert len(expressionroots)==1, "Can only handle 1 expression per file!"
	expressionroot = expressionroots[0]
	expressionroot = expressionroot.getchildren()[0]
	




