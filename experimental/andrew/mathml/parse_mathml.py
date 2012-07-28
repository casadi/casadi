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
	




