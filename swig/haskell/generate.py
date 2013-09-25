from lxml import etree
import sys

e = etree.parse(sys.argv[1])
r = e.getroot()

print filter(lambda x: 'IOSchemeVector' not in x,[a.attrib["value"] for a in r.findall('*//class/attributelist/attribute[@name="name"]')])

f = file('casadiTree.hs','w')

f.write("hello")
