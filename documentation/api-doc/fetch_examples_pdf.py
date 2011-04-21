"""
Grabs the pdf-compiled examples and puts them in the API tree.
"""

src = '../examples/'
xml = 'XML/index.xml'

import re
import shutil
from lxml import etree
import os

xmlData = etree.parse(xml)

examples = xmlData.findall("//compound[@kind='example']")


for example in examples:
  refid = example.attrib['refid']
  path = '/'.join(refid.split("/")[0:2])+'/'
  example.attrib['refid']
  name=example.findtext("name")
  m=re.search('^(.*)\.py$',name)
  if m:
    try:
      os.makedirs(os.path.dirname('html/'+path+m.group(1) + '.pdf'))
    except:
      pass
    shutil.copyfile(src+m.group(1) + '.pdf','html/'+path+m.group(1) + '.pdf')
    print m.group(1) + '.pdf'
