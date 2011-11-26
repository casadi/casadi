import sys
src = sys.argv[1]
import os

"""
Test cpp files that are made with CMAKE and put in the /build/bin directory
"""

from testsuite import TestSuite

from subprocess import *

t = TestSuite(dirname=src,
  suffix="cpp",
  command = lambda dir,fn:  ['./bin/'+fn.replace('.cpp','')],
  workingdir = lambda dir: '../build', 
  skipdirs=[".svn","ctemplate"],
  args=sys.argv[2:]
  )
  
t.run()




