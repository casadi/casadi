import sys
src = sys.argv[1]
import os

from testsuite import TestSuite

from subprocess import *

t = TestSuite(dirname=src,
  workingdir = lambda dir : os.path.join(dir,'..'),
  command = lambda dir,fn:  ['bin/'+fn],
  inputs = {'det_minor': "5"},
    args=sys.argv[2:]
  )
  
t.run()




