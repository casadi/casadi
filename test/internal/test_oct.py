import sys
src = sys.argv[1]
import os

from testsuite import TestSuite

t = TestSuite(dirname=src,
  suffix="m",
  command = lambda dir,fn:  ["octave",'--no-init-file','-p', os.getcwd() + '/../build/lib', fn],
  skipdirs=[".svn","ctemplate"],
  allowable_returncodes=[127],
    args=sys.argv[2:]
  )
  
t.run()
