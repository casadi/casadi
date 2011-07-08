import sys
src = sys.argv[1]
import os

from testsuite import TestSuite

t = TestSuite(dirname=src,
  suffix="m",
  command = lambda dir,fn:  ["octave", fn],
  skipdirs=[".svn","ctemplate"],
  allowable_returncodes=[127],
    args=sys.argv[2:]
  )
  
t.run()
