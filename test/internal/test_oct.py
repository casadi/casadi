
import sys
src = sys.argv[1]
import os

from testsuite import TestSuite

path = os.getcwd() + '/../build/lib'
if 'CASADILIBDIR' in os.environ:
  path = os.environ['CASADILIBDIR']

t = TestSuite(dirname=src,
  suffix="m",
  command = lambda dir,fn:  ["octave",'--no-init-file','-p', path, fn],
  skipdirs=[".svn","ctemplate"],
  allowable_returncodes=[127],
  stderr_trigger = ["error"],
    args=sys.argv[2:]
  )
  
t.run()
