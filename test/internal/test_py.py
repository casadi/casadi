import sys
src = sys.argv[1]
import os
import shutil

from testsuite import TestSuite

def setdummybackend(dir):
  """
  Set the matplotlib backend to a dummy value, so that no graphs are generated.
  
  Was originally using 'backend : Template' in matplotlibrc, but this leaked memory
  """
  shutil.copy('internal/pylab.py',dir)
  
def removedummybackend(dir):
  os.remove(dir+'/pylab.py')


t = TestSuite(dirname=src,
  suffix="py",
  preRun=setdummybackend,
  postRun=removedummybackend,
  command = lambda dir,fn:  ["python", fn],
  skipdirs=[".svn","ctemplate"],
    args=sys.argv[2:]
  )
  
t.run()
