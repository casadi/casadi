import sys
src = sys.argv[1]
import os

from testsuite import TestSuite

def setdummybackend(dir):
  """
  Set the matplotlib backend to a dummy value, so that no graphs are generated
  """
  f=file(dir+'/matplotlibrc','w')
  f.write('backend : Template')
  f.close()
  
def removedummybackend(dir):
  os.remove(dir+'/matplotlibrc')


t = TestSuite(dirname=src,
  suffix="py",
  preRun=setdummybackend,
  postRun=removedummybackend,
  command = lambda dir,fn:  ["python", fn],
  skipdirs=[".svn","ctemplate"]
  )
  
t.run()
