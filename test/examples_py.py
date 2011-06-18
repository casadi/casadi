import sys

src = sys.argv[1]

import os
import sys
from subprocess import *

def setdummybackend(dir):
  """
  Set the matplotlib backend to a dummy value, so that no graphs are generated
  """
  f=file(dir+'/matplotlibrc','w')
  f.write('backend : Template')
  f.close()
  
def removedummybackend(dir):
  os.remove(dir+'/matplotlibrc')

stats={'numtests':0,'numfails':0}

def test(dir,fn):
  print fn
  stats['numtests']+=1
  p=Popen(['python', fn],cwd=dir,stdout=PIPE, stderr=PIPE)
  stdoutdata, stderrdata = p.communicate()
  if not(p.returncode==0):
    stats['numfails']+=1
    print "In %s, %s failed:" % (dir,fn)
    print "="*30
    print p.returncode
    print stdoutdata
    print stderrdata
    print "="*30
    
for root, dirs, files in os.walk(src):
  setdummybackend(root)
  for name in files:
    if '.svn' in dirs:
      dirs.remove('.svn') # don't visit svn folders
    if 'ctemplate' in dirs:
      dirs.remove('ctemplate') # don't visit ctemplate folders 
    if name.endswith('.py'):
      test(root,name)
  removedummybackend(root)
  
print "Ran %d tests, %d fails." % (stats['numtests'],stats['numfails'])

if stats['numfails']>0:
  sys.exit(1)
