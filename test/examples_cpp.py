src = '../documentation/examples'

import os
from subprocess import *

stats={'numtests':0,'numfails':0}

p=Popen(['make', 'cpp'],cwd=src,stdout=PIPE, stderr=PIPE)
stdoutdata, stderrdata = p.communicate()
if not(p.returncode==0):
  print stdoutdata
  print stderrdata
  raise Exception("Was unable to compile.")
  
def test(dir,fn):
  print fn
  stats['numtests']+=1
  p=Popen(['./'+fn],cwd=dir,stdout=PIPE, stderr=PIPE)
  stdoutdata, stderrdata = p.communicate()
  if not(p.returncode==0):
    stats['numfails']+=1
    print "In %s, %S failed:" % (dir,fn)
    print "="*30
    print p.returncode
    print stdoutdata
    print stderrdata
    print "="*30
    
for root, dirs, files in os.walk(src):
  for name in files:
    if '.svn' in dirs:
      dirs.remove('.svn') # don't visit svn folders
    if 'ctemplate' in dirs:
      dirs.remove('ctemplate') # don't visit ctemplate folders 
    if name.endswith('.run'):
      test(root,name)
  
print "Ran %d tests, %d fails." % (stats['numtests'],stats['numfails'])
