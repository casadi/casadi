import sys
src = sys.argv[1]
import os

from testsuite import TestSuite

from subprocess import *

print "Compiling"
p=Popen(['make', 'cpp'],cwd=src,stdout=PIPE, stderr=PIPE)
stdoutdata, stderrdata = p.communicate()
if p.returncode==0:
  print "Done compiling"
else:
  print stdoutdata
  print stderrdata
  raise Exception("Was unable to compile.")

t = TestSuite(dirname=src,
  suffix="run",
  command = lambda dir,fn:  ['./'+fn],
  skipdirs=[".svn","ctemplate"],
    args=sys.argv[2:]
  )
  
t.run()




