import re
import os

builddir = 'build/'
topdir = '../../../'
makedir = 'documentation/examples/ctemplate/CMakeFiles/ctemplate.dir/'

cextension = 'cc'

absbuilddir = os.path.abspath(topdir+builddir)

# Write a linker script
link = file(topdir+builddir+makedir+'link.txt','r')
linker = link.readline()
linker = re.sub(' \.\./\.\./\.\.',' '+ absbuilddir,linker)
linker = re.sub('CMakeFiles/ctemplate.dir/ctemplate.o','"$1.o"',linker)
linker = re.sub('-o (.*?) ','-o "$1.run" ',linker)

linkerscript = file('linker.sh','w')
linkerscript.write('#!/bin/bash\n')
linkerscript.write(linker)
linkerscript.close()
link.close()

# Write a compiler script
myflags=''
mydefines=''

flags = file(topdir+builddir+makedir+'flags.make','r')
for l in flags:
  m = re.search("CXX_FLAGS = (.*)",l)
  if m:
    myflags = m.group(1)
  m = re.search("CXX_DEFINES = (.*)",l)
  if m:
    mydefines = m.group(1)
    
compiler = linker.split(' ')[0] + ' ' + myflags + mydefines + ' -o "$1.o" -c "$1.%s"' % cextension
flags.close()
compilerscript = file('compiler.sh','w')
compilerscript.write('#!/bin/bash\n')
compilerscript.write(compiler)
compilerscript.close()
link.close()

print compiler
print linker

print "Done"

