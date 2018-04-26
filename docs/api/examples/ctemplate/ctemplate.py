#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
import re
import os

builddir = 'build/'
topdir = '../../../../'
makedir = 'docs/api/examples/ctemplate/CMakeFiles/ctemplate.dir/'

cextension = 'cc'

absbuilddir = os.path.abspath(topdir+builddir)

# Write a linker script
link = file(topdir+builddir+makedir+'link.txt','r')
linker = link.readline()
linker = re.sub(' \.\./\.\./\.\./\.\.',' '+ absbuilddir,linker)
linker = re.sub('CMakeFiles/ctemplate.dir/ctemplate(\.cpp)?\.o','"$1.o"',linker)
linker = re.sub('-o (.*?) ','-o "$1.run" ',linker)

linkerscript = file('linker.sh','w')
linkerscript.write('#!/bin/bash\n')
linkerscript.write(linker)
linkerscript.close()
link.close()

# Write a compiler script
myflags=''
mydefines=''
myincludes=''

flags = file(topdir+builddir+makedir+'flags.make','r')
for l in flags:
  m = re.search("CXX_FLAGS = (.*)",l)
  if m:
    myflags = m.group(1)
  m = re.search("CXX_DEFINES = (.*)",l)
  if m:
    mydefines = m.group(1)
  m = re.search("CXX_INCLUDES = (.*)",l)
  if m:
    myincludes = m.group(1)
compiler = linker.split(' ')[0] + ' ' + myflags +" "+ mydefines +" "+ myincludes+ ' -o "$1.o" -c "$1.%s"' % cextension
flags.close()
compilerscript = file('compiler.sh','w')
compilerscript.write('#!/bin/bash\n')
compilerscript.write(compiler)
compilerscript.close()
link.close()

print compiler
print linker

print "Done"

