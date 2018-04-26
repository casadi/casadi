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
import glob
import subprocess
import os
import sys

# Python command
python_exec = sys.argv[1] if len(sys.argv)>1 else 'python'

# Process all snippets in the current directory
for f in glob.glob("pytex_*.py"):
  logfilename = f[:-3]+'.log'
  logfile = file(logfilename,'w')

  # Execute snippet
  p = subprocess.Popen([python_exec,f],stdout=logfile)
  p.wait()

  logfile = file(logfilename,'r')

  outfile = False
  contents = False
  for l in logfile.readlines():
    if l.startswith('pytex snippet::'):
      num=l[len('pytex snippet::'):].rstrip()

      # Remove logfiles that are empty
      #if outfile and not(contents):
      #  outfile.close()
      #  os.remove(outfile.name)

      outfile = file(f[:-3]+'_' + num + '.log','w')
      contents = False
    else:
      if outfile:
        contents = True
        outfile.write(l)
  os.remove(logfilename)
  # Remove logfiles that are empty
  #if outfile and not(contents):
  #  outfile.close()
  #  os.remove(outfile.name)
