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
src = 'src'
out = 'pdf'

import sys
args=sys.argv[1:]

try:
	from pyreport import *
except:
	raise Exception("Unable to compile.\n\nPlease install pyreport with 'sudo easy_install pyreport' to compile the tutorials")

try:
	import casadi
except:
	raise Exception("Couldn't load casadi's python module. Make sure you have compiled them with 'make python' and make sure you have set your PYTHONPATH to link to casadi's build_dir/swig_interface and build_dir/lib")


import os

def pyrep(indir,name,outdir):
	try:
		os.makedirs(outdir)
	except:
		pass
	(base,ext)=os.path.splitext(name)
	pdfname = base + ".pdf"
	try:
	  os.remove(pdfname)
	except:
	  pass
	
	cmd= 'cd %s && grep -v "^#[^\!]" %s > temp.py && pyreport -d -l temp.py -o %s' % (indir,name,pdfname)
	print cmd
	os.system(cmd)
	try:
	  os.remove(os.path.join(indir,"temp.py"))
	except:
	  pass

	try:
		os.rename(os.path.join(indir,pdfname), os.path.join(outdir,pdfname))
	except OSError:
		if os.path.exists(os.path.join(outdir,pdfname)):
			os.remove(os.path.join(outdir,pdfname))
			os.rename(os.path.join(indir,pdfname), os.path.join(outdir,pdfname))
		print "Couldn't move to output directory."
		pass 
		
for root, dirs, files in os.walk(src):
	for name in files:
		if '.svn' in dirs:
			dirs.remove('.svn') # don't visit svn folders
		if len(args)>0 and all(map(lambda a: name.find(a)==-1,args)):
			continue
		if name.endswith('.py') and name!="temp.py":
			relroot=os.path.relpath(root,src)
			pyrep(root,name,os.path.join(out,relroot))
			
			
			
