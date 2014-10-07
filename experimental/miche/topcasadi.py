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
from casadi import *
from numpy import *

# Static model parameters
nelx = 12
nely = 4
rmin = 3
penal = 3
volfrac = 0.5

class mbb2d:
	def __init__(self,nelx=120,nely=40,rmin=3,penal=3,volfrac=0.5):
		self.nelx   = nelx
		self.nely   = nely
		self.rmin   = rmin
		self.penal  = penal
		self.volfrac= volfrac

		self.nelem = nelem = nelx*nely  #total number of design variables

		# Prepare finite element analysis
		E0 = 1.0;
		Emin = 1.0E-9;
		nu = 0.3;

		ke = array([(0.5-nu/6) ,  (0.125+nu/8) , (-0.25-nu/12), (-0.125+3*nu/8),  (-0.25+nu/12), (-0.125-nu/8), (nu/6), (0.125-3*nu/8)])
		ke *= 1/(1-nu**2)

		# This should can be rephrased with a slice

		KE = array([[ke[0], ke[1], ke[2], ke[3], ke[4], ke[5], ke[6], ke[7]],
		[ke[1], ke[0], ke[7], ke[6], ke[5], ke[4], ke[3], ke[2]],
		[ke[2], ke[7], ke[0], ke[5], ke[6], ke[3], ke[4], ke[1]],
		[ke[3], ke[6], ke[5], ke[0], ke[7], ke[2], ke[1], ke[4]],
		[ke[4], ke[5], ke[6], ke[7], ke[0], ke[1], ke[2], ke[3]],
		[ke[5], ke[4], ke[3], ke[2], ke[1], ke[0], ke[7], ke[6]],
		[ke[6], ke[3], ke[4], ke[1], ke[2], ke[7], ke[0], ke[5]],
		[ke[7], ke[2], ke[1], ke[4], ke[3], ke[6], ke[5], ke[0]]])


		ndofs = 2*(nelx+1)*(nely+1)
		nfixeddofs = nely+2
		self.nfreedofs = nfreedofs = ndofs-nfixeddofs

		alldofs = r_[0:ndofs]

		fixeddofs = zeros(nfixeddofs) 
		for i in range(nfixeddofs):
			fixeddofs[i]=2*i
			
		fixeddofs[nely+1] = 2*(nelx+1)*(nely+1)-1

		freedofs = r_[nfreedofs]

		freedofs = list(set(alldofs) - set(fixeddofs))

		# (reduced) external force vector 
		F = SX.sparse(nfreedofs,1)
		
		F[0,0] = -1;

		#Filter matrix (row compressed format)

		crmin = int(ceil(rmin))
		nh =  (2 * (crmin-1)+1)*(2 * (crmin-1)+1)
		nH = nelem*nh
		self.pH = pH = r_[0:nelem+1]
		self.jH = jH = r_[0:nH]
		sH = r_[0:nH]
		Hs = r_[0:nelem]

		pH[0]=0

		k = 0
		print "Construction of H"
		for ielx in range(nelx):
			for iely in range(nely):
				ielem = iely + nely * ielx
				xmin = max(ielx-crmin+1,0)
				xmax = min(ielx+crmin-1,nelx-1)
				ymin = max(iely-crmin+1,0)
				ymax = min(iely+crmin-1,nely-1)
				fsum = 0.0
				for ixloc in range(xmin,xmax+1):
					for iyloc in range(ymin,ymax+1):	
						jH[k] = iyloc+ nely*ixloc;
						rx = ixloc-ielx
						ry = iyloc-iely
						hf = max(rmin - sqrt( rx * rx + ry * ry), 0.0)
						fsum += hf
						sH[k] = hf
						k+=1
				pH[ielem+1]=k
				Hs[ielem] = fsum
		nH = k

    
		self.sH = SX(sH)
		self.Hs = SX(Hs)
		H = SX()
		
	def dfilter(self,x):
		""" Density Filter """
		xf = SX(nelem+nfreedofs,1,0)
		for ielem in range(self.nelem):
			for t in range(self.pH[ielem],self.pH[ielem+1]):
				kelem = self.jH[t]
				xf[ielem,0] += self.sH[t,0] * x[kelem,0]/self.Hs[ielem,0]
		# copy state variables
		xf[nelem:,0] = x[nelem:,0]
		return xf
	
		

# Create a model instance
mbb = mbb2d(nelx,nely,rmin,penal,volfrac)

nelem = mbb.nelem
nfreedofs = mbb.nfreedofs

x = vertcat([ssym("x",nelem),ssym("u",nfreedofs)])
for i in range(nelem,nelem+nfreedofs-2,2):
	x[i]=-1;
	
print "filtering of H"
			
xf = mbb.dfilter(x)

print xf
