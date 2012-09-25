##
##	This file is part of qpOASES.
##
##	qpOASES -- An Implementation of the Online Active Set Strategy.
##	Copyright (C) 2007-2012 by Hans Joachim Ferreau, Andreas Potschka,
##	Christian Kirches et al. All rights reserved.
##
##	qpOASES is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	qpOASES is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with qpOASES; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##



INTRODUCTION
============

qpOASES is an open-source C++ implementation of the recently proposed 
online active set strategy (see [1], [2]), which was inspired by important 
observations from the field of parametric quadratic programming. It has 
several theoretical features that make it particularly suited for model 
predictive control (MPC) applications.

The software package qpOASES implements these ideas and has already been 
successfully used for many academic and industrial real-world applications, 
e.g. closed-loop control of a real-world Diesel engine testbench (see [3]). 
Recently, a couple of numerical modifications (as proposed in [4]) have been 
introduced that greatly increase qpOASES's reliability when solving
semi-definite, ill-posed or degenerated convex QPs.


References:

[1] H.J. Ferreau. An Online Active Set Strategy for Fast Solution of 
Parametric Quadratic Programs with Applications to Predictive Engine Control. 
Diplom thesis, University of Heidelberg, 2006.

[2] H.J. Ferreau, H.G. Bock, M. Diehl. An online active set strategy to 
overcome the limitations of explicit MPC. International Journal of Robust 
and Nonlinear Control, 18 (8), pp. 816-830, 2008.

[3] H.J. Ferreau, P. Ortner, P. Langthaler, L. del Re, M. Diehl. Predictive 
Control of a Real-World Diesel Engine using an Extended Online Active Set 
Strategy. Annual Reviews in Control, 31 (2), pp. 293-301, 2007.

[4] A. Potschka, C. Kirches, H.G. Bock, and J.P. Schloeder. Reliable solution 
of convex quadratic programs with parametric active set methods. Technical
report, IWR, Heidelberg University, Germany, November 2010.



GETTING STARTED
===============

1. For installation, usage and additional information on this software package 
   see the qpOASES User's Manual located at doc/manual.pdf or check its
   source code documentation!


2. The file LICENSE.txt contains a copy of the GNU Lesser General Public 
   License. Please read it carefully before using qpOASES!


3. The whole software package can be downloaded from 

       http://www.qpOASES.org/

   On this webpage you will also find a list of frequently asked questions.



CONTACT THE AUTHORS
===================

If you have got questions, remarks or comments on qpOASES, send them to

        support@qpOASES.org

Alternatively, you can directly contact one of the main authors:

        Hans Joachim Ferreau, joachim.ferreau@esat.kuleuven.be
        Andreas Potschka, potschka@iwr.uni-heidelberg.de
        Christian Kirches, christian.kirches@iwr.uni-heidelberg.de

Also bug reports and source code extensions are most welcome!



##
##	end of file
##
