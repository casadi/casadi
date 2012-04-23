%
%     This file is part of CasADi.
% 
%     CasADi -- A symbolic framework for dynamic optimization.
%     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
% 
%     CasADi is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 3 of the License, or (at your option) any later version.
% 
%     CasADi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public
%     License along with CasADi; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
% 
% 
% Example matrix from Davis p. 39
A = sparse(11,11);
A(3,2) = 1;
A(6,1) = 1;
A(6,4) = 1;
A(7,1) = 1;
A(8,2) = 1;
A(8,5) = 1;
A(9,6) = 1;
A(10,3) = 1;
A(10,4) = 1;
A(10,6) = 1;
A(10,8) = 1;
A(11,3) = 1;
A(11,5) = 1;
A(11,7) = 1;
A(11,8) = 1;
A(11,10) = 1;

% Make symmetric
A = A+A';

% Get the elimination tree
etr = etree(A,'col')
