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
%! Load casadi
casadi

x = SX.sym("x");     % A scalar (1-by-1 matrix) symbolic primitive
y = SX.sym("y",2);   % A vector (n-by-1 matrix) symbolic primitive
z = SX.sym("z",2,3); % An n-by-m matrix symbolic primitive

in = {x y}; % function inputs
disp('Function outputs are:')
out = {x,y,[x x; x x],y*x,0}

f = SXFunction(in,out);
f.init();

%! f now has two inputs and a 4 outputs:
number_in = f.getNumInputs()
number_out = f.getNumOutputs()

%! The outputs has the following string representation.
%! Note how all elements of out have been converted to SX by
%! automatic typecasting functionality

for i = 0:3
  disp(["Output:"])
  disp(f.outputExpr(i))
end
