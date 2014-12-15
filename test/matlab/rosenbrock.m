%
%     This file is part of CasADi.
%
%     CasADi -- A symbolic framework for dynamic optimization.
%     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
%                             K.U. Leuven. All rights reserved.
%     Copyright (C) 2011-2014 Greg Horn
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
% Load CasADi
import casadi.*
 
% Create NLP: Solve the Rosenbrock problem:
%     minimize    x^2 + 100*z^2
%     subject to  z + (1-x)^2 - y == 0
x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');
v = SX.casadi_vertcat({x,y,z}); % workaround!
f = x^2 + SX(100)*z^2;
g = z + (SX(1)-x)^2 - y;
nlp = SXFunction(nlpIn('x',v),nlpOut('f',f','g',g));
 
% Create IPOPT solver object
solver = NlpSolver('ipopt', nlp);
solver.init();
 
% Solution guess [2.5,3.0,0.75]
x0 = DMatrix([2.5;3.0;0.75]);
solver.setInput(x0,'x0');
 
% Set variable bounds
solver.setInput(0, 'lbg');
solver.setInput(0, 'ubg');
 
% Solve the NLP
solver.evaluate()
 
% Print the solution
f_opt = solver.getOutput('f')          % >> DMatrix(0)
x_opt = solver.getOutput('x')          % >> DMatrix([0, 1, 0]) 
lam_x_opt = solver.getOutput('lam_x')  % >> DMatrix([0, 0, 0])
lam_g_opt = solver.getOutput('lam_g')  % >> DMatrix(0)
 
% Unload CasADi (important!)
clear
