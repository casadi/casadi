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
addpath('defs')

% Create NLP: Solve the Rosenbrock problem:
%     minimize    x^2 + 100*z^2
%     subject to  z + (1-x)^2 - y == 0
x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');
v = [x;y;z];
f = x^2 + 100*z^2;
g = z + (1-x)^2 - y;
nlp = struct('x', v, 'f', f', 'g', g);

mycallback = NlpCallback('mycallback',3, 1)
opts = struct;
opts.iteration_callback = mycallback;

% Create IPOPT solver object
solver = nlpsol('solver', 'ipopt', nlp, opts);

x0 = [2.5 3.0 0.75]';
% Solve the NLP
res = solver('x0' , x0,... % solution guess
             'lbx', -inf,...           % lower bound on x
             'ubx',  inf,...           % upper bound on x
             'lbg',    0,...           % lower bound on g
             'ubg',    0);             % upper bound on g
 
data = mycallback.data;

assert(size(data,2)>5)
assert(all(data(:,end)==full(res.x)))
assert(all(data(:,1)==x0))

