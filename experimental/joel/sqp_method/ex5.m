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
%% Exercise 5 in Numerical Optimization
% Joel Andersson, ESAT, K.U. LEUVEN
% November 2009
% joel.andersson@esat.kuleuven.be

clear all

% make a parameter structure and set parameters
global par;
par = struct; 
% Starting position
par.s1 = 0;   par.v1 =0;   par.m1=1;

% Final position
par.s101=10; par.v101=0;

% Time step
par.dt=0.1;
par.alpha=0.05; par.beta=0.1;

% Constraints
par.vmax=1.3; par.umin = -1;  par.umax = 0.5;

% starting point
u0=0.4*ones(10,1); 

%% minimize fuel
par.obj = 'fuel';  x0=u0;
par.plot_traj = false; % no plotting
xopt_fuel = sqpmethod(@f,@g,[],x0);

return

figure(1);
set(gcf,'Name','Minimum fuel')
par.plot_traj = true; % turn on plotting
rocket(xopt_fuel);

%% minimize time
par.obj = 'time'; x0=[u0;par.dt];

par.plot_traj = false; % no plotting
xopt_time = sqpmethod(@f,@g,@h,x0);

figure(2);
set(gcf,'Name','Minimum time')
par.plot_traj = true; % turn on plotting
rocket(xopt_time);
