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
function [s,v,m] = rocket(x)
[u,dt]=get_u_and_dt(x);
global par;

% Calculate acceleration for all points
a = nan(100,1);
for i=1:10
    a( 10*(i-1)+1:10*(i-1)+10) = u(i);
end

% Loop over all k to get s and v and the endpoints
s = nan(101,1); s(1) = par.s1;
v = nan(101,1); v(1) = par.v1;
m = nan(101,1); m(1) = par.m1;

for k=1:100
    s(k+1) = s(k) + dt*v(k);
    v(k+1) = v(k) + dt / m(k) * (a(k) - par.alpha * v(k)^2);
    m(k+1) = m(k) - dt * par.beta * a(k)^2;
    assert(m(k+1)>=0); % error message if negative mass    
end

% Plot the trajectory
if ~isempty(par) && isfield(par,'plot_traj') && par.plot_traj==true    
    t = linspace(0,100*dt,101);

     %dirty: 
    u = repmat(u',10,1); u = u(:); u(101) = u(100);
    clf;
    subplot(2,2,1); plot(t,s); title('Position'); grid on;
    subplot(2,2,2); plot(t,v); title('Velocity'); grid on;
    subplot(2,2,3); plot(t,m); title('Mass');     grid on;
    subplot(2,2,4); plot(t,u); title('Control');  grid on;
end