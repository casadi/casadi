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
addpath('../../../build/experimental/greg');

opt_ms  = spring_multiple_shooting_out();
% opt_lqr = spring_lqr_out();
opt_lqr = spring_ddp_out();

figure(1)
subplot(221)
hold off
plot( opt_ms.time, opt_ms.states.x, '.' )
hold on
plot( opt_lqr.time, opt_lqr.states.x, 'r')
xlabel('time')
title('x')

subplot(223)
hold off
plot( opt_ms.time, opt_ms.states.v, '.' )
hold on
plot( opt_lqr.time, opt_lqr.states.v, 'r' )
xlabel('time')
title('v')

subplot(222)
hold off
plot( opt_ms.time, opt_ms.actions.u, '.' )
hold on
plot( opt_lqr.time, opt_lqr.actions.u, 'r' )
xlabel('time')
title('u')
