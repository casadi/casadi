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

opt  = raptor_out();

figure(1)
hold off
plot( opt.states.xh, opt.states.yh, 'r' )
hold on
plot( opt.states.x1, opt.states.y1, 'b' )
plot( opt.states.x2, opt.states.y2, 'b' )
plot( opt.states.x3, opt.states.y3, 'b' )
axis equal

figure(2)
subplot(211)
plot( opt.time, opt.actions.dtheta*180/pi )
title('dtheta (deg/s)')

subplot(212)
plot( opt.time, opt.states.theta*180/pi )
title('theta (deg/s)')

fprintf('survival time: %.8f\n', opt.time(end))