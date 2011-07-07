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