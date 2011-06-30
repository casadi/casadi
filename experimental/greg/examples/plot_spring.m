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
