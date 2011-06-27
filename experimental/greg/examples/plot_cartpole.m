addpath('../../../build/experimental/greg');

opt = cartpole_out();

figure(1)
subplot(231)
plot( opt.time, opt.states.x )
xlabel('time')
title('x')

subplot(234)
plot( opt.time, opt.states.theta )
xlabel('time')
title('theta')

subplot(232)
plot( opt.time, opt.states.vx )
xlabel('time')
title('vx')

subplot(235)
plot( opt.time, opt.states.vtheta )
xlabel('time')
title('vtheta')

subplot(233)
plot( opt.time, opt.actions.u )
xlabel('time')
title('u')


figure(2)
hold off
plot3( opt.outputs.cart_x, opt.time, opt.outputs.cart_y, 'b' );
hold on
plot3( opt.outputs.bob_x,  opt.time, opt.outputs.bob_y,  'r' );
xlabel('x')
ylabel('time')
zlabel('y')
grid on