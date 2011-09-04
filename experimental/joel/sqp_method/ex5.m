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
