% Car race along a track
% ----------------------
% An optimal control problem (OCP),
% solved with direct multiple-shooting.
%
% For more information see: http://labs.casadi.org/OCP

N = 100; % number of control intervals

opti = casadi.Opti(); % Optimization problem

% ---- decision variables ---------
X = opti.variable(2,N+1); % state trajectory
pos   = X(1,:);
speed = X(2,:);
U = opti.variable(1,N);   % control trajectory (throttle)
T = opti.variable();      % final time

% ---- objective          ---------
opti.minimize(T); % race in minimal time

% ---- dynamic constraints --------
f = @(x,u) [x(2);u-x(2)]; % dx/dt = f(x,u)

dt = T/N; % length of a control interval
for k=1:N % loop over control intervals
   % Runge-Kutta 4 integration
   k1 = f(X(:,k),         U(:,k));
   k2 = f(X(:,k)+dt/2*k1, U(:,k));
   k3 = f(X(:,k)+dt/2*k2, U(:,k));
   k4 = f(X(:,k)+dt*k3,   U(:,k));
   x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4); 
   opti.subject_to(X(:,k+1)==x_next); % close the gaps
end

% ---- path constraints -----------
limit = @(pos) 1-sin(2*pi*pos)/2;
opti.subject_to(speed<=limit(pos)); % track speed limit
opti.subject_to(0<=U<=1);           % control is limited

% ---- boundary conditions --------
opti.subject_to(pos(1)==0);   % start at position 0 ...
opti.subject_to(speed(1)==0); % ... from stand-still 
opti.subject_to(pos(N+1)==1); % finish line at position 1

% ---- misc. constraints  ----------
opti.subject_to(T>=0); % Time must be positive

% ---- initial values for solver ---
opti.set_initial(speed, 1);
opti.set_initial(T, 1);

% ---- solve NLP              ------
opti.solver('ipopt'); % set numerical backend
sol = opti.solve();   % actual solve

% ---- post-processing        ------

figure
hold on
plot(sol.value(speed));
plot(sol.value(pos));
plot(limit(sol.value(pos)),'r--');
stairs(1:N,sol.value(U),'k');
legend('speed','pos','speed limit','throttle','Location','northwest')


figure
spy(sol.value(jacobian(opti.g,opti.x)))
figure
spy(sol.value(hessian(opti.f+opti.lam_g'*opti.g,opti.x)))

