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
% SQPMETHOD:      SQP NLP solver
% This is a sequential quadratic programming solver for constrained
% nonlinear least squares problems of the form:
%                     minimize    f(x)
%                     subject to
%                                 g(x) == 0,
%                                 h(x) >= 0.
%
% It calculates Jacobians of the objective and constraint functions with
% the 'i-trick' method, therefore conjugate transposes (') must be replaced
% with non-conjugate transposes .' (dot-prime) or the 'transp(...)'
% function in the user supplied functions.
%
% The algorithm is a Quasi-Newton method with damped BFGS updating to that
% assures positive definitenes of the Hessian approximation. Line search is
% carried out via backtracking until with the Armijo condition applied to
% the T1 (in Nocedal phi1) merit function is satisfied.
%
% Syntax:
% [x_opt] = sqpmethod(f,g,h,x0)
%
% f:   function handles to the objective function
% g:   function handles to the equality constraint function
% h:   function handles to the inequality contraint function
% x0:  starting point for the iteration
%
% November 2008-2009
% Joel Andersson, ESAT, KULeuven
% joel.andersson@esat.kuleuven.be
%------------------------------------------------------

% Main function
function [x] = sqpmethod(ffun,gfun,hfun,x)
% Parameters in the algorithm
max_iter = 100; % maximum number of sqp iterations
toldx = 1e-12; % stopping criterion for the stepsize
tolgL = 1e-12; % stopping criterion for the lagrangian gradient
merit_mu = 0;  % current 'mu' in the T1 merit function

% Evaluate the functions
[Jfk,fk] = itrick(ffun,x);
[Jgk,gk] = itrick(gfun,x);
[Jhk,hk] = itrick(hfun,x);

% Get dimensions
m = length(gk); % Number of equality constraints
n = length(x);  % Number of variables
q = length(hk); % Number of inequality constraints

% Initial guess for the lagrange multipliers
lambda_k = zeros(m,1);
mu_k = zeros(q,1);

% Initial guess for the Hessian
Bk = eye(n);

% Header
disp(' k  nls | dx         gradL      eq viol    ineq viol')
k = 0;
while true    
    % Solve the QP subproblem
    % minimize:   1/2*p'*Bk*p + Jfk*p
    % subject to: Jgk*p + gk == 0
    %             Jhk*p + hk >= 0
    cvx_quiet(true); % Suppress output
    cvx_begin
    variable p(n)
    if ~isempty(gfun) dual variable lambda_hat;  end
    if ~isempty(hfun) dual variable mu_hat;      end

    minimize( 1/2*quad_form(p,Bk) + Jfk*p);
    subject to
    if ~isempty(gfun) lambda_hat : Jgk*p + gk == 0; end
    if ~isempty(hfun) mu_hat :     Jhk*p + hk >= 0; end
    cvx_end
        
    % Get the gradient of the Lagrangian
    gradL = Jfk';
    if ~isempty(gfun) gradL = gradL - Jgk'*lambda_hat; end
    if ~isempty(hfun) gradL = gradL - Jhk'*mu_hat;     end

    % Do a line search along p
    [tk,merit_mu,nlinsearch] = linesearch(ffun,gfun,hfun,x,...
        p,fk,gk,hk,Jfk,Jgk,Jhk,Bk,merit_mu);
    
    % Calculate the new step
    dx = p*tk;
    x = x + dx;
    if ~isempty(gfun) lambda_k = tk*lambda_hat + (1-tk)*lambda_k; end
    if ~isempty(hfun) mu_k     = tk*mu_hat + (1-tk)*mu_k;         end
    k = k+1;

    % Gather and print iteration information
    normdx = norm(dx); % step size
    normgradL = norm(gradL); % size of the Lagrangian gradient
    eq_viol = sum(abs(gk)); % equality constraint violation
    ineq_viol = sum(max(0,-hk)); % inequality constraint violation

    disp(sprintf('%3d %3d |%0.4e %0.4e %0.4e %0.4e',...
        k,nlinsearch,normdx,normgradL,eq_viol,ineq_viol));
    
    % Check convergence on dx
    if normdx < toldx
        disp('Convergence (small dx)');
        break;
    elseif normgradL < tolgL
        disp('Convergence (small gradL)');
        break
    end

    % Evaluate the functions
    [Jfk,fk] = itrick(ffun,x);
    [Jgk,gk] = itrick(gfun,x);
    [Jhk,hk] = itrick(hfun,x);
    
    % Check if maximum number of iterations reached
    if k >= max_iter
        disp('Maximum number of SQP iterations reached!');
        break;
    end

    % Complete the damped BFGS update (Procedure 18.2 in Nocedal)
    gradL_new = Jfk';
    if ~isempty(gfun) % Equality constraints
        gradL_new = gradL_new - Jgk'*lambda_k;
    end
    if ~isempty(hfun) % Inequality constraints
        gradL_new = gradL_new - Jhk'*mu_k;
    end

    yk = gradL_new - gradL;
    Bdx = Bk*dx;
    dxBdx = dx'*Bdx;
    ydx = dx'*yk;
    if ydx >= 0.2*dxBdx
        thetak = 1;
    else
        thetak = 0.8*dxBdx/(dxBdx - ydx);
    end
    rk = thetak*dx + (1-thetak)*Bdx; % rk replaces yk to assure Bk pos.def.
    Bk = Bk - Bdx*Bdx'/dxBdx + (rk*rk')/ (rk'*dx);

end
disp(sprintf('\nSQP algorithm terminated after %d iterations',k-1));

end % end of sqpmethod


% Line search
function [alpha,mu,iter] = linesearch(ffun,gfun,hfun,x,p,fk,gk,hk,Jfk,Jgk,Jhk,Bk,mu)
% parameters in the algorithm
sigma = 1; % Bk in BDGS is always pos.def.
rho = 0.5; % 
mu_safety = 1.1; % safety factor for mu (see below)
eta = 0.0001; % text to Noc 3.4
tau = 0.2; % 
max_iter = 100;

% 1-norm of the feasability violations
feasviol = sum(abs(gk)) + sum(max(0,-hk));

% Use a quadratic model of T1 to get a lower bound on mu (eq. 18.36
% in Nocedal)
mu_lb = (Jfk*p + sigma/2*p'*Bk*p)/(1-rho)/feasviol;

% Increase mu if it is below the lower bound
if mu < mu_lb
    mu = mu_lb*mu_safety;
end

% Calculate T1 at x (18.27 in Nocedal)
T1 = fk + mu*feasviol;

% Calculate the directional derivative of T1 at x (cf. 18.29 in Nocedal)
DT1 = Jfk*p;
if ~isempty(gfun)
    DT1 = DT1 - mu*sum(abs(gk));
end
if ~isempty(hfun)
    DT1 = DT1  + mu*(hk < 0)'*Jhk*p;
end

iter = 0;
alpha = 1;
while true
    % Evaluate prospective x
    x_new = x+alpha*p;
    fk_new = feval(ffun,x_new);

    % Evaluate gk, hk and get 1-norm of the feasability violations
    feasviol_new = 0;
    if ~isempty(gfun)
        gk_new = feval(gfun,x_new);
        feasviol_new = feasviol_new + sum(abs(gk_new));
    end
    if ~isempty(hfun)
        hk_new = feval(hfun,x_new);
        feasviol_new = feasviol_new + sum(max(0,-hk_new));
    end

    % New T1 function
    T1_new = fk_new + mu*feasviol_new;

    % Check Armijo condition, SQP version (18.28 in Nocedal)
    if T1_new <= T1 + eta*alpha*DT1
        break;
    end

    % Backtrace
    alpha = alpha*tau;

    % Go to next iteration
    iter = iter+1;
    if iter >= max_iter
        error('linesearch failed!');
    end
end

end

% Derivative calculation with the i-trick
function [J,f] = itrick(fun,x)
% This program computes the Jacobian J(x) of a function
% f(x) based on the i-trick.
% written by Joel Andersson, 31/10-2008, joel.andersson@esat.kuleuven.be

% quick return if empty
if isempty(fun)
    J = [];
    f = [];
    return;
end

f = feval(fun,x);

% Discretization step size
e = eps; % Mashine-epsilon

% Dimensions
n = length(x);
m = length(f);

% Jacobian
J = nan(m,n);

% Calculate the entries of the Jacobian
for j=1:n
    % Make a small perturbation in direction j
    x_pert = x;
    x_pert(j) = x_pert(j) + i*e;

    % Calculate the perturbed function value
    f_pert = feval(fun,x_pert);

    % Calculate row j of the Jacobian
    J(:,j) = imag(f_pert) / e;
end

end
