from casadi import *
from numpy import *
import matplotlib.pyplot as plt

# The algorithm is a Quasi-Newton method with damped BFGS updating to that
# assures positive definitenes of the Hessian approximation. Line search is
# carried out via backtracking until with the Armijo condition applied to
# the T1 (in Nocedal phi1) merit function is satisfied.

# Controls
u = symbolic("u",10)

# Starting position
s1=0
v1=0
m1=1

# Final position
s100=SX(10)
v100=SX(0)

# Time step
dt=SX(0.1)
alpha=SX(0.05)
beta=SX(0.1)

# Constraints
vmax=1.3

# starting point
u0=0.4*ones(10)
umin = -1*ones(10)
umax = 0.5*ones(10)

# minimize fuel
f = inner_prod(u,u)/2

# Calculate acceleration for all points
a = SXMatrix(100,1,0)
for i in range(10):
  a[10*i:10*(i+1)] = u[i]

# Loop over all k to get s and v and the endpoints
s = SXMatrix(101,1,0)
s[0] = s1
v = SXMatrix(101,1,0)
v[0] = v1
m = SXMatrix(101,1,0)
m[0] = m1

for k in range(100):
  s[k+1] = s[k] + dt*v[k]
  v[k+1] = v[k] + dt / m[k] * (a[k] - alpha * v[k]**2)
  m[k+1] = m[k] - dt * beta * a[k]**2

# Equality constraint function
g = [s[100]-s100, v[100]-v100]

# Create NLP functions
ffcn = SXFunction([u],[f])
ffcn.init()
gfcn = SXFunction([u],[g])
gfcn.init()

# Solve with IPOPT
solver = IpoptSolver(ffcn,gfcn)
solver.setOption("hessian_approximation", "limited-memory")

#solver.setOption("verbose",True)
solver.init()

# Set bounds and initial guess
solver.setInput(umin,     NLP_LBX)
solver.setInput(umax,     NLP_UBX)
solver.setInput(u0,       NLP_X_INIT)
solver.setInput(zeros(2), NLP_LBG)
solver.setInput(zeros(2), NLP_UBG)

# Solve the problem
solver.solve()

# Retrieve the solution
u_opt = array(solver.output(NLP_X_OPT))

# Get values at the beginning of each finite element
tgrid = linspace(0,10,10)

# Plot the results
plt.figure(1)
plt.clf()
#plt.plot(tgrid,x_opt,'--')
#plt.plot(tgrid,y_opt,'-')
plt.plot(tgrid,u_opt,'-.')
plt.title("Rocket trajectory optimization")
plt.xlabel('time')
#plt.legend(['x trajectory','y trajectory','u trajectory'])
plt.grid()

# Initial guess
x = u0

# Parameters in the algorithm
maxiter = 100 # maximum number of sqp iterations
toldx = 1e-12 # stopping criterion for the stepsize
tolgL = 1e-12 # stopping criterion for the lagrangian gradient
merit_mu = 0  # current 'mu' in the T1 merit function

# Get dimensions
m = gfcn.output().size() # Number of equality constraints
n = len(x)  # Number of variables
q = 0 # Number of inequality constraints

# Initial guess for the lagrange multipliers
lambda_k = zeros(m)
mu_k = zeros(q)

# Initial guess for the Hessian
Bk = eye(n)

# Jacobian function
jfcn = gfcn.jacobian()
jfcn.init()

# Allocate a QP solver
H_sparsity = sp_dense(n,n)
G_sparsity = sp_dense(n)
A_sparsity = jfcn.output().sparsity()
qp_solver = QPOasesSolver(H_sparsity,G_sparsity,A_sparsity)
#qp_solver = IpoptQPSolver(H_sparsity,G_sparsity,A_sparsity)
qp_solver.init()

# No bounds on the control
qp_solver.input(QP_LBX).setAll(-inf)
qp_solver.input(QP_UBX).setAll( inf)

# Header
print ' k  nls | dx         gradL      eq viol    ineq viol'
k = 0

# Evaluate the functions
#[Jfk,fk] = itrick(ffun,x);
#[Jgk,gk] = itrick(gfun,x);

while True:
  # Evaluate the constraint function
  gfcn.setInput(x)
  gfcn.evaluate()
  
  # Evaluate the Jacobian
  jfcn.setInput(x)
  jfcn.evaluate()
  
  # Evaluate the gradient of the objective function
  ffcn.setInput(x)
  ffcn.setAdjSeed(1.0)
  ffcn.evaluate(0,1)
  
  # Pass data to QP solver
  qp_solver.setInput(Bk,QP_H)
  qp_solver.setInput(jfcn.output(),QP_A)
  qp_solver.setInput(ffcn.adjSens(),QP_G)
  qp_solver.setInput(-gfcn.output(),QP_LBA)
  qp_solver.setInput(-gfcn.output(),QP_UBA)

  # Solve the QP subproblem
  qp_solver.evaluate()
  
  # Get the optimal solution
  p = qp_solver.output(QP_PRIMAL)
  
  # Get the dual solution
  




  break
  
    #% Get the gradient of the Lagrangian
          ##gradL = Jfk';
          #if ~isempty(gfun) gradL = gradL - Jgk'*lambda_hat; end
          #if ~isempty(hfun) gradL = gradL - Jhk'*mu_hat;     end
          
          #% Do a line search along p
          #[tk,merit_mu,nlinsearch] = linesearch(ffun,gfun,hfun,x,...
              #p,fk,gk,hk,Jfk,Jgk,Jhk,Bk,merit_mu);

          #% Calculate the new step
          #dx = p*tk;
          #x = x + dx;
          #if ~isempty(gfun) lambda_k = tk*lambda_hat + (1-tk)*lambda_k; end
          #if ~isempty(hfun) mu_k     = tk*mu_hat + (1-tk)*mu_k;         end
          #k = k+1;

          #% Gather and print iteration information
          #normdx = norm(dx); % step size
          #normgradL = norm(gradL); % size of the Lagrangian gradient
          #eq_viol = sum(abs(gk)); % equality constraint violation
          #ineq_viol = sum(max(0,-hk)); % inequality constraint violation

          #disp(sprintf('%3d %3d |%0.4e %0.4e %0.4e %0.4e',...
              #k,nlinsearch,normdx,normgradL,eq_viol,ineq_viol));

          #% Check convergence on dx
          #if normdx < toldx
              #disp('Convergence (small dx)');
              #break;
          #elseif normgradL < tolgL
              #disp('Convergence (small gradL)');
              #break
          #end

          #% Evaluate the functions
          #[Jfk,fk] = itrick(ffun,x);
          #[Jgk,gk] = itrick(gfun,x);
          #[Jhk,hk] = itrick(hfun,x);
          
          #% Check if maximum number of iterations reached
          #if k >= maxiter
              #disp('Maximum number of SQP iterations reached!');
              #break;
          #end

          #% Complete the damped BFGS update (Procedure 18.2 in Nocedal)
          #gradL_new = Jfk';
          #if ~isempty(gfun) % Equality constraints
              #gradL_new = gradL_new - Jgk'*lambda_k;
          #end
          #if ~isempty(hfun) % Inequality constraints
              #gradL_new = gradL_new - Jhk'*mu_k;
          #end
          
          #yk = gradL_new - gradL;
          #Bdx = Bk*dx;
          #dxBdx = dx'*Bdx;
          #ydx = dx'*yk;
          #if ydx >= 0.2*dxBdx
              #thetak = 1;
          #else
              #thetak = 0.8*dxBdx/(dxBdx - ydx);
          #end
          #rk = thetak*dx + (1-thetak)*Bdx; % rk replaces yk to assure Bk pos.def.
          #Bk = Bk - Bdx*Bdx'/dxBdx + (rk*rk')/ (rk'*dx);

      #end
      #disp(sprintf('\nSQP algorithm terminated after %d iterations',k-1));

      #end % end of sqpmethod


      #% Line search
      #function [alpha,mu,iter] = linesearch(ffun,gfun,hfun,x,p,fk,gk,hk,Jfk,Jgk,Jhk,Bk,mu)
      #% parameters in the algorithm
      #sigma = 1; % Bk in BDGS is always pos.def.
      #rho = 0.5; % 
      #mu_safety = 1.1; % safety factor for mu (see below)
      #eta = 0.0001; % text to Noc 3.4
      #tau = 0.2; % 
      #maxiter = 100;

      #% 1-norm of the feasability violations
      #feasviol = sum(abs(gk)) + sum(max(0,-hk));

      #% Use a quadratic model of T1 to get a lower bound on mu (eq. 18.36
      #% in Nocedal)
      #mu_lb = (Jfk*p + sigma/2*p'*Bk*p)/(1-rho)/feasviol;

      #% Increase mu if it is below the lower bound
      #if mu < mu_lb
          #mu = mu_lb*mu_safety;
      #end

      #% Calculate T1 at x (18.27 in Nocedal)
      #T1 = fk + mu*feasviol;

      #% Calculate the directional derivative of T1 at x (cf. 18.29 in Nocedal)
      #DT1 = Jfk*p;
      #if ~isempty(gfun)
          #DT1 = DT1 - mu*sum(abs(gk));
      #end
      #if ~isempty(hfun)
          #DT1 = DT1  + mu*(hk < 0)'*Jhk*p;
      #end
          
      #iter = 0;
      #alpha = 1;
      #while true
          #% Evaluate prospective x
          #x_new = x+alpha*p;
          #fk_new = feval(ffun,x_new);

          #% Evaluate gk, hk and get 1-norm of the feasability violations
          #feasviol_new = 0;
          #if ~isempty(gfun)
              #gk_new = feval(gfun,x_new);
              #feasviol_new = feasviol_new + sum(abs(gk_new));
          #end
          #if ~isempty(hfun)
              #hk_new = feval(hfun,x_new);
              #feasviol_new = feasviol_new + sum(max(0,-hk_new));
          #end

          #% New T1 function
          #T1_new = fk_new + mu*feasviol_new;

          #% Check Armijo condition, SQP version (18.28 in Nocedal)
          #if T1_new <= T1 + eta*alpha*DT1
              #break;
          #end

          #% Backtrace
          #alpha = alpha*tau;

          #% Go to next iteration
          #iter = iter+1;
          #if iter >= maxiter
              #error('linesearch failed!');
          #end
      #end

      #end

      #% Derivative calculation with the i-trick
      #function [J,f] = itrick(fun,x)
      #% This program computes the Jacobian J(x) of a function
      #% f(x) based on the i-trick.
      #% written by Joel Andersson, 31/10-2008, joel.andersson@esat.kuleuven.be

      #% quick return if empty
      #if isempty(fun)
          #J = [];
          #f = [];
          #return;
      #end

      #f = feval(fun,x);

      #% Discretization step size
      #e = eps; % Mashine-epsilon

      #% Dimensions
      #n = length(x);
      #m = length(f);

      #% Jacobian
      #J = nan(m,n);

      #% Calculate the entries of the Jacobian
      #for j=1:n
          #% Make a small perturbation in direction j
          #x_pert = x;
          #x_pert(j) = x_pert(j) + i*e;

          #% Calculate the perturbed function value
          #f_pert = feval(fun,x_pert);

          #% Calculate row j of the Jacobian
          #J(:,j) = imag(f_pert) / e;
      #end

      #end


























plt.show()










