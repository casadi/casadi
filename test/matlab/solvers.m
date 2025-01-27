import casadi.*

for pl = strsplit(CasadiMeta.plugins(),';')
  out  = strsplit(pl{:},'::');
  cls  = out{1};
  name = out{2};

if ~isempty(getenv(['SKIP_' upper(name) '_TESTS']))
  continue
end

if strcmp(cls,'XmlFile')
  % do nothing
elseif strcmp(cls,'Importer') || strcmp(cls,'Linsol')
  eval([cls '.load_plugin(''' name ''')'])
else
  eval(['load_' lower(cls) '(''' name ''')'])
end

end

plugins = strsplit(CasadiMeta.plugins(),';');

failure = false;
for i=1:length(plugins)
  plugin = plugins{i};
  parts = strsplit(plugin,'::');
  type = parts{1};
  name = parts{2};
  if ~isempty(getenv(['SKIP_' upper(name) '_TESTS']))
      continue
  end
  
  if strcmp(type,'Nlpsol')
      if ismember(name, {'scpgen', 'ampl', 'qrsqp', 'fatrop'})
         continue 
      end
      disp(['Solver:' name])
      
      x = MX.sym('x');
      nlp = struct;
      nlp.x = x;
      nlp.f = x;
      nlp.g = 2*x+1;
      
      try
        solver = nlpsol('solver',name, nlp);
        solver('lbg',0,'ubg',2);
      catch E
        disp(['Failure for nlpsol ' name])
        disp(E.message)
        failure = true;
      end
      
  end
  
  if strcmp(type,'Conic')
      if ismember(name, {'hpipm', 'nlpsol', 'qrqp', 'ipqp', 'fatrop', 'daqp'})
         continue 
      end
      disp(['Solver:' name])
      
      x = MX.sym('x');
      nlp = struct;
      nlp.x = x;
      nlp.f = x;
      nlp.g = 2*x+1;
      
      try
        solver = qpsol('solver',name, nlp);
        solver('lbg',0,'ubg',2);
      catch E
        disp(['Failure for qpsol ' name])
        disp(E.message)
        failure = true;
      end
      
  end
end

if ismember('Conic::hpipm',plugins)
try
    A=[1, 0.2, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0;
     -0.1, 0.4, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
     0.3, 0.2, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
     2, 0, 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     1, 1, 0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
     0, 0, 0, 1, 4, 2, 1, 0.3, -1, 0, 0, 0; 
     0, 0, 0, 3, 1, 0, 1, 0.2, 0, -1, 0, 0; 
     0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 2, 4, 0, -1; 
     0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 1, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3];

    A = DM(sparse(A));

    H=[7, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
     0, 7, 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
     0.2, 0.3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
     0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0; 
     0, 0, 0, 0, 2, 0.1, 0, 0.7, 0, 0, 0, 0; 
     0, 0, 0, 0, 0.1, 1, 0, 1, 0, 0, 0, 0; 
     0, 0, 0, 0, 0, 0, 1, 0.1, 0, 0, 0, 0; 
     0, 0, 0, 1, 0.7, 1, 0.1, 2, 0, 0, 0, 0; 
     0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 1, 0; 
     0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0; 
     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 4, 0; 
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9];

    H = DM(sparse(H));

    nx = [2,3,2,1];
    nu = [1, 2,1];
    ng = [2, 1, 1, 1];
    N = 3;


    options = struct;
    options.hpipm.iter_max = 100;
    options.hpipm.res_g_max = 1e-10;
    options.hpipm.res_b_max = 1e-10;
    options.hpipm.res_d_max = 1e-10;
    options.hpipm.res_m_max = 1e-10;

    solver = conic('solver', 'hpipm', struct('a',A.sparsity(), 'h', H.sparsity()),options);

    g = [1;1;0.2;0.4;1;0.5;0.3;1;0.6;1;1;0.7];
    lbg =[0;0;0;-2;-2;0;0;-2;0;-2;-2];
    ubg = [0;0;0;2;2;0;0;2;0;2;2];
    lbx = [0.5;0.2;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1];
    ubx = [0.5;0.2; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];

    sol = solver('a',A,'h',H,'lba',lbg,'uba',ubg,'g',g,'lbx',lbx,'ubx',ubx);
catch E
    disp(['Failure for hpipm'])
    disp(E.message)
    failure = true;
end
end

if ismember('Nlpsol::fatrop',plugins)
try

% Time horizon and number of control intervals
T = 10; 
N = 10;

% Declare model variables
import casadi.*
x1 = MX.sym('x1');
x2 = MX.sym('x2');
x = [x1; x2];
u = MX.sym('u');
p = MX.sym('p');

% Model equations
xdot = [(1-x2^2)*x1 - x2 + u + p, x1];

% Define integrator
opts = struct('simplify', true, 'number_of_finite_elements', 1);
F = casadi.integrator('F', 'rk', struct('x', x, 'p', p, 'u', u, 'ode', xdot), opts);

% Start with an empty NLP
w = {};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g = {};
lbg = [];
ubg = [];

% "Lift" initial conditions
Xk = MX.sym('X0', 2);
w = {w{:}, Xk};
lbw = [lbw, 0, 1];
ubw = [ubw, 0, 1];
w0 = [w0, 0.1, 0.2];

% Formulate the NLP
for k = 1:N
    % New NLP variable for the control
    Uk = MX.sym(['U_', num2str(k)]);
    w = {w{:}, Uk};
    lbw = [lbw, -1];
    ubw = [ubw, 1];
    w0 = [w0, 0.3];

    % Integrate till the end of the interval
    Fk = F('x0', Xk, 'u', Uk, 'p', p);
    Xk_end = Fk.xf;
    J = J + sumsqr(Xk) + sumsqr(Uk);

    % New NLP variable for state at end of interval
    Xk_next = MX.sym(['X_', num2str(k+1)], 2);
    w = {w{:}, Xk_next};
    lbw = [lbw, -0.25, -inf];
    ubw = [ubw, inf, inf];
    w0 = [w0, 0.1, 0.2];
        
    % Add equality constraint
    g = {g{:}, Xk_next-Xk_end};
    lbg = [lbg, 0, 0];
    ubg = [ubg, 0, 0];

    Xk = Xk_next;
end

nlp = struct;
nlp.f = J;
nlp.x = vertcat(w{:});
nlp.g = vertcat(g{:});
nlp.p = p;

args = struct;
args.x0 = w0;
args.lbx = lbw;
args.ubx = ubw;
args.lbg = lbg;
args.ubg = ubg;
args.p = 0;

f = nlpsol('solver', 'fatrop', nlp);
f.call(args);

catch E
    disp(['Failure for fatrop'])
    disp(E.message)
    failure = true;
end
end


assert(~failure)
disp('Search the log for Failure')
