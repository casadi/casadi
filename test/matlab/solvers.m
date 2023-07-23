import casadi.*

for pl = strsplit(CasadiMeta.plugins(),';')
  out  = strsplit(pl{:},'::');
  cls  = out{1};
  name = out{2};

if strcmp(cls,'Importer') || strcmp(cls,'XmlFile') || strcmp(cls,'Linsol')
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
      if ismember(name, {'scpgen', 'ampl', 'qrsqp'})
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
        disp(getReport(E))
        failure = true;
      end
      
  end
  
  if strcmp(type,'Conic')
      if ismember(name, {'hpipm', 'nlpsol', 'qrqp', 'ipqp'})
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
        disp(getReport(E))
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
    disp(getReport(E))
    failure = true;
end
end

assert(~failure)
