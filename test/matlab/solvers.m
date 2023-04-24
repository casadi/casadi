import casadi.*

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

assert(~failure)
