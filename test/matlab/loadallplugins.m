import casadi.*

for pl = strsplit(CasadiMeta.getPlugins(),';')
  out  = strsplit(pl{:},'::');
  cls  = out{1};
  name = out{2};
  if strcmp(cls, 'Integrator')
    eval(['Function.load_integrator(''' name ''')'])
  elseif strcmp(cls, 'NlpSolver')
    eval(['Function.load_nlpsol(''' name ''')'])
  elseif strcmp(cls, 'Qpsol')
    eval(['Function.load_qpsol(''' name ''')'])
  elseif strcmp(cls, 'Rootfinder')
    eval(['Function.load_nlsol(''' name ''')'])
  elseif strcmp(cls, 'LinearSolver')
    eval(['Function.load_linsol(''' name ''')'])
  else
    eval([cls '.loadPlugin(''' name ''')'])
  end


end

