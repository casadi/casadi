import casadi.*

for pl = strsplit(CasadiMeta.getPlugins(),';')
  out  = strsplit(pl{:},'::');
  cls  = out{1};
  name = out{2};
  if strcmp(cls, 'Integrator')
    eval(['Function.load_integrator(''' name ''')'])
  elseif strcmp(cls, 'NlpSolver')
    eval(['Function.load_nlp_solver(''' name ''')'])
  elseif strcmp(cls, 'QpSolver')
    eval(['Function.load_qp_solver(''' name ''')'])
  elseif strcmp(cls, 'ImplicitFunction')
    eval(['Function.load_rootfinder(''' name ''')'])
  else
    eval([cls '.loadPlugin(''' name ''')'])
  end


end

