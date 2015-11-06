import casadi.*

for pl = strsplit(CasadiMeta.getPlugins(),';')
  out  = strsplit(pl{:},'::');
  cls  = out{1};
  name = out{2};
  if strcmp(cls, 'Ivpsol')
    eval(['Function.load_ivpsol(''' name ''')'])
  elseif strcmp(cls, 'Nlpsol')
    eval(['Function.load_nlpsol(''' name ''')'])
  elseif strcmp(cls, 'Qpsol')
    eval(['Function.load_qpsol(''' name ''')'])
  elseif strcmp(cls, 'Nlsol')
    eval(['Function.load_nlsol(''' name ''')'])
  elseif strcmp(cls, 'Linsol')
    eval(['Function.load_linsol(''' name ''')'])
  else
    eval([cls '.loadPlugin(''' name ''')'])
  end


end

