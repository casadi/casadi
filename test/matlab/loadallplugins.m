import casadi.*

for pl = strsplit(CasadiMeta.getPlugins(),';')
  out  = strsplit(pl{:},'::');
  cls  = out{1};
  name = out{2};

if cls=='Importer':
  eval([cls '.load_plugin(''' name ''')'])
else
  eval(['load_' lower(cls) '(''' name ''')'])
end
    
end
