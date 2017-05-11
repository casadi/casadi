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
