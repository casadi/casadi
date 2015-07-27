function s = swig_typename_convertor_cpp2matlab(s)
  s = strrep(s,'C/C++ prototypes','Matlab usages');
  s = strrep(s,'casadi::','');
  s = strrep(s,'MXDict','struct:MX');
  s = strrep(s,'SXDict','struct:SX');
  s = strrep(s,'std::string','char');
  s = strrep(s,' const &','');
  
  s = regexprep(s,'(const )?Matrix< SXElement >( &)?','SX');
  s = regexprep(s,'(const )?GenericMatrix< ?(\w+) *>( ?&)?','$2');
  s = regexprep(s,'(const )?Matrix< ?(\w+) *>( ?&)?','array($2) ');
  s = regexprep(s,'(const )?GenericMatrix< ?([\w\(\)]+) *>( ?&)?','$2 ');
  s = regexprep(s,'const (\w+) &','$1 ');
  s = regexprep(s,'< [\w\(\)]+ +>\(','(');
  s = regexprep(s,'\<(\w+)(< \w+ >)?::\1','$1');
  for i=1:5
      s = regexprep(s,'(const )? ?std::pair< ?([\w\(\)\]\[\}\{: ]+?) ?, ?([\w\(\)\]\[\}\{: ]+?) ?> ?&?','{$2,$3} ');
      s = regexprep(s,'(const )? ?std::vector< ?([\w\(\)\[\]\}\{ ]+) ?(, ?std::allocator< ?\2 ?>)? ?> ?&?','{$2} ');
  end
  s = regexprep(s,'\<(\w+)(< \w+ >)?::\1','$1');

  s = strrep(s,'casadi::','');
  s = strrep(s,'IOInterface< Function >','Function');
  s = strrep(s,'::','.');
  s = strrep(s,'.operator ()','');
  s = regexprep(s,'([A-Z]\w+)Vector','{$1}');
end
