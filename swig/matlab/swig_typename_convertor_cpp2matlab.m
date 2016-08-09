function s = swig_typename_convertor_cpp2matlab(s)
  s = strrep(s,'C/C++ prototypes','Matlab usages');
  s = strrep(s,'casadi::','');
  s = strrep(s,'Dict','struct');
  s = strrep(s,'MXDict','struct:MX');
  s = strrep(s,'SXDict','struct:SX');
  s = strrep(s,'std::string','char');
  s = strrep(s,' const &','');
  s = strrep(s,'friendwrap_','');
  s = regexprep(s,'\<(\w+)(< \w+ >)::\1\>','$1$2');
  s = regexprep(s,'\<(\w+)::\1\>','$1');  
  s = regexprep(s,'(const )?Matrix< ?SXElem *>( &)?','SX');
  s = regexprep(s,'(const )?Matrix< ?double *>( &)?','DM');
  s = regexprep(s,'(const )?Matrix< ?int * >( &)?','IM');
  s = regexprep(s,'(const )?GenericMatrix< ?(\w+) *>( ?&)?','$2');
  s = regexprep(s,'(const )?GenericMatrix< ?([\w\(\)]+) *>( ?&)?','$2 ');
  s = regexprep(s,'const (\w+) &','$1 ');
  s = regexprep(s,'< [\w\(\)]+ +>\(','(');

  for i=1:5
      s = regexprep(s,'(const )? ?std::pair< ?([\w\(\)\]\[\}\{: ]+?) ?, ?([\w\(\)\]\[\}\{: ]+?) ?> ?&?','{$2,$3} ');
      s = regexprep(s,'(const )? ?std::vector< ?([\w\(\)\[\]\}\{ ]+) ?(, ?std::allocator< ?\2 ?>)? ?> ?&?','{$2} ');
      s = regexprep(s,'(const )? ?std::map< ?([\w\(\)\[\]\}\{ ]+) ?, ?() ?(, ?std::allocator< ?\2 ?>)? ?> ?&?','{$2} ');
      s = regexprep(s,'\{int\}','[int]');
      s = regexprep(s,'\{double\}','[double]');
  end
  s = regexprep(s,'\<(\w+)(< \w+ >)?::\1','$1');

  s = strrep(s,'casadi::','');
  s = strrep(s,'IOInterface< Function >','Function');
  s = strrep(s,'::','.');
  s = strrep(s,'.operator ()','');
  s = regexprep(s,'([A-Z]\w+)Vector','{$1}');
end
