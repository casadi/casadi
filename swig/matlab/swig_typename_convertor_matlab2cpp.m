function s = swig_typename_convertor_matlab2cpp(a)
    if iscell(a)
       s = ['{' strjoin(unique(cellfun(@swig_typename_convertor_matlab2cpp,a,'UniformOutput',false)),'|') '}'];
       s = regexprep(s,'\{int\}','[int]');
       s = regexprep(s,'\{double\}','[double]');
    elseif isstruct(a)
       names = fieldnames(a);
       values = cell(1,numel(names));
       for i=1:numel(names)
           values{i} = a.(names{i});
       end
       s = ['struct:' strjoin(unique(cellfun(@swig_typename_convertor_matlab2cpp,values,'UniformOutput',false)),'|')];
    else
       s = strrep(class(a),'casadi.','');
    end
end

