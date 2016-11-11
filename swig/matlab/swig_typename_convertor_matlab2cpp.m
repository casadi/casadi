function s = swig_typename_convertor_matlab2cpp(a)
    if iscell(a)
       s = ['{' strjoin(unique(cellfun(@swig_typename_convertor_matlab2cpp,a,'UniformOutput',false)),'|') '}'];
       s = regexprep(s,'\{int\}','[int]');
       s = regexprep(s,'\{double\}','[double]');
    else
       s = strrep(class(a),'casadi.','');
    end
end
