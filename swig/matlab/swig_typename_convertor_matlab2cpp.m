function s = swig_typename_convertor_matlab2cpp(a)
    if iscell(a)
       s = ['{' strjoin(unique(cellfun(@swig_typename_convertor_matlab2cpp,a,'UniformOutput',false)),'|') '}'];
    else
       s = strrep(class(a),'casadi.','');
    end
end