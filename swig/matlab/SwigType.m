function s = SwigType(varargin)
  s = strjoin(cellfun(@swig_typename_convertor_matlab2cpp,varargin,'UniformOutput',false),', ');
end
