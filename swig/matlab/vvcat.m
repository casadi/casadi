function r = vvcat(a)
%VVCAT vec-concatenate the contents of a cell array.
%
%  vvcat(a) is equivalent to veccat(a{:}). Returns an empty 0x1
%  native double when `a` is empty (pure-MATLAB-in / pure-MATLAB-out).
%  Wrap with casadi.DM(...) if you need a typed empty.
  if isempty(a)
    r = zeros(0, 1);
  else
    r = veccat(a{:});
  end
end
