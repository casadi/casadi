function r = dcat(a)
%DCAT Diagonally concatenate the contents of a cell array.
%
%  dcat(a) is equivalent to diagcat(a{:}). Returns an empty 0x0
%  native double when `a` is empty (pure-MATLAB-in / pure-MATLAB-out).
%  Wrap with casadi.DM(...) if you need a typed empty.
  if isempty(a)
    r = zeros(0, 0);
  else
    r = diagcat(a{:});
  end
end
