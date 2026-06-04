function r = hcat(a)
%HCAT Horizontally concatenate the contents of a cell array.
%
%  hcat(a) is equivalent to horzcat(a{:}). Returns an empty 1x0
%  native double when `a` is empty (pure-MATLAB-in / pure-MATLAB-out).
%  Wrap with casadi.DM(...) if you need a typed empty.
  if isempty(a)
    r = zeros(1, 0);
  else
    r = horzcat(a{:});
  end
end
