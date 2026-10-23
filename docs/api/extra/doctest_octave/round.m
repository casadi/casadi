function y = round(x, n)
  % Octave < 9 lacks round(x, n); minimal version for the CasADi doctests.
  if nargin < 2
    y = builtin('round', x);
  else
    y = builtin('round', x .* 10.^n) ./ 10.^n;
  end
end
