function tf = contains(str, pattern)
  % Octave < 9 lacks contains(); minimal version for the CasADi doctests.
  tf = ~isempty(strfind(str, pattern));
end
