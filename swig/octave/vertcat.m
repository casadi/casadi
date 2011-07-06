function y = vertcat( varargin )
    disp('vertcat called on swig')
    global casadi;
    y = casadi.__vertcat__(varargin);
end
