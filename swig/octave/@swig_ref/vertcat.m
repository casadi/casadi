function y = vertcat( varargin )
    global casadi;
    y = casadi.__vertcat__(varargin);
end
