function y = horzcat( varargin )
    disp('vertcat called on swig')
    global casadi;
    y = casadi.__horzcat__(varargin);
end
