classdef mx_function < fx

  methods

  % constructor
    function obj = mx_function(var,res)
        if ~isa(var,'mx_vec')
            var = mx_vec(var);
        end
        if ~isa(res,'mx_vec')
            res = mx_vec(res);
        end
        flag = calllib('libcasadi','casadi_mx_function',obj.ind, var.ind, res.ind);
        assert(flag==0);
    end

  end % methods

end % mx_function
