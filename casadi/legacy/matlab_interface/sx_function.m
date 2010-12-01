classdef sx_function < fx

  methods

  % constructor
    function obj = sx_function(var,res)
        if ~isa(var,'sx_matrix_vec')
            var = sx_matrix_vec(var);
        end
        if ~isa(res,'sx_matrix_vec')
            res = sx_matrix_vec(res);
        end
        flag = calllib('libcasadi','casadi_sx_function',obj.ind, var.ind, res.ind);
        assert(flag==0);
    end
    
    function res = evaluate_symbolic(obj,arg)
        if ~isa(arg,'sx_matrix_vec')
            arg = sx_matrix_vec(arg);
        end
        
        % return matrix
        res = sx_matrix();
        
        % evaluate
        flag = calllib('libcasadi','casadi_sx_function_evaluate_symbolically',obj.ind, arg.ind, res.ind);
        assert(flag==0);
    end
    
  end % methods
  
  

end % sx_function
