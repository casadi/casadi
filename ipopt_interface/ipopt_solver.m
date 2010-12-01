classdef ipopt_solver < fx

  methods

  % constructor
    function obj = ipopt_solver(ffcn,gfcn,hfcn,jfcn)
        % All functions
        f = {obj,ffcn,gfcn,hfcn,jfcn};
        
        % List of indices
        f_ind = int32(zeros(5,1));
        
        % Get the pointers
        for i=1:5
            sz = libpointer('int32Ptr', 0);
            calllib('libcasadi','casadi_ptr_to_int',f{i}.ind,sz);
            f_ind(i) = sz.Value;
        end
        
        % Create a solver
        ipopt_solver_mex(f_ind(1),f_ind(2),f_ind(3),f_ind(4),f_ind(5));
    end
        
  end % methods
  
  
end % sx_function
