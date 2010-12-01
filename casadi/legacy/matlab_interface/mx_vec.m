classdef mx_vec < handle
    properties (SetAccess = protected)
        ind = 0;
    end
    
    methods
        % constructor
        function obj = mx_vec(varargin)
            obj.ind = calllib('libcasadi','casadi_mx_vec_new');
            assert(obj.ind ~=0);
            
            for i=1:nargin
                push_back(obj,varargin{i});
            end
            
        end
        
        % destructor
        function delete(h)
            ret = calllib('libcasadi','casadi_mx_vec_delete',h.ind);
            assert(ret==0);
        end
        
        % add an element to the back
        function push_back(h,el)
            % recursive call if a cell array
            if iscell(el)
                for i=1:length(el)
                    push_back(h,el{i});
                end
                return;
            end
            
            % convert to mx if necessary
            if ~isa(el,'mx')
                el = mx(el);
            end
            
            % add to array
            ret = calllib('libcasadi','casadi_mx_vec_push_back',h.ind, el.ind);
            assert(ret==0);
        end

        % get the length
        function sz = length(h)
            sz = calllib('libcasadi','casadi_mx_vec_size',h.ind);
        end
        
        % display
        function disp(h)
	      flag = calllib('libcasadi','casadi_mx_vec_print_cout',h.ind);
              assert(flag==0);
        end

    end

end
