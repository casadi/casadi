classdef fx < handle
    properties (SetAccess = protected)
        ind = -1;
    end
    properties (Constant)        
        % id of the class
        FXID = 1;
        
        % ids of the methods
        FXDELETE = 0;
        FXDISP = 1;
        FXGETOPTION = 2;
        FXSETOPTION = 3;
        FXGETINPUT = 4;
        FXSETINPUT = 5;
        FXGETOUTPUT = 6;
        FXSETOUTPUT = 7;
        FXEVALUATE = 8;
        FXEVALUATEFWD = 9;
        FXEVALUATEADJ = 10;
    end
    methods
        % constructor
        function obj = fx()
	    obj.ind = calllib('libcasadi','casadi_fx_new');
        end
        % destructor
        function delete(h)
	    flag = calllib('libcasadi','casadi_fx_delete',h.ind);
	    assert(flag == 0);
        end
        % display
        function disp(h)
	      flag = calllib('libcasadi','casadi_fx_print_cout',h.ind);
              assert(flag==0);
        end
        % get options
        function res = get_option(h,opname)
	    error('not implemented');
            % res = casadi_matlab(fx.FXID,fx.FXGETOPTION,h.ind,opname);
%            casadi_matlab(fx.FXID,fx.FXGETOPTION,h.ind,opname);
        end
        % set options
        function set_option(h,opname,opval)
	    assert(ischar(opname));
	    if ischar(opval)
	      flag = calllib('libcasadi','casadi_fx_setoption_string',h.ind,opname,opval);
	      assert(flag==0);
        else
            assert(isnumeric(opval));
            flag = calllib('libcasadi','casadi_fx_setoption_double',h.ind,opname,opval,length(opval));
            assert(flag==0);
	    end;
           assert(nargin==3);
        end
        
        % get input (seed)
        function val = get_input(h,iind,ord)
            if nargin==1
                iind = 0;
                ord = 0;
            elseif nargin==2;
                ord = 0;
            else
                assert(nargin==3);
            end
	    sz = libpointer('int32Ptr', 0);
	    flag = calllib('libcasadi','casadi_fx_input_size',h.ind,iind,sz);
	    assert(flag==0);

	    v = libpointer('doublePtr', nan(sz.Value,1));
	    flag = calllib('libcasadi','casadi_fx_getinput',h.ind,iind,ord,v);
	    assert(flag==0);

	    val = v.Value;
        end
        
        % set input (seed)
        function h = set_input(h,val,iind,ord)
            if nargin==2
                iind = 0;
                ord = 0;
            elseif nargin==3;
                ord = 0;
            else
                assert(nargin==4);
            end
	    flag = calllib('libcasadi','casadi_fx_setinput',h.ind,iind,ord,val);
	    assert(flag==0);
        end
         
        % get output (seed)
        function val = get_output(h,oind,ord)
            if nargin==1
                oind = 0;
                ord = 0;
            elseif nargin==2;
                ord = 0;
            else
                assert(nargin==3);
            end
	    sz = libpointer('int32Ptr', 0);
	    flag = calllib('libcasadi','casadi_fx_output_size',h.ind,oind,sz);
	    assert(flag==0);

	    v = libpointer('doublePtr', nan(sz.Value,1));
	    flag = calllib('libcasadi','casadi_fx_getoutput',h.ind,oind,ord,v);
	    assert(flag==0);

	    val = v.Value;
        end
        
        % set output (seed)
        function h = set_output(h,val,oind,ord)
            if nargin==2
                oind = 0;
                ord = 0;
            elseif nargin==3;
                ord = 0;
            else
                assert(nargin==4);
            end
	    flag = calllib('libcasadi','casadi_fx_setoutput',h.ind,oind,ord,val);
	    assert(flag==0);
        end
        
        % initialize
        function h = init(h)
	      flag = calllib('libcasadi','casadi_fx_init',h.ind);
          assert(flag==0);
        end
        
        % evaluate
        function h = evaluate(h)
	    flag = calllib('libcasadi','casadi_fx_evaluate',h.ind);
	    assert(flag==0);
        end
        
        % evaluate adjoint
        function h = evaluate_adj(h)
	    flag = calllib('libcasadi','casadi_fx_evaluate_adj',h.ind);
	    assert(flag==0);
        end
        
        % evaluate forward
        function h = evaluate_fwd(h)
	    flag = calllib('libcasadi','casadi_fx_evaluate_fwd',h.ind);
	    assert(flag==0);
        end
        
        % jacobian of the expression
        function res = jacobian(obj)
            % return function
            res = fx();

            % evaluate
            iind = 0;
            oind = 0;
            flag = calllib('libcasadi','casadi_fx_jacobian',obj.ind, res.ind, iind, oind);
            assert(flag==0);
        end
        
        % hessian of the expression
        function res = hessian(obj)
            % return function
            res = fx();

            % evaluate
            iind = 0;
            oind = 0;
            flag = calllib('libcasadi','casadi_fx_hessian',obj.ind, res.ind, iind, oind);
            assert(flag==0);
        end
        
    end
end