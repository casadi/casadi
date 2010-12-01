classdef sx_matrix < handle
    properties (SetAccess = protected)
        ind = -1;
    end
    properties (Constant)        
        % elementary functions
        ADD_NODE = 0;
        SUB_NODE = 1;
        MUL_NODE = 2;
        DIV_NODE = 3;
        NEG_NODE = 4;
        EXP_NODE = 5;
        LOG_NODE = 6;
        POW_NODE = 7;
        SQRT_NODE = 8;
        SIN_NODE = 9;
        COS_NODE = 10;
        TAN_NODE = 11;
        ASIN_NODE = 12;
        ACOS_NODE = 13;
        ATAN_NODE = 14;
        STEP_NODE = 15;
        FLOOR_NODE = 16;
        CEIL_NODE = 17;
        EQUALITY_NODE = 18;
        ERF_NODE = 19;
        FMIN_NODE = 20;
        FMAX_NODE = 21;
    end
    methods
        % constructor
        function obj = sx_matrix(varargin)
	    obj.ind = calllib('libcasadi','casadi_sx_matrix_new');
	    if (nargin==0)
            return;
        end
        
        assert(nargin>=1);
        % symbolic variable if first argument is a string
        if(isa(varargin{1},'char'))
            % dimension
            switch(nargin)
                case 1
                    n=1; m=1;
                 case 2
                    n=varargin{2}; m=1;
                 case 3
                    n=varargin{2}; m=varargin{3};
                otherwise
                    error('too many arguments');
            end
            ret = calllib('libcasadi','casadi_sx_matrix_symbol',obj.ind,varargin{1},n,m);
            assert(ret==0);
        elseif (nargin==1)
                % argument is a constant
                data = varargin{1};
                assert(isnumeric(data));
                [nrow,ncol] = size(data);
                order = double('C');
                ret = calllib('libcasadi','casadi_sx_matrix_constant',obj.ind,data,nrow,ncol,order);
                assert(ret==0);
        else
                assert(nargin==2);
                % empty matrix
                nrow=varargin{1}; 
                ncol=varargin{2};
                data = zeros(nrow,ncol);
                order = double('C');
                ret = calllib('libcasadi','casadi_sx_matrix_constant',obj.ind,data,nrow,ncol,order);
                assert(ret==0);
        end
        end
        
        % destructor
        function delete(h)
	    ret = calllib('libcasadi','casadi_sx_matrix_delete',h.ind);
	    assert(ret==0);
        end
        % display
        function disp(h)
	      flag = calllib('libcasadi','casadi_sx_matrix_print_cout',h.ind);
              assert(flag==0);
        end

        % convert arguments, if necessary
        function [varargout] = convert_arg(varargin)
            assert(nargin == nargout);
            for i=1:nargin
                a = varargin{i};
                if(isa(a,'sx_matrix'))
                    varargout{i} = a;
                else
                    varargout{i} = sx_matrix(a);                    
                end
            end
        end

        % unary operation
        function res = unary(op,a)
	    res = sx_matrix();
	    flag = calllib('libcasadi','casadi_sx_matrix_unary',res.ind,op,a.ind);
	    assert(flag==0);
        end

        % binary operation
        function res = binary(op,a,b)
            [a,b] = convert_arg(a,b);
	    res = sx_matrix();
	    flag = calllib('libcasadi','casadi_sx_matrix_binary',res.ind,op,a.ind,b.ind);
	    assert(flag==0);
        end

        % addition
        function res = plus(a,b)            
            res = binary(sx_matrix.ADD_NODE,a,b);
        end
        
        % subtraction
        function res = minus(a,b)            
            res = binary(sx_matrix.SUB_NODE,a,b);
        end
        
        % elementwise multiplication
        function res = times(a,b)            
           res = binary(sx_matrix.MUL_NODE,a,b);
        end
        
        % matrix multiplication
        function res = mtimes(a,b) 
            if(length(a)==1 || length(b)==1)
                res = binary(sx_matrix.MUL_NODE,a,b);
            else
                [a,b] = convert_arg(a,b);
                res = sx_matrix();
                flag = calllib('libcasadi','casadi_sx_matrix_prod',res.ind,a.ind,b.ind);
                assert(flag==0);
            end
        end
        
        % elementwise division from the right
        function res = rdivide(a,b)            
            res = binary(sx_matrix.DIV_NODE,a,b);
        end

        % elementwise division from the left
        function res = ldivide(a,b)            
            res = binary(sx_matrix.DIV_NODE,b,a);
        end
        
        % matrix division from the right
        function res = mrdivide(a,b) 
            if(length(a)==1 || length(b)==1)
                res = binary(sx_matrix.DIV_NODE,a,b);
            else
                error('mrdivide: solve not implemented');
            end
        end
        
       % matrix division from the left
        function res = mldivide(a,b) 
            if(length(a)==1 || length(b)==1)
                res = binary(sx_matrix.DIV_NODE,b,a);
            else
                error('mrdivide: solve not implemented');
            end
        end
 
        % misc elementary functions
        function res = sin(a)            
            res = unary(sx_matrix.SIN_NODE,a);
        end
        function res = cos(a)            
            res = unary(sx_matrix.COS_NODE,a);
        end
        function res = tan(a)            
            res = unary(sx_matrix.TAN_NODE,a);
        end
        function res = asin(a)            
            res = unary(sx_matrix.ASIN_NODE,a);
        end
        function res = acos(a)            
            res = unary(sx_matrix.ACOS_NODE,a);
        end
        function res = atan(a)            
            res = unary(sx_matrix.ATAN_NODE,a);
        end
        function res = sqrt(a)            
            res = unary(sx_matrix.SQRT_NODE,a);
        end
        function res = uminus(a)            
            res = unary(sx_matrix.NEG_NODE,a);
        end
        function res = exp(a)            
            res = unary(sx_matrix.EXP_NODE,a);
        end
        function res = log(a)            
            res = unary(sx_matrix.LOG_NODE,a);
        end
        function res = erf(a)            
            res = unary(sx_matrix.ERF_NODE,a);
        end

        % elementwise power: .^
        function res = power(a,b)            
            res = binary(sx_matrix.MUL_NODE,a,b);
        end
        
        % minumum of two elements
        function res = min(a,b)            
            res = binary(sx_matrix.FMIN_NODE,a,b);
        end
        
        % maximum of two elements
        function res = max(a,b)            
            res = binary(sx_matrix.FMAX_NODE,a,b);
        end
        
        % matrix power: ^
        function res = mpower(a,b) 
            if(length(a)==1 || length(b)==1)
                res = binary(sx_matrix.MUL_NODE,a,b);
            else
                error('matrix power not implemented');
            end
        end
        
        % Horizontal concatenation
        function res = horzcat(varargin)
            v = sx_matrix_vec(varargin);
            res = sx_matrix();
            flag = calllib('libcasadi','casadi_sx_matrix_horzcat',res.ind,v.ind);
            assert(flag==0);
        end
        
        % Vertical concatenation
        function res = vertcat(varargin)
            v = sx_matrix_vec(varargin);
            res = sx_matrix();
            flag = calllib('libcasadi','casadi_sx_matrix_vertcat',res.ind,v.ind);
            assert(flag==0);
        end
        
        % number of elements
        function l = length(h)
            sz = libpointer('int32Ptr', 0);
            flag = calllib('libcasadi','casadi_sx_matrix_size',h.ind,sz);
            assert(flag==0);
            l = sz.Value;
        end
        
        % transpose
        function res = ctranspose(a)
            res = sx_matrix();
            flag = calllib('libcasadi','casadi_sx_matrix_transpose',res.ind,a.ind);
            assert(flag==0);
        end
        function res = transpose(a)
            res = sx_matrix();
            flag = calllib('libcasadi','casadi_sx_matrix_transpose',res.ind,a.ind);
            assert(flag==0);
        end
    end
end