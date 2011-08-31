classdef swig_ref < handle
    properties (SetAccess = protected)
        ind = libpointer; % a nullpointer
    end
    properties (Constant)
    end
        
    methods
        % constructor
        function h = mx(varargin)
          h.ind = libpointer('int32', 5);
        end

        function convert_input_swig(h)
            libname = 'libcasadi_matlab';
            calllib(libname,'convert_input_swig',h.ind);
        end

        % destructor
        function delete(h)
%            ret = calllib('libcasadi','casadi_mx_delete',h.ind);
%            assert(ret==0);
        end
    end
end
