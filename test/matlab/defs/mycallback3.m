classdef mycallback3 < casadi.Callback2

    methods
        function argout = paren(self,argin)
            argout{1} = argin{1}^2;
        end
        function out = options(self)
            out = struct('foobar',12)
        end
    end
    
end

