classdef mycallback < casadi.Callback2

    methods
        function argout = paren(self,argin)
            argout{1} = argin{1}^2;
        end
    end
    
end

