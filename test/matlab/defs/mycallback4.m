classdef mycallback4< casadi.Callback2
    properties
      fin_diff_eps
    end
    methods
        function self = mycallback4(fin_diff_eps)
          self.fin_diff_eps = fin_diff_eps
        end
        function argout = paren(self,argin)
            argout{1} = argin{1}^2;
        end
        function out = options(self)
            out = struct('fin_diff_eps',self.fin_diff_eps)
        end
    end
    
end

