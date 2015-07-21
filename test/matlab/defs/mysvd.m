classdef mysvd < casadi.Callback2

    properties
      n
      m
      fd
      k
      fwd
      adj
    end
    methods
        function self = mysvd(n,m,fd)
          self.n = n
          self.m = m
          self.fd = fd
          self.k = min(n,m)
          self.fwd = mydergensvd(true)
          self.adj = mydergensvd(false)
        end
        function out = nOut(self)
          out = 3;
        end
        function out = inputShape(self,i)
          out = [self.n,self.m];
        end
        function out = outputShape(self,i)
          if i==0
            out = [self.n, self.k];
          elseif i==1
            out = [self.k, 1];
          else
            out = [self.k, self.m];
          end
        end
        function [argout] = paren(self,argin)
            [u,s,v] = svd(full(argin{1}));
            
            argout{1} = u;
            argout{2} = diag(s);
            argout{3} = v;
        end
        function out = options(self)
          if self.fd
            out = struct();
          else
            out = struct('custom_forward',self.fwd.create(), 'custom_reverse',self.adj.create());
          end
        end
    end
    
end

