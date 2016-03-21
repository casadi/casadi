classdef MySVD < casadi.Callback
  properties
    n
    m
    fd
    k
    fwd
    adj
  end
  methods
    function self = MySVD(name, n, m, fd)
      self@casadi.Callback();
      self.n = n;
      self.m = m;
      self.fd = fd;
      self.k = min(n,m);
      construct(self, name);
    end
    function out = get_n_out(self)
      out = 3;
    end
    function out = get_sparsity_in(self, i)
      out = casadi.Sparsity.dense(self.n,self.m);
    end
    function out = get_sparsity_out(self, i)
    switch i
      case 0
        out = casadi.Sparsity.dense(self.n,self.k);
      case 1
        out = casadi.Sparsity.dense(self.k);
      case 2
        out = casadi.Sparsity.dense(self.k, self.m);
      otherwise
        error('No such index')
      end
    end
    function [res] = eval(self, arg)
      [u,s,v] = svd(full(arg{1}));

      res{1} = u;
      res{2} = diag(s);
      res{3} = v;
    end
  end
end
