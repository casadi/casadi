classdef Convolution < casadi.Callback
    properties
      A
      n
      m
      data
      flag_transp
    end

    methods
        function self = Convolution(A, flag_transp)
          self@casadi.Callback();
          self.A = A;
          if flag_transp
            [self.m,self.n] = size(A);
          else
            [self.n,self.m] = size(A);
          end
          self.data = {};
          self.flag_transp = flag_transp;
          construct(self,'Convolution');
        end
        function out = get_transpose(self)
          out = Convolution(self.A,~self.flag_transp);
          self.data = {self.data{:} out};
        end
        function out = get_sparsity_in(self,i)
          out = casadi.Sparsity.dense(self.m,1);
        end
        function out = get_sparsity_out(self,i)
          out = casadi.Sparsity.dense(self.n,1);
        end

        function out = eval(self, ins)
            x = ins{1};
            % Fill in smarter numerical code

            if self.flag_transp
                y = self.A'*x;
            else
                y = self.A*x;
            end

            out = {y};
        end
        function out = get_n_in(self)
            out = 1;
        end
        function out = get_n_out(self)
            out = 1;
        end
        function out = has_forward(self,nfwd)
            out = true;
        end
        function out = has_reverse(self,nadj)
            out = true;
        end
        function out = get_forward(self,nfwd,name,inames,onames,opts)
            mf = self.map(nfwd, 'serial');
            fwd_in = casadi.MX.sym('x',self.m,nfwd);
            fwd_out = mf(fwd_in);
            arg = casadi.MX.sym('x',self.m,1);
            dummy = casadi.MX.sym('x',casadi.Sparsity(self.n,1));
            out = casadi.Function(name, {arg, dummy, fwd_in}, {fwd_out},...
                                  inames, onames, opts);
        end
        function out = get_reverse(self,nadj,name,inames,onames,opts)
            mf = self.get_transpose().map(nadj, 'serial');
            adj_in = casadi.MX.sym('x',self.n,nadj);
            adj_out = mf(adj_in);
            arg = casadi.MX.sym('x',self.m,1);
            dummy = casadi.MX(self.n,1);
            out = casadi.Function(name, {arg, dummy, adj_in}, {adj_out},...
                                  inames, onames, opts);
        end
        function out = has_jacobian(self)
          out = false;
        end
    end
end
