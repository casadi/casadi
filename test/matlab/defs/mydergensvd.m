classdef mydergensvd < casadi.DerivativeGenerator2
  properties
    fwd
  end
  methods
    function self = mydergensvd(fwd)
      self.fwd = fwd;
    end
    function out = paren(self,fcn,ndir)
      import casadi.*
      % Obtain the symbols for nominal inputs/outputs
      nominal_in  = fcn.symbolicInput();
      nominal_out = fcn.symbolicOutput();
      % A growing list of inputs to the returned derivative function
      der_ins = {nominal_in{:}, nominal_out{:}};

      % A growing list of outputs to the returned derivative function
      der_outs = {};

      [A] = deal(nominal_in{:});
      [U,S,V] = deal(nominal_out{:});

      orth_U = (U'*U-eye(U.size2()));
      orth_V = (V*V'-eye(V.size1()));
      
      orth_Ut = orth_U.tril();
      orth_Vt = orth_V.tril();
      
      constr = [
        vec(U*diag(S)*V'-A);
        vec(orth_Ut.getNZ(false,0:orth_Ut.nnz()-1));
        vec(orth_Vt.getNZ(false,0:orth_Vt.nnz()-1))
      ];

      nominal_out_vec = {};
      for i=1:length(nominal_out)
         nominal_out_vec = {nominal_out_vec{:},vec(nominal_out{i})};
      end
      USV = vertcat(nominal_out_vec{:});

      f = MXFunction('f',{USV,A},{constr});

      impl = ImplicitFunction('impl','nlp.ipopt',f,struct('nlp_solver_options',struct('print_level',0,'print_time',false)));
      
      if self.fwd
        fd = impl.derivative(ndir,0);

        seeds = {};
        for i=1:ndir
          symin = fcn.symbolicInput();
          seeds = {seeds{:},symin{1}};
        end

        der_ins={der_ins{:},seeds{:}};

        allseeds = {};
        for s= seeds
          s = s{1};
          allseeds = {allseeds{:},0,s};
        end
        out = fd({USV,A,allseeds{:}});

        for s={out{2:end}}
          s = s{1};
          sp = vertsplit(s,{0,U.nnz(),U.nnz()+S.nnz(),U.nnz()+S.nnz()+V.nnz()});
          [du,ds,dv]=deal(sp{:});

          du = du.reshape(U.size1(),U.size2());
          ds = ds.reshape(S.size1(),S.size2());
          dv = dv.reshape(V.size1(),V.size2());
          der_outs = {der_outs{:},du,ds,dv};
        end

      else
        bd = impl.derivative(0,ndir);
        seeds = {};
        seedsflat = {};
        for i=1:ndir
          symout = fcn.symbolicOutput();
          seeds = {seeds{:},symout};
          der_ins = {der_ins{:},symout{:}};
          symout_flat = {};
          for s=symout
              symout_flat={symout_flat{:},vec(s{1})};
          end
          seedsflat = {seedsflat{:},vertcat(symout_flat{:})};
        end

        out = bd({USV,A,seedsflat{:}});
        
        der_outs = {out{3:2:end}};
      end
      out = MXFunction('my_derivative', der_ins, der_outs);
    end
  end
end
  
