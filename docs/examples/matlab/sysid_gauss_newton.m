%
%     MIT No Attribution
%
%     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
%
%     Permission is hereby granted, free of charge, to any person obtaining a copy of this
%     software and associated documentation files (the "Software"), to deal in the Software
%     without restriction, including without limitation the rights to use, copy, modify,
%     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
%     permit persons to whom the Software is furnished to do so.
%
%     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
%     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
%     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
%     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
%     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
%     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%
%

function solver = sysid_gauss_newton(e,nlp,V)
  % Helper funciton for sysid.m
  J = jacobian(e,V);
  H = triu(J'*J);
  sigma = casadi.MX.sym('sigma');
  
  io = struct;
  io.x = V;
  io.lam_f = sigma;
  io.hess_gamma_x_x = sigma*H;
  
  opts = struct;
  opts.jit = true;
  opts.compiler='shell';
  opts.jit_options.verbose = true;
  hessLag = casadi.Function('nlp_hess_l',io,{'x','p','lam_f','lam_g'}, {'hess_gamma_x_x'},opts);
  opts.hess_lag = hessLag;
  solver = casadi.nlpsol('solver','ipopt', nlp, opts);
end
