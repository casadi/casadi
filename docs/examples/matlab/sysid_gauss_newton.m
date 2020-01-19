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
