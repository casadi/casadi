  a = mx('a')
  cos(a) + 3
  f = sin(ans \ log(a*a))

  % create a function
  fcn = mx_function(a,f)

  % enable automatic differentiation
  fcn.set_option('ad_order',1);

  % initialize the function
  fcn.init()

  % set the argument
  fcn.set_input(5)

  % evaluate
  fcn.evaluate()

  % get the output
  out_val = fcn.get_output()

  % set the input seed
  fcn.set_input(1,0,1);

  % evaluate forward
  fcn.evaluate_fwd();

  % get the output seed
  out_seed = fcn.get_output(0,1)

  % set the output seed
  fcn.set_output(1,0,1);

  % evaluate adjoint
  fcn.evaluate_adj();

  % get the input seed
  in_seed = fcn.get_input(0,1)










