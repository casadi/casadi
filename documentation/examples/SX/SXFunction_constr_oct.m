%! Load casadi
casadi

x = SX("x"); % A scalar symbolic
y = ssym("y",2,1); % A matrix symbolic


in = {x y}; % function inputs
disp('Function outputs are:')
out = {x,y,[x x; x x],y*x,0}

f = SXFunction(in,out);
f.init();

%! f now has two inputs and a 4 outputs:
number_in = f.getNumInputs()
number_out = f.getNumOutputs()

%! The outputs has the following string representation.
%! Note how all elements of out have been converted to SXMatrix by
%! automatic typecasting functionality

for i = 0:3
  disp(["Output:"])
  disp(f.outputSX(i))
end
