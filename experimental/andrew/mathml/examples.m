% MATLAB script to generate some mathml files using matlab symbolic toolbox
syms x y z
f = cos(2*x+y)+sin(x);
s = feval(symengine,'generate::MathML',f);
s = char(s);
fid=fopen('example1.xml','w');
fprintf(fid,'%s',s);
fclose(fid);
