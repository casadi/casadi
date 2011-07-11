# This file makes pylab act inert and without leaks

def nullfun(*args,**kwargs):
  pass
  
plot = nullfun
xlabel = nullfun
ylabel = nullfun
show = nullfun
plot = nullfun
xlabel = nullfun
ylabel = nullfun
show = nullfun
figure = nullfun
legend = nullfun
spy = nullfun
title = nullfun
grid = nullfun
loglog = nullfun
