casadi_interface
global casadi = casadi_interface;

function [y] = casadi_sqrt(x)
  y = x.sqrt();
end
dispatch("sqrt","casadi_sqrt","swig_ref")
function [y] = casadi_sin(x)
  y = x.sin();
end
dispatch("sin","casadi_sin","swig_ref")
function [y] = casadi_cos(x)
  y = x.cos();
end
dispatch("cos","casadi_cos","swig_ref")
function [y] = casadi_tan(x)
  y = x.tan();
end
dispatch("tan","casadi_tan","swig_ref")
function [y] = casadi_atan(x)
  y = x.arctan();
end
dispatch("atan","casadi_atan","swig_ref")
function [y] = casadi_asin(x)
  y = x.arcsin();
end
dispatch("asin","casadi_asin","swig_ref")
function [y] = casadi_acos(x)
  y = x.arccos();
end
dispatch("acos","casadi_acos","swig_ref")
function [y] = casadi_exp(x)
  y = x.exp();
end
dispatch("exp","casadi_exp","swig_ref")
function [y] = casadi_log(x)
  y = x.log();
end
dispatch("log","casadi_log","swig_ref")
function [y] = casadi_floor(x)
  y = x.floor();
end
dispatch("floor","casadi_floor","swig_ref")
function [y] = casadi_ceil(x)
  y = x.ceil();
end
dispatch("ceil","casadi_ceil","swig_ref")
function [y] = casadi_abs(x)
  y = x.fabs();
end
dispatch("abs","casadi_abs","swig_ref")
function [y] = casadi_erf(x)
  y = x.erf();
end
dispatch("erf","casadi_abs","swig_ref")

function [y] = casadi_sparse(x)
  y = x.toSparse();
end
dispatch("sparse","casadi_sparse","swig_ref")

function [y] = casadi_full(x)
  y = full(x.toSparse());
end
dispatch("full","casadi_full","swig_ref")

function [y]=casadi_inv(x)
  global casadi;
  y = casadi.inv(x);
end

dispatch('inv','casadi_inv','swig_ref')

