%
%     This file is part of CasADi.
% 
%     CasADi -- A symbolic framework for dynamic optimization.
%     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
% 
%     CasADi is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 3 of the License, or (at your option) any later version.
% 
%     CasADi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public
%     License along with CasADi; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
% 
% 
1; % This is a script file

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

global op_any_div_scalar = @(a,b) a.__el_div__(b)

global op_scalar_add_any = @(a,b) b.__radd__(a)
global op_scalar_sub_any = @(a,b) b.__rsub__(a)
global op_scalar_mul_any = @(a,b) b.__rmul__(a)
global op_scalar_el_mul_any = @(a,b) b.__rel_mul__(a)
global op_scalar_el_div_any = @(a,b) b.__rel_div__(a)
global op_scalar_pow_any = @(a,b) b.__rpow__(a)
global op_scalar_el_pow_any = @(a,b) b.__rel_pow__(a)

global op_matrix_add_any = @(a,b) b.__radd__(a)
global op_matrix_sub_any = @(a,b) b.__rsub__(a)
global op_matrix_mul_any = @(a,b) b.__rmul__(a)
global op_matrix_el_mul_any = @(a,b) b.__rel_mul__(a)
global op_matrix_el_div_any = @(a,b) b.__rel_div__(a)
global op_matrix_pow_any = @(a,b) b.__rpow__(a)
global op_matrix_el_pow_any = @(a,b) b.__rel_pow__(a)
