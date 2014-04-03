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
casadi_main
casadi_integration
casadi_control
casadi_convex_programming
casadi_nonlinear_programming
casadi_optimal_control

failedimport = "Not implemented"

try
casadi_ipopt_interface
catch
casadi_ipopt_interface = failedimport
end_try_catch

try
casadi_sundials_interface
catch
casadi_sundials_interface = failedimport
end_try_catch

try
casadi_qpoases_interface
catch
casadi_qpoases_interface = failedimport
end_try_catch

try
casadi_dsdp_interface
catch
casadi_dsdp_interface = failedimport
end_try_catch

try
casadi_csparse_interface
catch
casadi_csparse_interface = failedimport
end_try_catch

try
casadi_knitro_interface
catch
casadi_knitro_interface = failedimport
end_try_catch

try
casadi_cplex_interface
catch
casadi_cplex_interface = failedimport
end_try_catch

try
casadi_ooqp_interface
catch
casadi_ooqp_interface = failedimport
end_try_catch

try
casadi_slicot_interface
catch
casadi_slicot_interface = failedimport
end_try_catch

try
casadi_worhp_interface
catch
casadi_worhp_interface = failedimport
end_try_catch

try
casadi_snopt_interface
catch
casadi_snopt_interface = failedimport
end_try_catch

try
casadi_lapack_interface
catch
casadi_lapack_interface = failedimport
end_try_catch

global casadi = struct();
names = {'casadi_main','casadi_integration','casadi_control', 'casadi_convex_programming','casadi_nonlinear_programming','casadi_optimal_control','casadi_ipopt_interface','casadi_sundials_interface','casadi_qpoases_interface','casadi_dsdp_interface','casadi_csparse_interface','casadi_knitro_interface','casadi_cplex_interface','casadi_ooqp_interface','casadi_slicot_interface','casadi_worhp_interface','casadi_snopt_interface','casadi_lapack_interface'};
interfaces = { casadi_main casadi_integration casadi_control casadi_convex_programming casadi_nonlinear_programming casadi_optimal_control casadi_ipopt_interface casadi_sundials_interface casadi_qpoases_interface casadi_dsdp_interface casadi_csparse_interface casadi_knitro_interface casadi_cplex_interface casadi_ooqp_interface casadi_slicot_interface casadi_worhp_interface casadi_snopt_interface casadi_lapack_interface};

for i=1:numel(names)
  name = names{i};
  s=completion_matches([name "."]);
  strlen = numel(name);
  for j=1:size(s,1)
    command = deblank(s(j,strlen+2:end));
    if ~isfield(casadi,command)
      casadi.(command) = subsref(interfaces{i},struct('type','.','subs',command));
    end
  end
end


casadi_helpers


